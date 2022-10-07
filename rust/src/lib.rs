// ============= Author: Benjamin Izadpanah =============
// This library is based on Matlab implementation of numgrid which implements all types included 'S', 'L', 'C', 'D', 'A', 'H' and 'B' except 'N'.
// Some optimization was done manually to speed up the calculation, like using temporary variables or pre-allocating border with 0, without -1<x<1 and -1<y<1 condition.
// Code description and its license is available at https://github.com/ben-izd/numgrid
// Most of the functions came from `numgrid_generator_float_range` super function except `numgrid_s` to increase performance.
// Only the parallel version of type 'B' was implemented as multi-threading in other types the overhead makes the speed-up useless

use rayon::prelude::*;
use std::mem::MaybeUninit;
use wolfram_library_link::{export, NumericArray, UninitNumericArray};

// This variable set the overall floating point calculation type
// Using f32, might give a little speed up
type Precision = f64;
type OutputArray = NumericArray<u32>;
type UninitArray = UninitNumericArray<u32>;
type Size = i64;

// TYPE 'B' single-threaded
#[export]
fn numgrid_b(n: Size) -> OutputArray {
    let output = generate_padded_array(n);
    numgrid_generator_float_range(n, true, output, &|x: Precision, y: Precision| {
        let t = y.atan2(x);
        (Precision::sqrt(x.powf(2.0) + y.powf(2.0)))
            >= Precision::sin(2.0 * t) + 0.2 * Precision::sin(8.0 * t)
    })
}

// TYPE 'B' multi-threaded
#[export]
fn numgrid_b_parallel(n: Size) -> OutputArray {
    let mut output: UninitArray = UninitNumericArray::from_dimensions(&[n as usize, n as usize]);

    let output_slice = output.as_slice_mut();

    let n_float = n as Precision;
    let mut g = vec![false; (n * n) as usize];
    let n_usize = n as usize;

    // skip first row and use .take to skip the last row
    g.iter_mut()
        .skip(n_usize)
        .take(n_usize * (n_usize - 2))
        .collect::<Vec<&mut bool>>()
        .par_chunks_mut(n as usize)
        .enumerate()
        .for_each(|(y_index, row)| {
            let y = (-(n_float - 3.0) + 2.0 * (y_index as Precision)) / (n_float - 1.0);
            let mut x: Precision;
            for x_index in 1..(n - 1) {
                x = (n_float - 1.0 - 2.0 * (x_index as Precision)) / (n_float - 1.0);
                {
                    let t = y.atan2(x);
                    *row[x_index as usize] = (Precision::sqrt(x.powf(2.0) + y.powf(2.0)))
                        >= Precision::sin(2.0 * t) + 0.2 * Precision::sin(8.0 * t);
                }
            }
        });

    let mut counter: u32 = 1;
    let mut index: usize;

    for x_index in 0..n {
        for y_index in 0..n {
            {
                index = (x_index + y_index * n) as usize;
                if g[index] {
                    output_slice[index].write(counter);
                    counter += 1;
                } else {
                    output_slice[index].write(0);
                }
            }
        }
    }

    unsafe { output.assume_init() }
}

// TYPE 'S'
// since it does not have any condition, for performance reason it will contain duplicate codes
#[export]
fn numgrid_s(n: Size) -> OutputArray {
    
    let mut output: UninitArray = UninitNumericArray::from_dimensions(&[n as usize, n as usize]);

    let output_slice = output.as_slice_mut();

    let (n_usize, mut counter): (usize, u32) = (n as usize, 1);

    set_array_border_zero(n, output_slice);

    for x in 1..(n_usize - 1) {
        for y in 1..(n_usize - 1) {
            output_slice[x + y * n_usize].write(counter);
            counter += 1;
        }
    }

    unsafe { output.assume_init() }
}

// TYPE 'L'
#[export]
fn numgrid_l(n: Size) -> OutputArray {
    numgrid_generator_integer_range_with_padding(n, &|x, y| {

        // n-2 for the padding
        let width = ((n - 2) / 2) as usize;

        // x > (width + 1) || y <= width
        x > width || y <= width
    })
}

// TYPE 'C'
#[export]
fn numgrid_c(n: Size) -> OutputArray {
    let output = generate_padded_array(n);
    numgrid_generator_float_range(n, true, output, &|x: Precision, y: Precision| {
        ((x - 1.0).powf(2.0) + (1.0 - y).powf(2.0)) > 1.0
    })
}

// TYPE 'D'
#[export]
fn numgrid_d(n: Size) -> OutputArray {
    numgrid_generator_float_range_without_padding(n, &|x: Precision, y: Precision| {
        (x.powf(2.0) + (y).powf(2.0)) < 1.0
    })
}

// TYPE 'A'
#[export]
fn numgrid_a(n: Size) -> OutputArray {
    numgrid_generator_float_range_without_padding(n, &|x: Precision, y: Precision| {
        let temp = x.powf(2.0) + y.powf(2.0);
        temp < 1.0 && temp > (1.0 / 3.0)
    })
}

// TYPE 'H'
const RHO: Precision = 0.75;
const SIGMA: Precision = 0.75;

#[export]
fn numgrid_h(n: Size) -> OutputArray {
    numgrid_generator_float_range_without_padding(n, &|x: Precision, y: Precision| {

        let y = -y;
        let x_power_2 = x.powf(2.0);
        let temp = x_power_2 + y.powf(2.0);

        temp * (temp - SIGMA * y) < RHO * x_power_2
    })
}

// bean-like shape - since it's not in original numgrid, it's not included by default
// #[export]
// fn numgrid_h(n: Size) -> OutputArray {
//     numgrid_generator_float_range_without_padding(n, &|x: Precision, y: Precision| {
//         let temp = x.powf(2.0) + y.powf(2.0);
//         temp * temp - SIGMA * y < RHO * x.powf(2.0) - EPSILON
//     })
// }

#[inline(always)]
fn numgrid_generator_float_range_without_padding(
    n: Size,
    predicate: &dyn Fn(Precision, Precision) -> bool,
    ) -> OutputArray {

    let output: UninitArray = UninitNumericArray::from_dimensions(&[n as usize, n as usize]);
    numgrid_generator_float_range(n, false, output, predicate)
}

// A super function to build patterns using float ranges from -1 to 1, step size is dependant on n
// Arguments:
// - n : size of array
// - predicate: function which determine which cell to include with (x:Precision,y:Precision) -> bool signature
#[inline(always)]
fn numgrid_generator_float_range(
    n: Size,
    padded: bool,
    mut output: UninitArray,
    predicate: &dyn Fn(Precision, Precision) -> bool,
    ) -> OutputArray {

    let offset = if padded { 1 } else { 0 };
    let output_slice = output.as_slice_mut();

    let (n_usize, n_float, mut counter): (usize, Precision, u32) = (n as usize, n as Precision, 1);
    let (mut index, mut x, mut y): (usize, Precision, Precision);

    // iterating over grid, set the grid column by column
    for x_index in offset..(n_usize - offset) {

        x = (n_float - 1.0 - 2.0 * (x_index as Precision)) / (n_float - 1.0);

        for y_index in offset..(n_usize - offset) {

            y = (-(n_float - 1.0) + 2.0 * (y_index as Precision)) / (n_float - 1.0);
            index = x_index + y_index * n_usize;

            if predicate(x, y) {
                output_slice[index].write(counter);
                counter += 1;
            } else {
                output_slice[index].write(0);
            }
        }
    }

    unsafe { output.assume_init() }
}

// generate a padded array with size n*n
// effectively set the borders (first row, last row, first column, last column) to zero
#[inline(always)]
fn generate_padded_array(n: Size) -> UninitArray {
    let mut output: UninitArray = UninitNumericArray::from_dimensions(&[n as usize, n as usize]);

    set_array_border_zero(n, output.as_slice_mut());

    output
}

// A super function to build patterns using integer ranges 0..n
// Arguments:
// - n : size of array
// - predicate: function which determine which cell to include with (x:usize,y:usize) -> bool signature
#[inline(always)]
fn numgrid_generator_integer_range_with_padding(
    n: Size,
    predicate: &dyn Fn(usize, usize) -> bool,
    ) -> OutputArray {

    let mut output: UninitArray = generate_padded_array(n);
    let output_slice = output.as_slice_mut();
    let n_usize = n as usize;
    let mut counter: u32 = 1;

    let mut index: usize;

    for x in 1..(n_usize - 1) {
        for y in 1..(n_usize - 1) {
            index = x + y * n_usize;
            if predicate(x, y) {
                output_slice[index].write(counter);
                counter += 1;
            } else {
                output_slice[index].write(0);
            }
        }
    }

    unsafe { output.assume_init() }
}

// setting array border (first row, last row, first column, last column) to zero
fn set_array_border_zero(n: Size, array_slice: &mut [MaybeUninit<u32>]) {

    let n_usize = n as usize;

    for i in 0..n_usize {

        // top border
        array_slice[i].write(0);

        // left border
        array_slice[i * n_usize].write(0);

        // right border
        array_slice[n_usize * (i + 1) - 1].write(0);

        // bottom border
        array_slice[n_usize * (n_usize - 1) + i].write(0);

    }
}
