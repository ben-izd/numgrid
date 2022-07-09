use wolfram_library_link::{export, NumericArray, UninitNumericArray};
use rayon::prelude::*;

const EPSILON:f32=0.00001;
type Precision = f32;

#[export]
fn numgrid_b(n:i64) -> NumericArray<u32>{
    let mut output: UninitNumericArray<u32> =
        UninitNumericArray::from_dimensions(&[n as usize,n as usize]);

    let output_slice = output.as_slice_mut();

    let mut x:Precision;
    let mut y:Precision;
    let mut t:Precision;
    let mut counter:u32 = 1;

    let step:Precision = 2.0 / (n as Precision - 1.0);
    x = -1.0;

    for x_index in 0..n{

        y = 1.0;
        for y_index in 0..n{
            {
                t = y.atan2(x);

                if (x.abs()<1.0-EPSILON) & (y.abs()<1.0-EPSILON) & ((Precision::sqrt(x.powf(2.0) + y.powf(2.0))) >= Precision::sin(2.0*t) + 0.2 * Precision::sin(8.0*t) )
                {
                    output_slice[(x_index + y_index * n) as usize].write(counter);
                    counter += 1;
                }else{
                    output_slice[(x_index + y_index * n) as usize].write(0);
                }

            }
            y -=step;
        }
        x += step;
    }

    unsafe { output.assume_init() }
}


#[export]
fn numgrid_b_parallel(n:i64) -> NumericArray<u32>{
    let mut output: UninitNumericArray<u32> =
        UninitNumericArray::from_dimensions(&[n as usize,n as usize]);

    let output_slice = output.as_slice_mut();

    let step:Precision = 2.0 / (n as Precision - 1.0);

    let mut g = vec![false;(n*n) as usize];

    g.par_chunks_mut(n as usize).enumerate().for_each(|(x_index,row)|{
        let x = -1.0-step + x_index as f32 * step;
        let mut y:f32 = 1.0;
        for y_index in 0..n{
            {
                let t = y.atan2(x);

                row[y_index as usize] = ((x.abs() < 1.0 - EPSILON) & (y.abs() < 1.0 - EPSILON) & ((Precision::sqrt(x.powf(2.0) + y.powf(2.0))) >= Precision::sin(2.0 * t) + 0.2 * Precision::sin(8.0 * t)));

            }
            y -= step;
        }
    });

    let mut counter:u32 = 1;
    let mut index:usize;

    for x_index in 0..n{
        for y_index in 0..n{
            {
                index= (x_index + y_index * n) as usize;
                if g[index]
                {
                    output_slice[index].write(counter);
                    counter += 1;
                }else{
                    output_slice[index].write(0);
                }
            }
        }
    }


    unsafe { output.assume_init() }
}