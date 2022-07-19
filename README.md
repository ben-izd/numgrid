# numgrid
Matlab numgrid implementation in Rust for Mathematica


### Instruction 
1. Download the `.dll` file from [Release Section](https://github.com/ben-izd/numgrid/releases/tag/Main).
2. In Mathematica set ```NumGrid`$libraryPath``` variable to the path you downloaded `.dll` file.
3. Run [`NumGrid.wl`](https://github.com/ben-izd/numgrid/blob/main/NumGrid.wl) file which define a `NumGrid` function that include all the interfaces.


`NumGrid` supports `S`, `L`, `C`, `D`, `A`, `H` and `B` except `N` types. You can `NumGrid` as follows:
```
(* Support lower case *)
NumGrid["a", 5]
  
(* Support upper case *)
NumGrid["A", 5]

(* Only "B" type can run in parallel *)
NumGrid["B", 5, Parallelization -> True]
```

### Non-Windows users
Install [Rust](https://www.rust-lang.org/) and build the project on your OS, then follow from step 2 of instruction.

<hr>

You can find solutions implemeneting `B` type in Mathematica, in [this Mathematica Stackexchange question](https://mathematica.stackexchange.com/q/270516/77079).
