# numgrid
Matlab numgrid implementation in Rust for Mathematica

Matlab `numgrid`, has multiple types, for now, only the 'B' type is implemented.

You can find more solutions implemeneted in Mathematica, in [this Mathematica Stackexchange question](https://mathematica.stackexchange.com/q/270516/77079).

Download the `.dll` file from release section and use the following command in Mathematica:
```
NumGridBCompiled = 
 LibraryFunctionLoad["C:\\numgrid.dll", "numgrid_b", {Integer}, 
  LibraryDataType["NumericArray", "UnsignedInteger32", 2]]
  
NumGridBParallelCompiled = 
 LibraryFunctionLoad["C:\\numgrid.dll", "numgrid_b_parallel", {Integer}, 
  LibraryDataType["NumericArray", "UnsignedInteger32", 2]]
```
`NumGridBCompiled` is single threaded, and `NumGridBParallelCompiled` is multi-thread. Both return  `NumericArray` with "UnsignedInteger32" type.
