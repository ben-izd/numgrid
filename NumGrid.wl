BeginPackage["NumGrid`"];

ClearAll[NumGrid];

Begin["`Private`"];

ClearAll[Unload];

Unload[]:=(LibraryUnload[NumGrid`$libraryPath];LibraryUnload[NumGrid`$libraryPath];)

If[FileExistsQ[NumGrid`$libraryPath],
$Types={"A","B","C","D","H","L","S"};

Block[{library=LibraryLoad[NumGrid`$libraryPath]},

{NumGridA,NumGridB,NumGridC,NumGridD,NumGridH,NumGridL,NumGridS}=LibraryFunctionLoad[library,"numgrid_"<>#,{Integer},LibraryDataType["NumericArray","UnsignedInteger32",2]]&/@ToLowerCase[$Types];
ParallelNumGridB=LibraryFunctionLoad[library,"numgrid_b_parallel",{Integer},LibraryDataType["NumericArray","UnsignedInteger32",2]];
]

Options[NumGrid]={Parallelization->False};

NumGrid::invalidType="Type `1` is not valid. Valid types are `2` .";
NumGrid[type_String,n_Integer,OptionsPattern[]]:=Switch[ToUpperCase@type,
"S",NumGridS[n],
"L",NumGridL[n],
"C",NumGridC[n],
"D",NumGridD[n],
"A",NumGridA[n],
"H",NumGridH[n],
"B",If[TrueQ[OptionValue[Parallelization]],ParallelNumGridB[n],NumGridB[n]],
_,Message[NumGrid::invalidType,type,StringRiffle[$Types,", "]];$Failed]
,Print["Variable \"NumGrid`$libraryPath\" is not set. Can not initialize."];$Failed]

End[];
EndPackage[];
