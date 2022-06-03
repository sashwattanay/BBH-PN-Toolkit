(* ::Package:: *)

BeginPackage["pkg`"];



EndPackage[];

If[$VersionNumber < 10.,
  Print["This package requires Mathematica version 10 or greater"];
  Abort[]];


Get["pkg`freq`"];
Get["pkg`Jfl`"];
Get["pkg`Lfl`"];
Get["pkg`Jzfl`"];
Get["pkg`SeffLfl`"];
Get["pkg`Hfl`"];
Get["pkg`Jflnum`"];
Get["pkg`Lflnum`"];
Get["pkg`Jzflnum`"];
Get["pkg`SeffLflnum`"];
Get["pkg`Hflnum`"];
