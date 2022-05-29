(* ::Package:: *)

BeginPackage["pkg`"];



EndPackage[];

If[$VersionNumber < 10.,
  Print["This package requires Mathematica version 10 or greater"];
  Abort[]];


Get["pkg`freq`"];
Get["pkg`Jsqfl`"];
Get["pkg`Lsqfl`"];
Get["pkg`Jzfl`"];
Get["pkg`SeffLfl`"];
