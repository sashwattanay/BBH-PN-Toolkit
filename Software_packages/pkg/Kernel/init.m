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
Get["pkg`Hfl`"];
Get["pkg`Jsqflnum`"];
Get["pkg`Lsqflnum`"];
Get["pkg`Jzflnum`"];
Get["pkg`SeffLflnum`"];
Get["pkg`Hflnum`"];
