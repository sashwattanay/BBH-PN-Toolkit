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
Get["pkg`action1`"];
Get["pkg`action2`"];
Get["pkg`action3`"];
Get["pkg`action4`"];
Get["pkg`action5`"];
Get["pkg`nmac1`"];
Get["pkg`nmac2`"];
Get["pkg`nmac3`"];
Get["pkg`nmac4`"];
Get["pkg`nmac5`"];
Get["pkg`AAflow`"];
Get["pkg`AAflownum`"];
