(* ::Package:: *)

BeginPackage["pkg`action3`", 
  {"pkg`Lfl`"}];
  Jac3::usage = 
	"Jac3 implements flow along third action J3=L in phase space"

              Begin[ "`Private`"]

               Jac3[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{fvec},
                

                fvec=Re[Lflow[m1,m2,Rinit,Pinit,S1init,S2init,\[Lambda]mx,\[Epsilon]]]//N;
                
                Return[fvec];
                
                 ]
                 End[]

                 EndPackage[]
