(* ::Package:: *)

BeginPackage["BBHpnToolkit`action1`", 
  {"BBHpnToolkit`Jfl`"}];
  Jac1::usage = 
	"Jac1 implements flow along first action J1=J in phase space"

              Begin[ "`Private`"]

               Jac1[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{fvec},
                

                fvec=Re[Jflow[m1,m2,Rinit,Pinit,S1init,S2init,\[Lambda]mx,\[Epsilon]]]//N;
                
                Return[fvec];
                
                 ]
                 End[]

                 EndPackage[]
