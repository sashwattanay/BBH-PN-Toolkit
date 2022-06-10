(* ::Package:: *)

BeginPackage["BBHpnToolkit`action2`", 
  {"BBHpnToolkit`Jzfl`"}];
  Jac2::usage = 
	"Jac2 implements flow along second action J2=Jz in phase space"

              Begin[ "`Private`"]

               Jac2[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{fvec},
                

                fvec=Re[Jzflow[m1,m2,Rinit,Pinit,S1init,S2init,\[Lambda]mx,\[Epsilon]]]//N;
                
                Return[fvec];
                
                 ]
                 End[]

                 EndPackage[]
