(* ::Package:: *)

BeginPackage["BBHpnToolkit`nmac3`", 
  {"BBHpnToolkit`Lflnum`"}];
  NmJac3::usage = 
	"NmJac3 implements flow along third action J3=L in phase space numerically"

              Begin[ "`Private`"]

               NmJac3[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{fvec},
                

                fvec=Re[NmLflow[m1,m2,Rinit,Pinit,S1init,S2init,\[Lambda]mx,\[Epsilon]]]//N;
                
                Return[fvec];
                
                 ]
                 End[]

                 EndPackage[]
