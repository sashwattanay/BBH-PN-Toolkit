(* ::Package:: *)

BeginPackage["BBHpnToolkit`nmac1`", 
  {"BBHpnToolkit`Jflnum`"}];
  NmJac1::usage = 
	"NmJac1 implements flow along first action J1=J in phase space numerically"

              Begin[ "`Private`"]

               NmJac1[G_,m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{fvec},
                

                fvec=Re[NmJflow[G,m1,m2,Rinit,Pinit,S1init,S2init,\[Lambda]mx,\[Epsilon]]]//N;
                
                Return[fvec];
                
                 ]
                 End[]

                 EndPackage[]
