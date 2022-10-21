(* ::Package:: *)

BeginPackage["BBHpnToolkit`nmac2`", 
  {"BBHpnToolkit`Jzflnum`"}];
  NmJac2::usage = 
	"NmJac2 implements flow along second action J2=Jz in phase space numerically"

              Begin[ "`Private`"]

               NmJac2[G_,m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{fvec},
                

                fvec=Re[NmJzflow[m1,m2,Rinit,Pinit,S1init,S2init,\[Lambda]mx,\[Epsilon]]]//N;
                
                Return[fvec];
                
                 ]
                 End[]

                 EndPackage[]
