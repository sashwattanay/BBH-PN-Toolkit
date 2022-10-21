(* ::Package:: *)

BeginPackage["BBHpnToolkit`AAflow`", 
  {"BBHpnToolkit`freq`", "BBHpnToolkit`action1`","BBHpnToolkit`action2`","BBHpnToolkit`action3`","BBHpnToolkit`action4`","BBHpnToolkit`action5`"}];
  Acflow::usage = 
	"Acflow implements Hamiltonian flow by flowing along the five actions by appropriate flow amounts dictated by frequency"

              Begin[ "`Private`"]

               Acflow[G_,m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{j1out,j2out,j3out,j4out,j5out,fr, frq},
                
                       fr= frequency[G,m1,m2,Rinit,Pinit,S1init,S2init,\[Epsilon]];
                       frq=fr[[2]];
                     j1out=Jac1[G,m1,m2,Rinit,Pinit,S1init,S2init, \[Lambda]mx * frq[[1]] ,\[Epsilon]];
                     j2out=Jac2[G,m1,m2,j1out[[1]],j1out[[2]],j1out[[3]],j1out[[4]], \[Lambda]mx* frq[[2]] ,\[Epsilon]];
                     j3out=Jac3[G,m1,m2,j2out[[1]],j2out[[2]],j2out[[3]],j2out[[4]],\[Lambda]mx * frq[[3]] ,\[Epsilon]];
                     j4out=Jac4[G,m1,m2,j3out[[1]],j3out[[2]],j3out[[3]],j3out[[4]],\[Lambda]mx * frq[[4]] ,\[Epsilon]];
                     j5out=Jac5[G,m1,m2,j4out[[1]],j4out[[2]],j4out[[3]],j4out[[4]],\[Lambda]mx * frq[[5]] ,\[Epsilon]];

                
                Return[j5out];
                
                 ]
                 End[]

                 EndPackage[]
