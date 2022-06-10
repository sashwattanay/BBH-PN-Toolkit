(* ::Package:: *)

BeginPackage["BBHpnToolkit`AAflownum`", 
  {"BBHpnToolkit`freq`", "BBHpnToolkit`nmac1`","BBHpnToolkit`nmac2`","BBHpnToolkit`nmac3`","BBHpnToolkit`nmac4`","BBHpnToolkit`nmac5`"}];
  NmAcflow::usage = 
	"NmAcflow numerically implements Hamiltonian flow by flowing along the 
            five actions by appropriate flow amounts dictated by frequency"

              Begin[ "`Private`"]

               NmAcflow[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{j1out,j2out,j3out,j4out,j5out,fr, frq},
                
                       fr= frequency[m1,m2,Rinit,Pinit,S1init,S2init,\[Epsilon]];
                       frq=fr[[2]];
                     j1out=NmJac1[m1,m2,Rinit,Pinit,S1init,S2init, \[Lambda]mx * frq[[1]] ,\[Epsilon]];
                     j2out=NmJac2[m1,m2,j1out[[1]],j1out[[2]],j1out[[3]],j1out[[4]], \[Lambda]mx* frq[[2]] ,\[Epsilon]];
                     j3out=NmJac3[m1,m2,j2out[[1]],j2out[[2]],j2out[[3]],j2out[[4]],\[Lambda]mx * frq[[3]] ,\[Epsilon]];
                     j4out=NmJac4[m1,m2,j3out[[1]],j3out[[2]],j3out[[3]],j3out[[4]],\[Lambda]mx * frq[[4]] ,\[Epsilon]];
                     j5out=NmJac5[m1,m2,j4out[[1]],j4out[[2]],j4out[[3]],j4out[[4]],\[Lambda]mx * frq[[5]] ,\[Epsilon]];

                
                Return[j5out];
                
                 ]
                 End[]

                 EndPackage[]
