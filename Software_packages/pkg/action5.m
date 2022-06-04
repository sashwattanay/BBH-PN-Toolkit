(* ::Package:: *)

BeginPackage["pkg`action5`", 
  {"pkg`freq`", "pkg`Jfl`","pkg`Jzfl`","pkg`Lfl`","pkg`Hfl`","pkg`SeffLfl`"}];
  Jac5::usage = 
	"Jac5 implements flow along fifth action J5 in phase space"

              Begin[ "`Private`"]

               Jac5[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{Jout,Jzout,Lout,Hout,Seffout,fr, jcm},
                
                       fr= frequency[m1,m2,Rinit,Pinit,S1init,S2init,\[Epsilon]];
                       jcm=fr[[3]];

                Jout=Jflow[m1,m2,Rinit,Pinit,S1init,S2init,\[Lambda]mx *jcm[[5]][[1]] ,\[Epsilon]];
               Jzout=Jzflow[m1,m2,Jout[[1]],Jout[[2]],Jout[[3]],Jout[[4]], \[Lambda]mx*jcm[[5]][[2]],\[Epsilon]];
                Lout=Lflow[m1,m2,Jzout[[1]],Jzout[[2]],Jzout[[3]],Jzout[[4]],\[Lambda]mx *jcm[[5]][[3]],\[Epsilon]];
                Hout=Hflow[m1,m2,Lout[[1]],Lout[[2]],Lout[[3]],Lout[[4]],\[Lambda]mx *jcm[[5]][[4]],\[Epsilon]];
             Seffout=SeffLflow[m1,m2,Hout[[1]],Hout[[2]],Hout[[3]],Hout[[4]],\[Lambda]mx *jcm[[5]][[5]],\[Epsilon]];
                
                Return[Seffout];
                
                 ]
                 End[]

                 EndPackage[]
