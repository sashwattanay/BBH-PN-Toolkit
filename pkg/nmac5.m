(* ::Package:: *)

BeginPackage["pkg`nmac5`", 
  {"pkg`freq`", "pkg`Jflnum`","pkg`Jzflnum`","pkg`Lflnum`","pkg`Hflnum`","pkg`SeffLflnum`"}];
  NmJac5::usage = 
	"NmJac5 implements flow along first action J5 in phase space numerically"

              Begin[ "`Private`"]

               NmJac5[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]mx_,\[Epsilon]_]:=
               Module[{fvec,Jout,Jzout,Lout,Hout,Seffout,fr, jcm},
                
                       fr= frequency[m1,m2,Rinit,Pinit,S1init,S2init,\[Epsilon]];
                       jcm=fr[[3]];

                Jout=NmJflow[m1,m2,Rinit,Pinit,S1init,S2init,\[Lambda]mx *jcm[[5]][[1]] ,\[Epsilon]];
               Jzout=NmJzflow[m1,m2,Jout[[1]],Jout[[2]],Jout[[3]],Jout[[4]], \[Lambda]mx*jcm[[5]][[2]],\[Epsilon]];
                Lout=NmLflow[m1,m2,Jzout[[1]],Jzout[[2]],Jzout[[3]],Jzout[[4]],\[Lambda]mx *jcm[[5]][[3]],\[Epsilon]];
                Hout=NmHflow[m1,m2,Lout[[1]],Lout[[2]],Lout[[3]],Lout[[4]],\[Lambda]mx *jcm[[5]][[4]],\[Epsilon]];
             Seffout=NmSeffLflow[m1,m2,Hout[[1]],Hout[[2]],Hout[[3]],Hout[[4]],\[Lambda]mx *jcm[[5]][[5]],\[Epsilon]];
                
                Return[Seffout];
                
                 ]
                 End[]

                 EndPackage[]
