(* ::Package:: *)

BeginPackage[ "BBHpnToolkit`Lflnum`"]

     NmLflow::usage = 
	"NmLflow implements flow along L in phase space numerically"

  Begin[ "`Private`"]

    NmLflow[G_,m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]max_,\[Epsilon]_]:=
    Module[{c,Linit,LN,Jinit,zhat,R,P,S1,S2,eqa,eqb,eqc,eqd,system0,initCond,
            sol,\[Lambda],finalvec,Rx,Ry,Rz,Px,Py,Pz,S1x,S1y,S1z,S2x,S2y,S2z,\[Lambda]0},

                      c = 1/Sqrt[\[Epsilon]]   ;   \[Lambda]0=0; (*set initial time to 0 such that \[Lambda]max is the flow amount*) 

                 Linit=Cross[Rinit, Pinit];
                 Jinit=Linit+S1init+S2init; (* J  stays conserved *)
                 LN=Norm[Linit];


              R={Rx[\[Lambda]],Ry[\[Lambda]],Rz[\[Lambda]]};P={Px[\[Lambda]],Py[\[Lambda]],Pz[\[Lambda]]};
              S1={S1x[\[Lambda]],S1y[\[Lambda]],S1z[\[Lambda]]};S2={S2x[\[Lambda]],S2y[\[Lambda]],S2z[\[Lambda]]};
              L=Cross[R, P];


              eqa=D[R,\[Lambda]]-1/(2 LN)( Cross[2 (*L*) Linit,R]);
              eqb=D[P,\[Lambda]]-1/(2 LN)( Cross[2 (*L*) Linit,P]);
              eqc=D[S1,\[Lambda]];
              eqd=D[S2,\[Lambda]];


              system0={eqa[[1]]==0,eqa[[2]]==0,eqa[[3]]==0, eqb[[1]]==0,eqb[[2]]==0,eqb[[3]]==0,
              eqc[[1]]==0,eqc[[2]]==0,eqc[[3]]==0,eqd[[1]]==0,eqd[[2]]==0,eqd[[3]]==0}//Simplify;
              initCond = {Rx[\[Lambda]0] == Rinit[[1]], Ry[\[Lambda]0] == Rinit[[2]], Rz[\[Lambda]0] == Rinit[[3]], 
              Px[\[Lambda]0] == Pinit[[1]], Py[\[Lambda]0] == Pinit[[2]], Pz[\[Lambda]0] == Pinit[[3]] , 
              S1x[\[Lambda]0] ==S1init[[1]], S1y[\[Lambda]0] ==S1init[[2]], S1z[\[Lambda]0] == S1init[[3]], 
              S2x[\[Lambda]0] == S2init[[1]], S2y[\[Lambda]0] ==S2init[[2]], S2z[\[Lambda]0] == S2init[[3]]}  ;
              
              sol=NDSolve[   system0~Join~initCond  ,  {Rx, Ry, Rz, Px, Py,Pz,S1x, S1y, S1z,S2x, S2y, S2z},{\[Lambda],\[Lambda]0,\[Lambda]max}][[1]];

              finalvec=Re[{{Rx[\[Lambda]],Ry[\[Lambda]],Rz[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},{Px[\[Lambda]],Py[\[Lambda]],Pz[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},
              {S1x[\[Lambda]],S1y[\[Lambda]],S1z[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},{S2x[\[Lambda]],S2y[\[Lambda]],S2z[\[Lambda]]}/. sol/.{\[Lambda]->\[Lambda]max}}]//N;
              
              Return[finalvec];
               

]
  End[]

  EndPackage[]

