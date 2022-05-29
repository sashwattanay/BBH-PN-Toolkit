(* ::Package:: *)

BeginPackage[ "pkg`Lsqfl`"]

     Lsqflow::usage = 
	"Lsqflow implements flow along Lsquared in phase space"

  Begin[ "`Private`"]

    Lsqflow[m1N_, m2N_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]max_,\[Lambda]0_]:=
    Module[{Linit,Jinit,zhat,R,P,S1,S2,eqa,eqb,eqc,eqd,system0,initCond,
            sol,\[Lambda],finalvec},

                  

                 Linit=Cross[Rinit, Pinit];
                 Jinit=Linit+S1init+S2init; (* J  stays conserved *)


              R={Rx[\[Lambda]],Ry[\[Lambda]],Rz[\[Lambda]]};P={Px[\[Lambda]],Py[\[Lambda]],Pz[\[Lambda]]};
              S1={S1x[\[Lambda]],S1y[\[Lambda]],S1z[\[Lambda]]};S2={S2x[\[Lambda]],S2y[\[Lambda]],S2z[\[Lambda]]};



              eqa=D[R,\[Lambda]]-( Cross[2 Linit,R]);
              eqb=D[P,\[Lambda]]-( Cross[2 Linit,P]);
              eqc=D[S1,\[Lambda]];
              eqd=D[S2,\[Lambda]];


              system0={eqa[[1]]==0,eqa[[2]]==0,eqa[[3]]==0, eqb[[1]]==0,eqb[[2]]==0,eqb[[3]]==0,
              eqc[[1]]==0,eqc[[2]]==0,eqc[[3]]==0,eqd[[1]]==0,eqd[[2]]==0,eqd[[3]]==0}//Simplify;
              initCond = {Rx[\[Lambda]0] == Rinit[[1]], Ry[\[Lambda]0] == Rinit[[2]], Rz[\[Lambda]0] == Rinit[[3]], 
              Px[\[Lambda]0] == Pinit[[1]], Py[\[Lambda]0] == Pinit[[2]], Pz[\[Lambda]0] == Pinit[[3]] , 
              S1x[\[Lambda]0] ==S1init[[1]], S1y[\[Lambda]0] ==S1init[[2]], S1z[\[Lambda]0] == S1init[[3]], 
              S2x[\[Lambda]0] == S2init[[1]], S2y[\[Lambda]0] ==S2init[[2]], S2z[\[Lambda]0] == S2init[[3]]}  ;
              
              sol=DSolve[   system0~Join~initCond  ,  {Rx, Ry, Rz, Px, Py,Pz,S1x, S1y, S1z,S2x, S2y, S2z},\[Lambda]][[1]];

              finalvec=Re[{{Rx[\[Lambda]],Ry[\[Lambda]],Rz[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},{Px[\[Lambda]],Py[\[Lambda]],Pz[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},
              {S1x[\[Lambda]],S1y[\[Lambda]],S1z[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},{S2x[\[Lambda]],S2y[\[Lambda]],S2z[\[Lambda]]}/. sol/.{\[Lambda]->\[Lambda]max}}]//N;
              
              Print[finalvec];


Show[Graphics3D[{{Red,Arrowheads[0.03],Arrow[{{0,0,0}, Linit}]},{Blue,Arrowheads[0.03],Arrow[{{0,0,0},Rinit}]},
{Green,Arrowheads[0.03],Arrow[{{0,0,0},Pinit}]},{Brown,Arrowheads[0.03],Arrow[{{0,0,0},S1init}]},
{Magenta,Arrowheads[0.03],Arrow[{{0,0,0},S2init}]}}],ParametricPlot3D[Evaluate[{Rx[\[Lambda]],Ry[\[Lambda]],Rz[\[Lambda]]}/. sol],{\[Lambda],\[Lambda]0,\[Lambda]max},
PlotStyle->Blue],ParametricPlot3D[Evaluate[{Px[\[Lambda]],Py[\[Lambda]],Pz[\[Lambda]]}/. sol],{\[Lambda],\[Lambda]0,\[Lambda]max},PlotStyle->Green],
ParametricPlot3D[Evaluate[{S1x[\[Lambda]],S1y[\[Lambda]],S1z[\[Lambda]]}/. sol],{\[Lambda],\[Lambda]0,\[Lambda]max},PlotStyle->Brown],
ParametricPlot3D[Evaluate[{S2x[\[Lambda]],S2y[\[Lambda]],S2z[\[Lambda]]}/. sol],{\[Lambda],\[Lambda]0,\[Lambda]max},PlotStyle->Magenta],Boxed->False]

]
  End[]

  EndPackage[]

