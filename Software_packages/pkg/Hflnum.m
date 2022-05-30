(* ::Package:: *)

BeginPackage[ "pkg`Hflnum`"]

     NmHflow::usage = 
	"NmHflow implements flow induced by H in phase space numerically"

  Begin[ "`Private`"]

    NmHflow[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]max_,\[Lambda]0_]:=
    Module[{precisionGoal,Q1,Q2,Linit,R,P,S1,S2,eqa,eqb,eqc,eqd,system0,initCond,
                        sol,\[Lambda],finalvec,Rx,Ry,Rz,Px,Py,Pz,S1x,S1y,S1z,S2x,S2y,S2z,Seff,L,spmg},
                                \[Epsilon]=3/1000;  G=1;

                            precisionGoal=(*Automatic*)   100  ;

                             Q1=2+3 m2/(2 m1);  Q2=2+3 m1/(2 m2);
                              Linit=Cross[Rinit,Pinit];



                             R={Rx[\[Lambda]],Ry[\[Lambda]],Rz[\[Lambda]]}; P={Px[\[Lambda]],Py[\[Lambda]],Pz[\[Lambda]]};
                             S1={S1x[\[Lambda]],S1y[\[Lambda]],S1z[\[Lambda]]}; S2={S2x[\[Lambda]],S2y[\[Lambda]],S2z[\[Lambda]]};
                             Seff= (Q1 S1 + Q2 S2) ;
                             L = Cross[R,P] ;



                eqa=D[R,\[Lambda]]-((1/m1 +1/m2)P - ( G \[Epsilon] (R \[Cross] Seff))/(R . R)^(3/2))  + 1/(2  m1^3  m2^3 (R . R)^(3/2)) \[Epsilon](2 G m1^3 m2^3 R  (P . R)+P (2 G  m1^2  m2^2 (3 m1^2+7 m1 m2+3 m2^2)+(m1^3+m2^3) (P . P) Sqrt[R . R]) (R . R)) ;     
                eqb=D[P,\[Lambda]] +1/(R . R)^(5/2) G (R(m1 m2 (R . R)+3 \[Epsilon]( - L . Seff))+\[Epsilon]  (R . R)  Cross[P,Seff])   -  1/(2  m1 m2 (R . R)^(5/2)) G \[Epsilon] (2  m1 m2  P (P . R)   (R . R)+R (-3 m1 m2 (P . R)^2+2 G  m1^2  m2^2 ( m1 + m2) Sqrt[R . R]-(3 m1^2+7 m1 m2+3  m2^2) (P . P)  (R . R))) ;
                eqc=D[S1,\[Lambda]]- (\[Epsilon]  G Q1 (-(P . S1) R +(R . S1) P))/(R . R)^(3/2);
                eqd=D[S2,\[Lambda]]-(\[Epsilon] G Q2 (-(P . S2) R +(R . S2) P))/(R . R)^(3/2);


             system0={eqa[[1]]==0,eqa[[2]]==0,eqa[[3]]==0, eqb[[1]]==0,eqb[[2]]==0,eqb[[3]]==0,eqc[[1]]==0,eqc[[2]]==0,eqc[[3]]==0,eqd[[1]]==0,eqd[[2]]==0,eqd[[3]]==0}//Simplify;
initCond = {Rx[0] == Rinit[[1]], Ry[0] == Rinit[[2]], Rz[0] == Rinit[[3]], Px[0] == Pinit[[1]], Py[0] == Pinit[[2]], Pz[0] == Pinit[[3]] , S1x[0] ==S1init[[1]], S1y[0] == S1init[[2]], S1z[0] ==S1init[[3]], S2x[0] == S2init[[1]], S2y[0] == S2init[[2]], S2z[0] == S2init[[3]]}  ;
sol=NDSolve[   system0~Join~initCond  ,  {Rx, Ry, Rz, Px, Py,Pz,S1x, S1y, S1z,S2x, S2y, S2z},{\[Lambda],0,\[Lambda]max}, PrecisionGoal->precisionGoal];

                finalvec=Re[{{Rx[\[Lambda]],Ry[\[Lambda]],Rz[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},{Px[\[Lambda]],Py[\[Lambda]],Pz[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},
                {S1x[\[Lambda]],S1y[\[Lambda]],S1z[\[Lambda]]}/.sol/.{\[Lambda]->\[Lambda]max},{S2x[\[Lambda]],S2y[\[Lambda]],S2z[\[Lambda]]}/. sol/.{\[Lambda]->\[Lambda]max}}]//N;
                (*Print["The final state is"];*)
                Print[finalvec];


                           (*Print[Plot[{Rx[\[Lambda]]/.sol,Ry[\[Lambda]]/.sol,Rz[\[Lambda]]/.sol},{\[Lambda],0,\[Lambda]max}]];
                Print[Plot[{Px[\[Lambda]]/.sol,Py[\[Lambda]]/.sol,Pz[\[Lambda]]/.sol},{\[Lambda],0,\[Lambda]max}]];
                Print[Plot[{S1x[\[Lambda]]/.sol,S1y[\[Lambda]]/.sol,S1z[\[Lambda]]/.sol},{\[Lambda],0,\[Lambda]max}]];
                Print[Plot[{S2x[\[Lambda]]/.sol,S2y[\[Lambda]]/.sol,S2z[\[Lambda]]/.sol},{\[Lambda],0,\[Lambda]max}]];*)

              (* spmg=50;
                Show[Graphics3D[{{Blue,Arrowheads[0.03],Arrow[{{0,0,0},Rinit}]},
                {Green,Arrowheads[0.03],Arrow[{{0,0,0},Pinit}]},{Brown,Arrowheads[0.03],Arrow[{{0,0,0},spmg S1init}]},
                {Magenta,Arrowheads[0.03],Arrow[{{0,0,0},spmg S2init}]}}],ParametricPlot3D[Evaluate[{Rx[\[Lambda]],Ry[\[Lambda]],Rz[\[Lambda]]}/. sol],
                {\[Lambda],\[Lambda]0,\[Lambda]max},PlotStyle->Blue],ParametricPlot3D[Evaluate[{Px[\[Lambda]],Py[\[Lambda]],Pz[\[Lambda]]}/. sol],{\[Lambda],\[Lambda]0,\[Lambda]max},PlotStyle->Green],
                ParametricPlot3D[Evaluate[spmg {S1x[\[Lambda]],S1y[\[Lambda]],S1z[\[Lambda]]}/. sol],{\[Lambda],\[Lambda]0,\[Lambda]max},PlotStyle->Brown],
                ParametricPlot3D[Evaluate[ spmg {S2x[\[Lambda]],S2y[\[Lambda]],S2z[\[Lambda]]}/. sol],{\[Lambda],\[Lambda]0,\[Lambda]max},PlotStyle->Magenta],Boxed->False]*)
                 ]
  End[]

  EndPackage[]

