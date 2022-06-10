(* ::Package:: *)

BeginPackage[ "BBHpnToolkit`SeffLfl`"]

     SeffLflow::usage = 
	"SeffLflow implements flow along Seff.L in phase space"

  Begin[ "`Private`"]

   SeffLflow[m1_, m2_, Rin_,Pin_,S1in_, S2in_,\[Lambda]max_,\[Epsilon]_]:=Module[{G,c,\[Mu],M,Q1,Q2,Linit,Jinit,RN,PN,S1N,S2N,LN,JN,
   SeffLN,Rn,Pn,S1n,S2n,Ln,Jn,SeffL,sign,f0,\[Phi]L0,acos, \[Phi]R0, \[Phi]P0,X1,Y1,R1init,P1init,R1n,P1n,R1N,P1N,\[Phi]S10,\[Phi]R10, \[Phi]P10,X2,Y2,
   R2init,P2init,R2N,P2N,\[Phi]S20,\[Phi]R20, \[Phi]P20, En,H,\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]21,\[CapitalSigma]1,\[CapitalSigma]2,a3,a2,a1,a0,A,p,q,f1,f2,f3,k,\[Alpha],B1,B2,D1,D2,\[Alpha]1sq,\[Alpha]2sq,
   B1s1,B2s1,D1s1,D2s1,\[Alpha]1s1sq,\[Alpha]2s1sq,B1s2,B2s2,D1s2,D2s2,\[Alpha]1s2sq,\[Alpha]2s2sq,f,\[Phi]L,\[Phi]R,\[Phi]P,Emat,Rfinal,Pfinal,Lfinal,\[Phi]S1,\[Phi]R1,
   \[Phi]P1,EmatS1,R1final,P1final,S1final,S2final,\[Lambda],rx,ry,rz,px,py,pz,finalvec,sphericalAngles,vectorComponents,\[Xi]2J, \[Xi]1J,
   EulMat,Rinit,Pinit,S1init,S2init,\[Lambda]0},
     G=1 ;    c = 1/Sqrt[\[Epsilon]]   ;   \[Lambda]0=0; (*set initial time to 0 such that \[Lambda]max is the flow amount*) 
    \[Mu] =m1  m2 /(m1+m2); M=m1+m2;
Q1=(1+3 m2/(4 m1));Q2=(1+3 m1/(4 m2));

Linit=Cross[Rin, Pin];
Jinit=Linit+S1in+S2in ;(* J aligned along z and stays conserved *)
RN=Norm[Rin];
PN=Norm[Pin];
S1N = Norm[S1in] ;    
S2N = Norm[S2in] ;
LN=Norm[Linit] ;
JN=Norm[Jinit];
SeffLN=(Q1 S1in+ Q2 S2in) . Linit;
sign=If[Linit . Cross[S1in,S2in]>0,1,-1];

sphericalAngles[V_]:=Module[{polarAngJ,azimuthAngJ},
polarAngJ = ToSphericalCoordinates[V][[2]];
azimuthAngJ =  ToSphericalCoordinates[V][[3]];
Return[{polarAngJ,azimuthAngJ}] ; ];
vectorComponents[{V_, \[Theta]_, \[Phi]_}]:=Module[{Vx,Vy,Vz},
Vz = V Cos[\[Theta]];
Vx=V Sin[\[Theta]]Cos[\[Phi]];  
Vy=V Sin[\[Theta]]Sin[\[Phi]];
Return[{Vx, Vy, Vz}] ; ];
{\[Xi]2J, \[Xi]1J}=sphericalAngles[Jinit]+{0,\[Pi]/2};
EulMat ={ {Cos[\[Xi]1J], Sin[\[Xi]1J],0},{-Sin[\[Xi]1J]Cos[\[Xi]2J], Cos[\[Xi]1J]Cos[\[Xi]2J], Sin[\[Xi]2J]},{Sin[\[Xi]1J]Sin[\[Xi]2J], -Cos[\[Xi]1J]Sin[\[Xi]2J], Cos[\[Xi]2J]}}      ;

Rinit= EulMat . Rin;
Pinit= EulMat . Pin;
S1init= EulMat . S1in;
S2init= EulMat . S2in;
Linit= EulMat . Linit;
Jinit= EulMat . Jinit;



(*RP space initial angles *)
f0=S1init . S2init/(Q1-Q2);
\[Phi]L0=ArcTan[Linit[[1]],Linit[[2]]];
acos[k_,a_,b_]:= If[k . (a\[Cross]b)>0,VectorAngle[a,b],2\[Pi]-VectorAngle[a,b]];
\[Phi]R0=acos[Linit, Cross[Jinit , Linit],Rinit];
\[Phi]P0=acos[Linit, Cross[Jinit , Linit],Pinit];

(*subspin sector : here we find a candidate (R1, P1) and (R2,P2) whose cross products give initial S1 and S2 vectors*)

(*subspin space initial angles for S1 *)
X1=Cross[Jinit,S1init];
Y1=Cross[S1init,X1];
Quiet[R1init={rx,ry,rz}/.Solve[Cross[{rx,ry,rz},{px,py,pz}]==S1init,{px,py,pz,rx,ry,rz}][[2]]/.{ry->1,py->1}];
Quiet[P1init={px,py,pz}/.Solve[Cross[{rx,ry,rz},{px,py,pz}]==S1init,{px,py,pz,rx,ry,rz}][[2]]/.{ry->1,py->1}];
R1N=Norm[R1init]; P1N=Norm[P1init];


\[Phi]S10=ArcTan[S1init[[1]],S1init[[2]]];
\[Phi]R10=acos[S1init, X1,R1init];
\[Phi]P10=acos[S1init, X1,P1init];

(*subspin space initial angles for S2 *)
X2=Cross[Jinit,S2init];
Y2=Cross[S2init,X2];
Quiet[R2init={rx,ry,rz}/.Solve[Cross[{rx,ry,rz},{px,py,pz}]==S2init,{px,py,pz,rx,ry,rz}][[2]]/.{ry->1,py->1}]; Quiet [P2init={px,py,pz}/.Solve[Cross[{rx,ry,rz},{px,py,pz}]==S2init,{px,py,pz,rx,ry,rz}][[2]]/.{ry->1,py->1}];
R2N=Norm[R2init];  P2N=Norm[P2init];


\[Phi]S20=ArcTan[S2init[[1]],S2init[[2]]];
\[Phi]R20=acos[S2init, X2,R2init];
\[Phi]P20=acos[S2init, X2,P2init];

(*En=\[Mu] ((Pinit.Pinit)/(2 \[Mu]^2)-(G M)/Norm[Rinit]) + ( G \[Epsilon] SeffLN)/(Norm[Rinit])^3 ;
Print["The energy for the initial data is"];Print[En];*)


(*some quantities computed from initial data*)
\[CapitalDelta]1=(1/(Q1-Q2))((1/2)(Jn^2-Ln^2-S1n^2-S2n^2)-SeffL/Q2 );
\[CapitalDelta]2=(1/(Q1-Q2))((1/2)(Jn^2-Ln^2-S1n^2-S2n^2)-SeffL/Q1 );
\[CapitalDelta]21=SeffL/(Q1 Q2 );
\[CapitalSigma]1=((Q1-Q2) \[CapitalDelta]1)/(S1n * S2n);
\[CapitalSigma]2= SeffL/( Q2 * Ln* S2n );a3 =2 Q1 Q2(Q2-Q1);a2=2(\[CapitalDelta]1+\[CapitalDelta]2)(Q1-Q2)Q1 Q2- Ln^2 (Q1-Q2)^2 -Q1^2 S1n^2-Q2^2 S2n^2;
a1=2(Q1^2 S1n^2 \[CapitalDelta]2+Q2^2 S2n^2 \[CapitalDelta]1+Q1 Q2 \[CapitalDelta]1 \[CapitalDelta]2 (Q2-Q1));
a0=Ln^2 S1n^2 S2n^2-Q1^2 S1n^2 \[CapitalDelta]2^2-Q2^2 S2n^2 \[CapitalDelta]1^2;


A=a3;
p=(3 a1 a3 -a2^2)/(3 a3^2); 
q=(2 a2^3-9 a1 a2 a3+27 a0 a3^2)/(27 a3^3);
f1=-a2/(3 a3)+2 Sqrt[-p/3]Cos[1/3 ArcCos[(3 q)/(2p) Sqrt[-3/p]]+(2 \[Pi] )/3];
f2=-a2/(3 a3)+2 Sqrt[-p/3]Cos[1/3 ArcCos[(3 q)/(2p) Sqrt[-3/p]]+(2*2 \[Pi] )/3] ;
f3=-a2/(3 a3)+2 Sqrt[-p/3]Cos[1/3 ArcCos[(3 q)/(2p) Sqrt[-3/p]]+(3*2 \[Pi] )/3];

k=(f2-f1)/(f3-f1);
\[Alpha]=sign 2/Sqrt[A(f3-f1)] EllipticF[ArcSin[Sqrt[(f0-f1)/(f2-f1)]],k];




B1=(1/2)((SeffL +Ln^2 (Q1+Q2))(Jn+Ln)+Ln(Q1  S1n^2+Q2  S2n^2+(Q1+Q2)(\[CapitalDelta]2 Q1 -\[CapitalDelta]1 Q2)));B2=(1/2)((SeffL +Ln^2 (Q1+Q2))(Jn-Ln)-Ln(Q1  S1n^2+Q2  S2n^2+(Q1+Q2)(\[CapitalDelta]2 Q1 -\[CapitalDelta]1 Q2)));D1=Ln(Ln+Jn)+(\[CapitalDelta]2 Q1 -\[CapitalDelta]1 Q2);D2=Ln(Ln-Jn)+(\[CapitalDelta]2 Q1 -\[CapitalDelta]1 Q2);


\[Alpha]1sq=((-f1+f2) (Q1-Q2))/(D1+f1 (-Q1+Q2)); \[Alpha]2sq=((-f1+f2) (Q1-Q2))/(D2+f1 (-Q1+Q2));


(*subspin space for spin S1*)
B1s1=1/2 (-S1n Q1 (Ln^2 -Jn S1n + S1n^2 +\[CapitalDelta]2 Q1)+(Jn-S1n)^2 S1n Q2-(Jn-2S1n)\[CapitalDelta]1 Q1 Q2 +(Jn-S1n)\[CapitalDelta]1 Q2^2);
B2s1=1/2 (S1n Q1 (Ln^2 +Jn S1n + S1n^2 +\[CapitalDelta]2 Q1)-(Jn+S1n)^2 S1n Q2-(Jn+2S1n)\[CapitalDelta]1 Q1 Q2 +(Jn+S1n)\[CapitalDelta]1 Q2^2);
D1s1=(S1n-Jn)S1n - \[CapitalDelta]1 Q2;
D2s1=(S1n+Jn)S1n - \[CapitalDelta]1 Q2;


\[Alpha]1s1sq=((-Q1)(f2-f1) )/(D1s1+f1 Q1); \[Alpha]2s1sq=((-Q1)(f2-f1))/(D2s1+f1 Q1);




(*orbital sector*)

Rn=RN; Pn=PN;S1n=S1N;S2n=S2N; Ln=LN;Jn=JN; SeffL=SeffLN;

f=f1+(f2-f1)(JacobiSN[1/2 Sqrt[A(f3-f1)](\[Alpha]+(\[Lambda]-\[Lambda]0)),k])^2;





(*\[Gamma]=JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Lambda]-\[Lambda]0)+EllipticF[ArcSin[Sqrt[(-f0+f1)/(f1-f2)]],k],k];*) (*in paper we call it \[Phi]p but changed it to \[Beta] to avoid confusion with \[Phi] azimuth angle of P vector*)
(*Following expressions are worked out in phiexpressions_seffl.nb*)

\[Phi]L=\[Phi]L0-1/Sqrt[A (-f1+f3)] 2 ((B1 EllipticPi[\[Alpha]1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D1+f1 (-Q1+Q2))+(B2 EllipticPi[\[Alpha]2sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D2+f1 (-Q1+Q2)))+1/Sqrt[A (-f1+f3)] 2 ((B1 EllipticPi[\[Alpha]1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k])/(D1+f1 (-Q1+Q2))+1/(D2+f1 (-Q1+Q2)) B2 EllipticPi[\[Alpha]2sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k]);






\[Phi]R=-(((Ln^2 (Q1+Q2)+SeffL+Q1 Q2 (\[CapitalDelta]1-\[CapitalDelta]2)) (\[Lambda]-\[Lambda]0))/Ln)+\[Phi]R0-1/Sqrt[A (-f1+f3)] 2 ((B1 EllipticPi[\[Alpha]1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D1+f1 (-Q1+Q2))-(B2 EllipticPi[\[Alpha]2sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D2+f1 (-Q1+Q2)))+1/Sqrt[A (-f1+f3)] 2 ((B1 EllipticPi[\[Alpha]1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k])/(D1+f1 (-Q1+Q2))-1/(D2+f1 (-Q1+Q2)) B2 EllipticPi[\[Alpha]2sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k]);

\[Phi]P=-(((Ln^2 (Q1+Q2)+SeffL+Q1 Q2 (\[CapitalDelta]1-\[CapitalDelta]2)) (\[Lambda]-\[Lambda]0))/Ln)+\[Phi]P0-1/Sqrt[A (-f1+f3)] 2 ((B1 EllipticPi[\[Alpha]1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D1+f1 (-Q1+Q2))-(B2 EllipticPi[\[Alpha]2sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D2+f1 (-Q1+Q2)))+1/Sqrt[A (-f1+f3)] 2 ((B1 EllipticPi[\[Alpha]1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k])/(D1+f1 (-Q1+Q2))-1/(D2+f1 (-Q1+Q2)) B2 EllipticPi[\[Alpha]2sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k]);


(*Show[Plot[f,{\[Lambda],0,\[Lambda]max},PlotStyle\[Rule]Red],fnumeric];
Show[Plot[Mod[\[Phi]L,2\[Pi]],{\[Lambda],0,\[Lambda]max},PlotStyle\[Rule]Red],phiLnumeric];
Show[Plot[Mod[\[Phi]R,2\[Pi]],{\[Lambda],0,\[Lambda]max},PlotStyle\[Rule]Red],phinumeric];*)


Emat={{-Sin[\[Phi]L],Cos[\[Phi]L],0},{-(((Ln^2+f (-Q1+Q2)+S1n S2n \[CapitalSigma]1+Ln S2n \[CapitalSigma]2) Cos[\[Phi]L])/(Jn Ln)),-(((Ln^2+f (-Q1+Q2)+S1n S2n \[CapitalSigma]1+Ln S2n \[CapitalSigma]2) Sin[\[Phi]L])/(Jn Ln)),Sqrt[(Jn^2 Ln^2-(Ln^2+f (-Q1+Q2)+S1n S2n \[CapitalSigma]1+Ln S2n \[CapitalSigma]2)^2)/(Jn^2 Ln^2)]},{Sqrt[(Jn^2 Ln^2-(Ln^2+f (-Q1+Q2)+S1n S2n \[CapitalSigma]1+Ln S2n \[CapitalSigma]2)^2)/(Jn^2 Ln^2)] Cos[\[Phi]L],Sqrt[(Jn^2 Ln^2-(Ln^2+f (-Q1+Q2)+S1n S2n \[CapitalSigma]1+Ln S2n \[CapitalSigma]2)^2)/(Jn^2 Ln^2)] Sin[\[Phi]L],(Ln^2-f Q1+f Q2+S1n S2n \[CapitalSigma]1+Ln S2n \[CapitalSigma]2)/(Jn Ln)}};


Rfinal=Flatten[Inverse[EulMat] . Inverse[Emat] . ({
 {Rn Cos[\[Phi]R]},
 {Rn Sin[\[Phi]R]},
 {0}
})] (* in inertial frame*);
Pfinal=Flatten[Inverse[EulMat] . Inverse[Emat] . ({
 {Pn Cos[\[Phi]P]},
 {Pn Sin[\[Phi]P]},
 {0}
})] (* in inertial frame*);
Lfinal=Cross[Rfinal,Pfinal];



(*spin sector*)
R1n=R1N; P1n=P1N;
\[Phi]S1=Jn Q2 (\[Lambda]-\[Lambda]0)+\[Phi]S10-1/Sqrt[A (-f1+f3)] 2 ((B1s1 EllipticPi[\[Alpha]1s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D1s1+f1 Q1)+(B2s1 EllipticPi[\[Alpha]2s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D2s1+f1 Q1))+(2 (B1s1 (D2s1+f1 Q1) EllipticPi[\[Alpha]1s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k]+B2s1 (D1s1+f1 Q1) EllipticPi[\[Alpha]2s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k]))/(Sqrt[A (-f1+f3)] (D1s1+f1 Q1) (D2s1+f1 Q1));






\[Phi]R1=(-Q1+Q2) S1n (\[Lambda]-\[Lambda]0)+\[Phi]R10+1/Sqrt[A (-f1+f3)] 2 ((B1s1 EllipticPi[\[Alpha]1s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D1s1+f1 Q1)-(B2s1 EllipticPi[\[Alpha]2s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D2s1+f1 Q1))+(-2 B1s1 (D2s1+f1 Q1) EllipticPi[\[Alpha]1s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k]+2 B2s1 (D1s1+f1 Q1) EllipticPi[\[Alpha]2s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k])/(Sqrt[A (-f1+f3)] (D1s1+f1 Q1) (D2s1+f1 Q1));

\[Phi]P1=(-Q1+Q2) S1n (\[Lambda]-\[Lambda]0)+\[Phi]P10+1/Sqrt[A (-f1+f3)] 2 ((B1s1 EllipticPi[\[Alpha]1s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D1s1+f1 Q1)-(B2s1 EllipticPi[\[Alpha]2s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] \[Alpha],k],k])/(D2s1+f1 Q1))+(-2 B1s1 (D2s1+f1 Q1) EllipticPi[\[Alpha]1s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k]+2 B2s1 (D1s1+f1 Q1) EllipticPi[\[Alpha]2s1sq,JacobiAmplitude[1/2 Sqrt[A (-f1+f3)] (\[Alpha]+\[Lambda]-\[Lambda]0),k],k])/(Sqrt[A (-f1+f3)] (D1s1+f1 Q1) (D2s1+f1 Q1));



EmatS1={{-Sin[\[Phi]S1],Cos[\[Phi]S1],0},{-(((f Q1 (Q1-Q2)+S1n (Q1 S1n-Q2 (S1n+S2n \[CapitalSigma]1))) Cos[\[Phi]S1])/(Jn (Q1-Q2) S1n)),-(((f Q1 (Q1-Q2)+S1n (Q1 S1n-Q2 (S1n+S2n \[CapitalSigma]1))) Sin[\[Phi]S1])/(Jn (Q1-Q2) S1n)),Sqrt[1-(f Q1 (Q1-Q2)+S1n (Q1 S1n-Q2 (S1n+S2n \[CapitalSigma]1)))^2/(Jn^2 (Q1-Q2)^2 S1n^2)]},{Sqrt[1-(f Q1 (Q1-Q2)+S1n (Q1 S1n-Q2 (S1n+S2n \[CapitalSigma]1)))^2/(Jn^2 (Q1-Q2)^2 S1n^2)] Cos[\[Phi]S1],Sqrt[1-(f Q1 (Q1-Q2)+S1n (Q1 S1n-Q2 (S1n+S2n \[CapitalSigma]1)))^2/(Jn^2 (Q1-Q2)^2 S1n^2)] Sin[\[Phi]S1],((f Q1)/S1n+(Q1 S1n-Q2 (S1n+S2n \[CapitalSigma]1))/(Q1-Q2))/Jn}};


R1final=Flatten[Inverse[EulMat] . Inverse[EmatS1] . ({
 {R1n Cos[\[Phi]R1]},
 {R1n Sin[\[Phi]R1]},
 {0}
})] (* in inertial frame*);
P1final=Flatten[Inverse[EulMat] . Inverse[EmatS1] . ({
 {P1n Cos[\[Phi]P1]},
 {P1n Sin[\[Phi]P1]},
 {0}
})] (* in inertial frame*);
S1final=Cross[R1final,P1final];
S2final=Inverse[EulMat] . Jinit-Lfinal-S1final;

finalvec={Rfinal, Pfinal, S1final, S2final}/.{\[Lambda]->\[Lambda]max}//N;
Return[finalvec];

(*Print[Plot[{Rfinal[[1]],Rfinal[[2]],Rfinal[[3]]},{\[Lambda],\[Lambda]0,\[Lambda]max}]];
Print[Plot[{Pfinal[[1]],Pfinal[[2]],Pfinal[[3]]},{\[Lambda],\[Lambda]0,\[Lambda]max}]];
Print[Plot[{S1final[[1]],S1final[[2]],S1final[[3]]},{\[Lambda],\[Lambda]0,\[Lambda]max}]];
Print[Plot[{S2final[[1]],S2final[[2]],S2final[[3]]},{\[Lambda],\[Lambda]0,\[Lambda]max}]];*)
(* Initial vectors rotated back to standard frame for plotting purposes*)



(*Rinit= Inverse[EulMat] . Rinit;
Pinit= Inverse[EulMat] . Pinit;
S1init= Inverse[EulMat] . S1init;
S2init= Inverse[EulMat] . S2init;
Show[Graphics3D[{{Blue,Arrowheads[0.03],Arrow[{{0,0,0},Rinit}]},
{Green,Arrowheads[0.03],Arrow[{{0,0,0},Pinit}]},{Brown,Arrowheads[0.03],Arrow[{{0,0,0},S1init}]},
{Magenta,Arrowheads[0.03],Arrow[{{0,0,0},S2init}]}}],ParametricPlot3D[{Rfinal[[1]],Rfinal[[2]],Rfinal[[3]]},{\[Lambda],0,\[Lambda]max},
PlotStyle->Blue],ParametricPlot3D[{Pfinal[[1]],Pfinal[[2]],Pfinal[[3]]},{\[Lambda],0,\[Lambda]max},
PlotStyle->Green],
ParametricPlot3D[{S1final[[1]],S1final[[2]],S1final[[3]]},{\[Lambda],0,\[Lambda]max},
PlotStyle->Brown],
ParametricPlot3D[{S2final[[1]],S2final[[2]],S2final[[3]]},{\[Lambda],0,\[Lambda]max},
PlotStyle->Magenta],Boxed->False]
*)




]

  End[]

  EndPackage[]

