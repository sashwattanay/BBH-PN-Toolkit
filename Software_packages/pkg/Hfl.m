(* ::Package:: *)

BeginPackage[ "pkg`Hfl`"]

     Hflow::usage = 
	"Hflow implements flow induced by H in phase space "

  Begin[ "`Private`"]

   Hflow[m1_, m2_, Rinit_,Pinit_,S1in_, S2in_,\[Lambda]max_,\[Epsilon]_]:=
               Module[{G,c, precisionGoal,kGoldstein, M, \[Mu], \[Nu],spinScalingFactor, S1init, S2init,S1ninit,S2ninit,
               rinit,pinit,Linit,Lninit,\[Delta]1,\[Delta]2,rminit, SeffdLinit, Hinit,m1N,m2N, sign, Jvec, sphericalAngles,
               vectorComponents, \[Xi]2, \[Xi]1, EulMat, R, Rx, Ry, Rz, P, Px, Py,Pz,S1, S1x, S1y, S1z, S2, S2x, S2y, S2z,
               S1n, S2n, Seff, r, p, L, Ln, Jn, SeffL, \[Sigma]1, \[Sigma]2,A, rtsol, cubiceq, x2, x3,x1,x, er2, et2, ec,\[Beta]et, ar2, n2,
               drMagnitudedtInit, tempU0, Cosu2, u0, v2, u2, findu2,solPiece0, \[CapitalUpsilon] , \[Alpha], \[Beta], cos\[Kappa]1, fcos\[Kappa]1, LinitJframe, 
               L\[Xi]10,\[Alpha]1,\[Alpha]2, \[Beta]1,\[Beta]2, solPiece1, solPiece2, L\[Xi]1sol, cos\[Kappa]2, LvecAzimuthAngle,LvecPolarAngle, Lsol, S1initJframe,
               S1\[Xi]10 ,\[Alpha]1S1, \[Alpha]2S1, \[Beta]1S1, \[Beta]2S1,solPiece1S1,solPiece2S1, S1\[Xi]1sol, cos\[Gamma], S1vecAzimuthAngle, S1vecPolarAngle,
               S1sol, LsolJframe, S2sol, JtoLframeEulMat, EulMatJtoL, rinitLframe,r\[Phi]0, IntReciOfr2,IntReciOfr3,IntReciOfr4,
               IntReciOfr5, \[Alpha]1r,\[Alpha]2r,\[Beta]1r,\[Beta]2r, solPiece1r,solPiece2r,solPiece3r,solPiece4r, r\[Phi]solLframe,rvecAzimuthAngleLframe,
               Rsol, rMagnitude, \[Phi]OffsetPvecNIFTemp, pdn, pMagnitude, \[Phi]OffsetPvecNIF,pvecAzimuthAngleLframe, spmg,
               Psol,t0, t,finalvec},


G=1 ;   c = 1/Sqrt[\[Epsilon]]   ;   \[Lambda]0=0; (*set initial time to 0 such that \[Lambda]max is the flow amount*) 

precisionGoal=Automatic (* 50*)  ;
 kGoldstein = G M \[Mu] ;  M = m1+m2 ;  \[Mu] = m1 m2/M ; \[Nu]=\[Mu]/M; spinScalingFactor = G M \[Mu] ;

S1init=S1in/spinScalingFactor  ;         
S2init=S2in/spinScalingFactor ;
(*These spin numbers (times spinSuppressFac) entered by user are unscaled quantities. But S1init and S2init are scaled (reduced) spins*)
S1ninit = Norm[S1init];
S2ninit = Norm[S2init];
rinit=Rinit/(G M)     ;     pinit=Pinit/\[Mu] ;
Linit=Cross[rinit,pinit] ;     Lninit = Norm[Linit] ;
\[Delta]1=2 m1 m2/M^2 (1+3/4 m2/m1)  ;\[Delta]2=2 m2 m1/M^2 (1+3/4 m1/m2) ;     rminit = Norm[rinit];
SeffdLinit =   Linit . (\[Delta]1  S1init  +  \[Delta]2 S2init)  ;
Hinit = Norm[pinit]^2/2 - 1/rminit + \[Epsilon]/rminit^3 SeffdLinit + (\[Epsilon] (-1)/8 (1-3 \[Nu])Norm[pinit]^4-\[Epsilon]/2  1/rminit ((3+\[Nu])Norm[pinit]^2+ \[Nu] ((rinit/rminit) . pinit)^2)+ \[Epsilon]/(2rminit^2))  ; 
m1N =m1; m2N=m2;

(*Prepare  a sign function for later use in solution*)

If [ \[Delta]2/(c^2 Norm[rinit]^3 S1ninit Lninit) Linit . Cross[S1init, S2init]   > 0, sign =1, sign =-1 ] ;





Jvec= Linit + S1init + S2init  ;  J = Norm[Jvec] ;

sphericalAngles[V_]:=Module[{polarAngJ,azimuthAngJ},
polarAngJ = ToSphericalCoordinates[Re[V]][[2]];
azimuthAngJ =  ToSphericalCoordinates[Re[V]][[3]];
Return[{polarAngJ,azimuthAngJ}] ; ];

vectorComponents[{V_, \[Theta]_, \[Phi]_}]:=Module[{Vx,Vy,Vz},
Vz = V Cos[\[Theta]];
Vx=V Sin[\[Theta]]Cos[\[Phi]];  
Vy=V Sin[\[Theta]]Sin[\[Phi]];
Return[{Vx, Vy, Vz}] ; ];

{\[Xi]2, \[Xi]1}=sphericalAngles[Jvec]+{0,\[Pi]/2};
EulMat ={ {Cos[\[Xi]1], Sin[\[Xi]1],0},{-Sin[\[Xi]1]Cos[\[Xi]2], Cos[\[Xi]1]Cos[\[Xi]2], Sin[\[Xi]2]},{Sin[\[Xi]1]Sin[\[Xi]2], -Cos[\[Xi]1]Sin[\[Xi]2], Cos[\[Xi]2]}}      ;
  (* Frame A := the one in which the above components are given.

Euler matrix which when multiplies with a column containing components in Frame A yields components in the inertial frame whose z-axis is along J vector (total angular momentum).

Gives components: Basic frame --> J-vector centered frame
*)


(*rescaled/unrescaled qunatities*)
R = {Rx, Ry,Rz} ; P ={Px,Py,Pz};
S1= {S1x, S1y,S1z} ;    S2= {S2x, S2y,S2z} ;
S1n = Norm[S1init] ;     S2n = Norm[S2init] ;

Seff= (\[Delta]1 S1+ \[Delta]2 S2);
     r=R/(G M);       p=P/\[Mu];
L=Cross[r,p] ;
Ln=Norm[Linit] ;
Jn=Norm[Linit+S1init+S2init];
SeffL= Seff . L ;


(*Evaluating \[Sigma]1 and \[Sigma]2*)
{Rx,Ry,Rz} = Rinit  ;
{Px,Py,Pz}=Pinit  ;
{S1x,S1y,S1z} = S1init ; 
{S2x,S2y,S2z} = S2init  ;

\[Sigma]1=S1 . S2/(S1n S2n) -Ln/S2n (\[Delta]1-\[Delta]2)/\[Delta]2  S1 . L/(Ln S1n) ;   (*\[Sigma]1=Cos \[Gamma] -L/S2(Q1-Q2)/Q2Cos \[Kappa]1 *)
\[Sigma]2=  S2 . L/(Ln S2n) + (\[Delta]1 S1n)/(\[Delta]2 S2n)  S1 . L/(Ln S1n) ;           (*\[Sigma]2= Cos \[Kappa]2 +(Q1 S1)/(Q2 S2) Cos \[Kappa]1*)
Clear[Rx, Ry, Rz, Px, Py, Pz, S1x,S1y,S1z,S2x,S2y,S2z];


(*finding roots of cubic*)
A=2 Ln S1n \[Delta]1 (\[Delta]2-\[Delta]1);cubiceq=A x^3 -(Ln^2 (\[Delta]1-\[Delta]2)^2+2 \[Delta]2 Ln \[Sigma]2 S2n  (\[Delta]2-\[Delta]1)+S1n^2 \[Delta]1^2+2 \[Delta]1 \[Delta]2 \[Sigma]1 S1n S2n +\[Delta]2^2  S2n^2  )x^2+ (2  \[Delta]2 S2n (Ln \[Sigma]1(\[Delta]2-\[Delta]1) +\[Sigma]2( \[Delta]1 S1n+ \[Delta]2 \[Sigma]1 S2n) )) x+-S2n^2 \[Delta]2^2 (\[Sigma]1^2+\[Sigma]2^2-1);(*=A(x -x1)(x -x2)(x -x3)*)
 
(*evaluate the cubic and find the roots*)
rtsol=Solve[cubiceq==0,x] ;
x2=x/.rtsol[[1]] ; 
x3=x/.rtsol[[2]] ;
x1=x/.rtsol[[3]];
{x2, x3, x1} = Sort[{x2, x3, x1}] ;



(*cos \[Kappa]1 solution*)
er2 = (1 + 2Hinit Lninit^2   +     2 (6-\[Nu])(-Hinit)\[Epsilon] - 5 (3-\[Nu])Hinit^2 Lninit^2 \[Epsilon]    + 8 (1+ Hinit Lninit^2) SeffdLinit/Lninit^2 Hinit \[Epsilon]  )^(1/2)  ;
et2 = ( 1 + 2 Hinit Lninit^2   +  4 (1-\[Nu])Hinit \[Epsilon] +(17-7\[Nu])Hinit^2 Lninit^2 \[Epsilon]    +4 SeffdLinit/Lninit^2 Hinit \[Epsilon]   )^(1/2) ;
ec =  (-8 Hinit Lninit^2-16 Hinit^2 Lninit^4+2 Hinit SeffdLinit+4 Hinit^2 Lninit^2 SeffdLinit+3 Hinit Lninit^2 \[Nu]+6 Hinit^2 Lninit^4 \[Nu])/(Lninit^2 Sqrt[1+2 Hinit Lninit^2]) ;  (*er = et + \[Epsilon] ec *)
\[Beta]et = et2/(1+(1-et2^2)^(1/2)) ;
ar2 = -1/(2 Hinit) (1 -1/2 (7-\[Nu])  (-Hinit)\[Epsilon] - 2 SeffdLinit/Lninit^2 Hinit \[Epsilon]  ) ;
n2 = (-2 Hinit)^(3/2)   (1 +  (-2 Hinit)/8 \[Epsilon] (-15+\[Nu])  )/(G M);    (*unscaled n*)


drMagnitudedtInit = Pinit . Rinit(2+\[Epsilon] (-1+3 \[Nu]) Norm[Pinit]^2-(2 \[Epsilon] (3+2 \[Nu]))/Norm[Rinit]) ;
tempU0 = ArcCos[Cosu2/.(Solve[ Norm[rinit ] == ar2 (1-er2 Cosu2),Cosu2]    //N)[[1]] ]    ;
If[drMagnitudedtInit  > 0,
u0 = tempU0 ; ,
u0 = -tempU0 ;];
t0 =- (t0/.Solve[n2 t0 == u0 - et2 Sin[u0], t0] [[1]]   );          (*This is t0 in n(t-t0) = u - et sin[u]*)





v2 = u2 + 2 ArcTan[ (\[Beta]et Sin[u2])/(1-\[Beta]et Cos[u2])];            (* PN extension of v= 2 ArcTan[Sqrt[(1+e)/(1-e)] Tan[u/2]] but without the ArcTan issues*)
findu2[t_]:=FindRoot[n2 (t-t0) == u2 - et2 Sin[u2],{u2, n2 (t-t0), n2 (t-t0) - 2 \[Pi], n2 (t-t0)+2 \[Pi]}, PrecisionGoal->precisionGoal];

\[Beta] = ((x3-x2)/(x1-x2))^(1/2);
\[Alpha] =  sign 2/Sqrt[A(x1-x2)] EllipticF[ArcSin[Sqrt[( S1init . Linit/(Lninit S1ninit)-x2)/(x3-x2)]],\[Beta]^2]   ;
solPiece0 = ((-2+2 et2^2-9 ec et2 \[Epsilon]) (-((Sqrt[-1+et2^2] v2)/(2 Sqrt[1-et2^2]))))/(-1+et2^2)^(5/2)+((2 et2-2 et2^3+6 ec \[Epsilon]+3 ec et2^2 \[Epsilon]+et2 (-2 et2+2 et2^3-3 ec \[Epsilon]-6 ec et2^2 \[Epsilon]) Cos[u2]) Sin[u2])/(2 (-1+et2^2)^2 (-1+et2 Cos[u2])^2)  ;
\[CapitalUpsilon] = 1/2 Sqrt[A(x1-x2)](\[Alpha]+solPiece0/(c^2 (n2 G M) ar2^3)) ;
cos\[Kappa]1=x2+(x3-x2)JacobiSN[\[CapitalUpsilon],\[Beta]^2]^2;
fcos\[Kappa]1[t_]:=cos\[Kappa]1/.findu2[t];


(*L vector solution*)

LinitJframe = EulMat . Linit      ;
L\[Xi]10 =   sphericalAngles[LinitJframe][[2]]+ \[Pi]/2 ;
\[Alpha]1 =-((\[Delta]2 (J+Lninit+S2ninit \[Sigma]2))/(S1ninit (\[Delta]1-\[Delta]2))) ;
\[Alpha]2 =-((\[Delta]2 (-J+Lninit+S2ninit \[Sigma]2))/(S1ninit (\[Delta]1-\[Delta]2))) ;
\[Beta]1 =-(1/(2 S1ninit (\[Delta]1-\[Delta]2)))\[Delta]2 (Lninit^2 \[Delta]1+J^2 \[Delta]2+S1ninit (\[Delta]1-\[Delta]2) (S1ninit+S2ninit \[Sigma]1)+Lninit S2ninit \[Delta]1 \[Sigma]2+J (Lninit (\[Delta]1+\[Delta]2)+S2ninit \[Delta]2 \[Sigma]2)) ;
\[Beta]2 =-(1/(2 S1ninit (\[Delta]1-\[Delta]2)))\[Delta]2 (Lninit^2 \[Delta]1+J^2 \[Delta]2+S1ninit (\[Delta]1-\[Delta]2) (S1ninit+S2ninit \[Sigma]1)+Lninit S2ninit \[Delta]1 \[Sigma]2-J (Lninit (\[Delta]1+\[Delta]2)+S2ninit \[Delta]2 \[Sigma]2)) ;
solPiece1 =(2/(A(x1-x2))^(1/2)) (( \[Beta]1  EllipticPi[(x2-x3)/(\[Alpha]1+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]1+x2)-( \[Beta]2 EllipticPi[(x2-x3)/(\[Alpha]2+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]2+x2))  ; 
solPiece2 =L\[Xi]10 -(solPiece1/.u2->u0 )    ;
L\[Xi]1sol=solPiece1 + solPiece2  ;
cos\[Kappa]2[t_]:= \[Sigma]2 -  \[Delta]1 S1ninit/(\[Delta]2 S2ninit) fcos\[Kappa]1[t];
LvecAzimuthAngle[t_]:=-\[Pi]/2+L\[Xi]1sol/.findu2[t];
LvecPolarAngle[t_]:= ArcCos[(Lninit^2+ Lninit S1ninit  fcos\[Kappa]1[t] + Lninit S2ninit  cos\[Kappa]2[t]  )/(Lninit J)]; 
Lsol[t_]:=Module[{Lx,Ly,Lz},
{Lx,Ly,Lz} = vectorComponents[{Lninit, LvecPolarAngle[t],LvecAzimuthAngle[t]}];
Return[ kGoldstein(Inverse[EulMat] . {Lx,Ly,Lz}) ]; ];


(*S1 and S2 solution*)
S1initJframe = EulMat . S1init      ;
S1\[Xi]10 =   sphericalAngles[S1initJframe][[2]]+ \[Pi]/2 ;

\[Alpha]1S1 =  (\[Delta]2 (-J+S1ninit+S2ninit \[Sigma]1))/(Lninit \[Delta]1)  ;
\[Alpha]2S1 =(\[Delta]2 (J+S1ninit+S2ninit \[Sigma]1))/(Lninit \[Delta]1) ;
\[Beta]1S1 =-(1/(2 Lninit \[Delta]1))\[Delta]2 (Lninit^2 \[Delta]1+(S1ninit (\[Delta]1-\[Delta]2)+J \[Delta]2) (-J+S1ninit+S2ninit \[Sigma]1)+Lninit S2ninit \[Delta]1 \[Sigma]2)  ;
\[Beta]2S1 =1/(2 Lninit \[Delta]1) \[Delta]2 (-Lninit^2 \[Delta]1+(-S1ninit \[Delta]1+(J+S1ninit) \[Delta]2) (J+S1ninit+S2ninit \[Sigma]1)-Lninit S2ninit \[Delta]1 \[Sigma]2) ;

solPiece1S1 =(2/(A(x1-x2))^(1/2)) (( \[Beta]1S1  EllipticPi[(x2-x3)/(\[Alpha]1S1+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]1S1+x2)-( \[Beta]2S1 EllipticPi[(x2-x3)/(\[Alpha]2S1+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]2S1+x2))   + solPiece0/((n2 G M) ar2^3)  (J \[Delta]2)/c^2  ;
solPiece2S1 =S1\[Xi]10 -(solPiece1S1/.u2->u0 )    ;
S1\[Xi]1sol=solPiece1S1 + solPiece2S1  ;
cos\[Gamma][t_]:= \[Sigma]1 +   Lninit (\[Delta]1-\[Delta]2)/(\[Delta]2 S2ninit) fcos\[Kappa]1[t];
S1vecAzimuthAngle[t_]:=-\[Pi]/2+S1\[Xi]1sol/.findu2[t];
S1vecPolarAngle[t_]:= ArcCos[(S1ninit^2+ Lninit S1ninit  fcos\[Kappa]1[t] + S1ninit S2ninit  cos\[Gamma][t]  )/(S1ninit J)]; 
S1sol[t_]:=Module[{S1xx,S1yy,S1zz},
{S1xx,S1yy,S1zz} = vectorComponents[{S1ninit, S1vecPolarAngle[t],S1vecAzimuthAngle[t]}];
Return[ spinScalingFactor(Inverse[EulMat] . {S1xx,S1yy,S1zz}) ]; ];
S2sol[t_]:=spinScalingFactor(Jvec - Lsol[t]/kGoldstein-S1sol[t]/spinScalingFactor);

(*R solution*)
LsolJframe[t_]:=Module[{Lx,Ly,Lz},
{Lx,Ly,Lz} = vectorComponents[{Lninit, LvecPolarAngle[t],LvecAzimuthAngle[t]}];
Return[ kGoldstein({Lx,Ly,Lz}) ]; ];

JtoLframeEulMat[t_]:=Module[{Lx,Ly,Lz, \[Xi]1L, \[Xi]2L},
{Lx,Ly,Lz} = LsolJframe[t];
{\[Xi]2L, \[Xi]1L}=sphericalAngles[{Lx,Ly,Lz}]+{0,\[Pi]/2};
EulMatJtoL ={ {Cos[\[Xi]1L], Sin[\[Xi]1L],0},{-Sin[\[Xi]1L]Cos[\[Xi]2L], Cos[\[Xi]1L]Cos[\[Xi]2L], Sin[\[Xi]2L]},{Sin[\[Xi]1L]Sin[\[Xi]2L], -Cos[\[Xi]1L]Sin[\[Xi]2L], Cos[\[Xi]2L]}}      ;
Return[EulMatJtoL];];
rinitLframe=JtoLframeEulMat[0] . EulMat . rinit       ;
r\[Phi]0 =   sphericalAngles[rinitLframe ][[2]] ;


IntReciOfr2 = ((v2 (-1+et2^2-2 ec et2 \[Epsilon]))/Sqrt[1-et2^2]+(2 ec \[Epsilon] Sin[u2])/(-1+et2 Cos[u2]))/(ar2^2 (-1+et2^2) G M n2)     (*Integral of 1/r^2*)  ;
IntReciOfr3 =((v2 (2-2 et2^2+9 ec et2 \[Epsilon]))/Sqrt[1-et2^2]+((2 et2-2 et2^3+6 ec \[Epsilon]+3 ec et2^2 \[Epsilon]+et2 (-2 et2+2 et2^3-3 ec \[Epsilon]-6 ec et2^2 \[Epsilon]) Cos[u2]) Sin[u2])/(-1+et2 Cos[u2])^2)/(2 ar2^3 (-1+et2^2)^2 G M n2) (*Integral of 1/r^3*)  ;

IntReciOfr4 =1/(6 ar2^4 G M n2) (-((3 v2 (-2+et2^2+et2^4-16 ec et2 \[Epsilon]-4 ec et2^3 \[Epsilon]))/(1-et2^2)^(7/2))+(8 ec \[Epsilon] Sin[u2])/((-1+et2^2) (-1+et2 Cos[u2])^3)+((3 et2-3 et2^3+8 ec \[Epsilon]+12 ec et2^2 \[Epsilon]) Sin[u2])/((-1+et2)^2 (1+et2)^2 (-1+et2 Cos[u2])^2)+((9 et2-9 et2^3+8 ec \[Epsilon]+52 ec et2^2 \[Epsilon]) Sin[u2])/((-1+et2^2)^3 (-1+et2 Cos[u2])))     ;
IntReciOfr5 =   1/(ar2^5 (G M n2)) (1/(96 (-1+et2^2)^4) (-((12 v2 (-8-4 et2^2+12 et2^4-100 ec et2 \[Epsilon]-75 ec et2^3 \[Epsilon]))/Sqrt[1-et2^2])+1/(-1+et2 Cos[u2])^4 (288 et2-64 et2^3-136 et2^5-88 et2^7+480 ec \[Epsilon]+1920 ec et2^2 \[Epsilon]+2620 ec et2^4 \[Epsilon]+230 ec et2^6 \[Epsilon]+5 et2 (-144 et2+124 et2^3+4 et2^5+16 et2^7-144 ec \[Epsilon]-1182 ec et2^2 \[Epsilon]-169 ec et2^4 \[Epsilon]-80 ec et2^6 \[Epsilon]) Cos[u2]+2 et2^2 (152 et2-124 et2^3-28 et2^5+120 ec \[Epsilon]+1360 ec et2^2 \[Epsilon]+95 ec et2^4 \[Epsilon]) Cos[2 u2]-44 et2^4 Cos[3 u2]+28 et2^6 Cos[3 u2]+16 et2^8 Cos[3 u2]-30 ec et2^3 \[Epsilon] Cos[3 u2]-415 ec et2^5 \[Epsilon] Cos[3 u2]-80 ec et2^7 \[Epsilon] Cos[3 u2]) Sin[u2]) )  ;
\[Beta]1r =  -(1/(2 (S1ninit \[Delta]1-S1ninit \[Delta]2) \[Epsilon]))\[Delta]2 (-J Lninit \[Delta]1 \[Epsilon]+Lninit^2 \[Delta]1 \[Epsilon]+S1ninit^2 \[Delta]1 \[Epsilon]+J^2 \[Delta]2 \[Epsilon]-J Lninit \[Delta]2 \[Epsilon]-S1ninit^2 \[Delta]2 \[Epsilon]+S1ninit S2ninit \[Delta]1 \[Epsilon] \[Sigma]1-S1ninit S2ninit \[Delta]2 \[Epsilon] \[Sigma]1+Lninit S2ninit \[Delta]1 \[Epsilon] \[Sigma]2-J S2ninit \[Delta]2 \[Epsilon] \[Sigma]2);
\[Beta]2r=-(1/(2 (S1ninit \[Delta]1-S1ninit \[Delta]2) \[Epsilon]))\[Delta]2 (J Lninit \[Delta]1 \[Epsilon]+Lninit^2 \[Delta]1 \[Epsilon]+S1ninit^2 \[Delta]1 \[Epsilon]+J^2 \[Delta]2 \[Epsilon]+J Lninit \[Delta]2 \[Epsilon]-S1ninit^2 \[Delta]2 \[Epsilon]+S1ninit S2ninit \[Delta]1 \[Epsilon] \[Sigma]1-S1ninit S2ninit \[Delta]2 \[Epsilon] \[Sigma]1+Lninit S2ninit \[Delta]1 \[Epsilon] \[Sigma]2+J S2ninit \[Delta]2 \[Epsilon] \[Sigma]2) ;
\[Alpha]1r = (J \[Delta]2-Lninit \[Delta]2-S2ninit \[Delta]2 \[Sigma]2)/(S1ninit \[Delta]1-S1ninit \[Delta]2) ;
\[Alpha]2r = (-J \[Delta]2-Lninit \[Delta]2-S2ninit \[Delta]2 \[Sigma]2)/(S1ninit \[Delta]1-S1ninit \[Delta]2);
solPiece1r =(2/(A(x1-x2))^(1/2)) (( \[Beta]1r  EllipticPi[(x2-x3)/(\[Alpha]1r+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]1r+x2)+( \[Beta]2r EllipticPi[(x2-x3)/(\[Alpha]2r+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]2r+x2))  ;
solPiece2r =-(1/2) IntReciOfr5 Lninit \[Epsilon]^2 (-1+3 \[Nu]) (2 SeffdLinit+Lninit^2 \[Nu])+1/2 IntReciOfr4 Lninit \[Epsilon]^2 (-6+17 \[Nu]+3 \[Nu]^2)+1/2 IntReciOfr2 Lninit (2-Hinit^2 \[Epsilon]^2 (1-3 \[Nu])^2+2 Hinit \[Epsilon] (-1+3 \[Nu]))-1/Lninit IntReciOfr3 \[Epsilon] (-SeffdLinit+Lninit^2 (4+\[Delta]1+\[Delta]2+4 Hinit \[Epsilon]-2 \[Nu]-13 Hinit \[Epsilon] \[Nu]+3 Hinit \[Epsilon] \[Nu]^2)+Lninit S2ninit \[Delta]2 \[Sigma]2);
solPiece3r =  solPiece1r + solPiece2r ;
solPiece4r =r\[Phi]0 -(solPiece3r/.u2->u0 )    ;
r\[Phi]solLframe=solPiece3r + solPiece4r  ;
rvecAzimuthAngleLframe[t_]:=r\[Phi]solLframe/.findu2[t];
rMagnitude[t_] := ar2 (1-er2 Cos[u2])/.findu2[t];
Rsol[t_]:=Module[{rx,ry,rz},
{rx,ry,rz} = vectorComponents[{rMagnitude[t], \[Pi]/2,rvecAzimuthAngleLframe[t]}];
Return[ G M(Inverse[EulMat] . Inverse[JtoLframeEulMat[t]] . {rx,ry,rz}) ]; ];








(* P solution*)

pdn[t_] :=  ((2 Hinit + \[Epsilon](1- 3 \[Nu])Hinit^2) + (2(1 + \[Epsilon] (4-\[Nu])Hinit))/rMagnitude[t] + (-Lninit^2+\[Epsilon] (6+\[Nu]))/rMagnitude[t]^2 + (-\[Epsilon] \[Nu] (Lninit^2+ 2SeffdLinit/\[Nu]) )/rMagnitude[t]^3  )^(1/2) ;
pMagnitude[t_]:=( pdn[t]^2+Lninit^2/rMagnitude[t]^2)^(1/2) ;
\[Phi]OffsetPvecNIFTemp[t_] :=  ArcSin[Lninit/(rMagnitude[t]pMagnitude[t])] ;
\[Phi]OffsetPvecNIF[t_]:=Module[{tempVar} ,
If[pdn[t] >0,
tempVar = \[Phi]OffsetPvecNIFTemp[t] ,
tempVar = \[Phi]OffsetPvecNIFTemp[t]+2(\[Pi]/2 - \[Phi]OffsetPvecNIFTemp[t]) ;];
Return[tempVar] ;  ];
pvecAzimuthAngleLframe[t_]:= rvecAzimuthAngleLframe[t]+\[Phi]OffsetPvecNIF[t];
Psol[t_]:=Module[{px,py,pz},
{px,py,pz} = vectorComponents[{pMagnitude[t], \[Pi]/2,pvecAzimuthAngleLframe[t]}];
Return[ \[Mu](Inverse[EulMat] . Inverse[JtoLframeEulMat[t]] . {px,py,pz}) ]; ];


                

               
                finalvec=Re[{Rsol[\[Lambda]max-\[Lambda]0],Psol[\[Lambda]max-\[Lambda]0],
                S1sol[\[Lambda]max-\[Lambda]0],S2sol[\[Lambda]max-\[Lambda]0]}]//N;
                (*Print["The final state is"];*)
                Return[finalvec];


(*Print[Plot[{Rsol[\[Lambda]][[1]],Rsol[\[Lambda]][[2]],Rsol[\[Lambda]][[3]]},{\[Lambda],0,\[Lambda]max-\[Lambda]0}]];
Print[Plot[{Psol[\[Lambda]][[1]],Psol[\[Lambda]][[2]],Psol[\[Lambda]][[3]]},{\[Lambda],0,\[Lambda]max-\[Lambda]0}]];
Print[Plot[{S1sol[\[Lambda]][[1]],S1sol[\[Lambda]][[2]],S1sol[\[Lambda]][[3]]},{\[Lambda],0,\[Lambda]max-\[Lambda]0}]];
Print[Plot[{S2sol[\[Lambda]][[1]],S2sol[\[Lambda]][[2]],S2sol[\[Lambda]][[3]]},{\[Lambda],0,\[Lambda]max-\[Lambda]0}]];
*)
(*spmg=50;

            Show[Graphics3D[{{Blue,Arrowheads[0.03],Arrow[{{0,0,0},Rinit}]},
{Green,Arrowheads[0.03],Arrow[{{0,0,0},Pinit}]},{Brown,Arrowheads[0.03],Arrow[{{0,0,0}, spmg S1in}]},
{Magenta,Arrowheads[0.03],Arrow[{{0,0,0}, spmg S2in}]}}],
ParametricPlot3D[{Rsol[t][[1]],Rsol[t][[2]],Rsol[t][[3]]},{t,0,\[Lambda]max-\[Lambda]0},
PlotStyle->Blue],ParametricPlot3D[{Psol[t][[1]],Psol[t][[2]],Psol[t][[3]]},{t,0,\[Lambda]max-\[Lambda]0},
PlotStyle->Green],
ParametricPlot3D[spmg {S1sol[t][[1]],S1sol[t][[2]],S1sol[t][[3]]},{t,0,\[Lambda]max-\[Lambda]0},
PlotStyle->Brown],
ParametricPlot3D[ spmg {S2sol[t][[1]],S2sol[t][[2]],S2sol[t][[3]]},{t,0,\[Lambda]max-\[Lambda]0},
PlotStyle->Magenta],Boxed->False]*)
                 ]
  End[]

  EndPackage[]

