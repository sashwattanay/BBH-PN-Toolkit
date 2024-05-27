BeginPackage[ "BBHpnToolkit`Hfl`"]

     Hflow::usage = 
	"Hflow implements flow induced by H in phase space "

  Begin[ "`Private`"]

   Hflow[G_,m1_, m2_, Rinit_,Pinit_,S1in_, S2in_,\[Lambda]max_,\[Epsilon]_]:=
               Module[{c,\[Lambda]0,precisionGoal,M,\[Mu],\[Nu],s1vec,s2vec,s1,s2,rvec,r,
               pvec,p,lvec,l,\[Delta]1,\[Delta]2,dotlseff,h,\[CapitalEpsilon],subslist,sign,
               jvec,j,sphericalAngles,vectorComponents,\[Xi]1,\[Xi]2,EulMat,\[CapitalSigma]1,\[CapitalSigma]\[CapitalSigma]1,\[CapitalSigma]2,\[CapitalSigma]\[CapitalSigma]2,A,cubiceq,x,x1,
               x2,x3,rtsol,er,et,eN,\[Beta]eN,arN,ar,nN,n,nunscaled,drMagnitudedtInit,tempU0,Cosu2,
               u0,t0,v,u,findu,\[Beta],\[Alpha],cos\[Kappa]10,\[CapitalUpsilon],cos\[Kappa]1,fcos\[Kappa]1,lvecJframe,l\[Xi]10,\[Alpha]1,\[Alpha]2,\[Beta]1,\[Beta]2,solPiece1,
               solPiece2,l\[Xi]1sol,cos\[Kappa]2,lvecAzimuthAngle,lvecPolarAngle,Lsol,s1vecJframe,s2vecJframe,
               s1\[Xi]10,s2\[Xi]10,\[Alpha]1s1,\[Alpha]1s2,\[Alpha]2s1,\[Alpha]2s2,\[Beta]1s1,\[Beta]1s2,\[Beta]2s1,\[Beta]2s2,solPiece1s1,solPiece1s2,solPiece2s1,
               solPiece2s2,s1\[Xi]1sol,s2\[Xi]1sol,cos\[Gamma],s1vecAzimuthAngle,s2vecAzimuthAngle,s1vecPolarAngle,
               s2vecPolarAngle,S1sol,S2sol,LsolJframe,JtoLframeEulMat,EulMatJtoL,rLframe,r\[Phi]0,
               vN,e\[Theta],\[Beta]e\[Theta],v\[Theta],e\[Theta]prime,\[Beta]e\[Theta]prime,v\[Theta]prime,solPiece1r,solPiece2r,solPiece3r,solPiece4r,
               r\[Phi]solLframe,rvecAzimuthAngleLframe,rMagnitude,Rsol,pdn,pMagnitude,\[Phi]OffsetPvecNIFTemp,
               \[Phi]OffsetPvecNIF,pvecAzimuthAngleLframe,Psol,t,finalvec,pdninit,periodpdn,findpdnsign,first0},

c = 1/Sqrt[\[Epsilon]];   
\[Lambda]0=0; (*set initial time to 0 such that \[Lambda]max is the flow amount*) 
precisionGoal = Automatic (* 50*)  ;
M = m1+m2 ;  \[Mu] = m1 m2/M ; \[Nu]=\[Mu]/M;
s1vec = S1in/(\[Mu] G M);
s2vec = S2in/(\[Mu] G M);
(*These spin bers (times spinSuppressFac) entered by user are unscaled quantities. But s1 and s2 are scaled (reduced) spins*)
s1 = Norm[s1vec];
s2 = Norm[s2vec];
rvec = Rinit/(G M);
r = Norm[rvec];
pvec = Pinit/\[Mu];
p = Norm[pvec];
lvec = Cross[rvec,pvec];
l = Norm[lvec] ;
\[Delta]1 = 2 m1 m2/M^2 (1+3/4 m2/m1);
\[Delta]2 = 2 m2 m1/M^2 (1+3/4 m1/m2);
dotlseff = lvec . (\[Delta]1  s1vec  +  \[Delta]2 s2vec);
\[CapitalEpsilon] = -h;
h = p^2/2 - 1/r + \[Epsilon]/r^3 dotlseff + \[Epsilon](1/8 (3\[Nu]-1)p^4+1/(2r^2)-1/(2r) ((3+\[Nu])p^2+\[Nu] Dot[rvec/r,pvec]^2))  ; 

(*Prepare  a sign function for later use in solution*)

If [ \[Delta]2/(c^2 r^3 s1 l) lvec . Cross[s1vec, s2vec] > 0, sign =1, sign =-1 ] ;

jvec = lvec + s1vec + s2vec  ;  j = Norm[jvec] ;

sphericalAngles[V_]:=Module[{polarAngJ,azimuthAngJ},
polarAngJ = ToSphericalCoordinates[Re[V]][[2]];
azimuthAngJ =  ToSphericalCoordinates[Re[V]][[3]];
Return[{polarAngJ,azimuthAngJ}] ; ];

vectorComponents[{V_, \[Theta]_, \[Phi]_}]:= Module[{Vx,Vy,Vz},
Vz = V Cos[\[Theta]];
Vx=V Sin[\[Theta]]Cos[\[Phi]];  
Vy=V Sin[\[Theta]]Sin[\[Phi]];
Return[{Vx, Vy, Vz}] ; ];

{\[Xi]2, \[Xi]1}=sphericalAngles[jvec]+{0,\[Pi]/2};
EulMat ={{Cos[\[Xi]1], Sin[\[Xi]1],0},{-Sin[\[Xi]1]Cos[\[Xi]2], Cos[\[Xi]1]Cos[\[Xi]2], Sin[\[Xi]2]},{Sin[\[Xi]1]Sin[\[Xi]2], -Cos[\[Xi]1]Sin[\[Xi]2], Cos[\[Xi]2]}}      ;
(* Frame A := the one in which the above components are given.
Euler matrix which when multiplies with a column containing components in Frame A yields components in the inertial frame whose z-axis is along J vector (total angular momentum).
Gives components: Basic frame --> J-vector centered frame
*)

(*Evaluating \[Sigma]1 and \[Sigma]2*)
\[CapitalSigma]1 = N[s1vec . s2vec/(s1 s2) - l/s2 (\[Delta]1-\[Delta]2)/\[Delta]2  s1vec . lvec/(l s1)];   (*\[CapitalSigma]1 = Cos \[Gamma] -L/S2(Q1-Q2)/Q2Cos \[Kappa]1 *)
\[CapitalSigma]2 = N[s2vec . lvec/(l s2) + (\[Delta]1 s1)/(\[Delta]2 s2)  s1vec . lvec/(l s1)];     (*\[CapitalSigma]2 = Cos \[Kappa]2 +(Q1 S1)/(Q2 S2) Cos \[Kappa]1*)

(*finding roots of cubic*)
A = 2 l s1 \[Delta]1 (\[Delta]2-\[Delta]1);
cubiceq = A x^3 -(l^2 (\[Delta]1-\[Delta]2)^2+2 \[Delta]2 l \[CapitalSigma]2 s2(\[Delta]2-\[Delta]1)+s1^2 \[Delta]1^2+2 \[Delta]1 \[Delta]2 \[CapitalSigma]1 s1 s2 +\[Delta]2^2  s2^2  )x^2+ (2 \[Delta]2 s2 (l \[CapitalSigma]1(\[Delta]2-\[Delta]1) +\[CapitalSigma]2(\[Delta]1 s1+ \[Delta]2 \[CapitalSigma]1 s2) )) x+-s2^2 \[Delta]2^2 (\[CapitalSigma]1^2+\[CapitalSigma]2^2-1);(*=A(x -x1)(x -x2)(x -x3)*)
 
(*evaluate the cubic and find the roots*)
rtsol = Solve[cubiceq==0,x];
x2 = x/.rtsol[[1]]; 
x3 = x/.rtsol[[2]];
x1 = x/.rtsol[[3]];
{x2, x3, x1} = Sort[{x2, x3, x1}];

(*cos \[Kappa]1 solution*)
(*qkp OPN*)
eN = Sqrt[1-2 l^2 \[CapitalEpsilon]];
\[Beta]eN = eN/(1+(1-eN^2)^(1/2));
arN = 1/(2 \[CapitalEpsilon]);
nN = (2\[CapitalEpsilon])^(3/2)/(G M);

(*qkp 1.5PN*)
ar = 1/(2\[CapitalEpsilon])+ (dotlseff \[Epsilon])/l^2+\[Epsilon]/4 (-7+\[Nu]);
er = (1 - 2 l^2 \[CapitalEpsilon] + \[Epsilon](-((8 dotlseff \[CapitalEpsilon])/l^2) + 8dotlseff \[CapitalEpsilon]^2) + \[Epsilon](12 \[CapitalEpsilon]-15 l^2 \[CapitalEpsilon]^2-2 \[CapitalEpsilon] \[Nu]+5 l^2 \[CapitalEpsilon]^2 \[Nu]))^(1/2);
et = (1 - 2 l^2 \[CapitalEpsilon] - (4 dotlseff \[Epsilon] \[CapitalEpsilon])/l^2 + \[Epsilon](l^2 \[CapitalEpsilon]^2 (17 - 7\[Nu]) + 4\[CapitalEpsilon](\[Nu]-1)))^(1/2);
n = ((2\[CapitalEpsilon])^(3/2)+ \[Epsilon] \[CapitalEpsilon]^(5/2) (-(15/Sqrt[2])+\[Nu]/Sqrt[2]));
nunscaled = ((2\[CapitalEpsilon])^(3/2)+ \[Epsilon] \[CapitalEpsilon]^(5/2) (-(15/Sqrt[2])+\[Nu]/Sqrt[2]))/(G M);  (*unscaled n*)

v = u + 2ArcTan[ (\[Beta]eN Sin[u])/(1-\[Beta]eN Cos[u])];            (* PN extension of v= 2 ArcTan[Sqrt[(1+e)/(1-e)] Tan[u/2]] but without the ArcTan issues*)

drMagnitudedtInit = Pinit . Rinit(2+\[Epsilon] (-1+3 \[Nu]) Norm[Pinit]^2-(2 \[Epsilon] (3+2 \[Nu]))/Norm[Rinit]);
tempU0 = ArcCos[Cosu2/.(Solve[r == ar (1-er Cosu2),Cosu2]//N)[[1]]];
If[drMagnitudedtInit  > 0,
u0 = tempU0;,
u0 = -tempU0;];
t0 =- (t0/.Solve[ nunscaled t0 == u0 - et Sin[u0], t0] [[1]]);          (*This is t0 in n(t-t0) = u - et sin[u]*)

findu[t_]:=FindRoot[nunscaled (t-t0) == u - et Sin[u],{u, nunscaled (t-t0), nunscaled (t-t0) - 2 \[Pi], nunscaled (t-t0)+2 \[Pi]}, PrecisionGoal->precisionGoal]; (*u is 1.5PN*)

cos\[Kappa]10 = Dot[lvec,s1vec]/(Norm[lvec]Norm[s1vec]);

\[Beta] = ((x3-x2)/(x1-x2))^(1/2);
\[Alpha] =  sign 2/Sqrt[A(x1-x2)] EllipticF[ArcSin[Sqrt[(cos\[Kappa]10-x2)/(x3-x2)]],\[Beta]^2];
\[CapitalUpsilon] = 1/2 Sqrt[A(x1-x2)](\[Alpha]+\[Epsilon] (v+eN Sin[v])/l^3) ;
cos\[Kappa]1 = x2+(x3-x2)JacobiSN[\[CapitalUpsilon],\[Beta]^2]^2;
fcos\[Kappa]1[t_]:= cos\[Kappa]1/.findu[t];

(*L vector solution*)
lvecJframe = EulMat . lvec;
l\[Xi]10 = sphericalAngles[lvecJframe][[2]]+ \[Pi]/2;
\[Alpha]1 = -((\[Delta]2 ( j+l+s2 \[CapitalSigma]2))/(s1 (\[Delta]1-\[Delta]2)));
\[Alpha]2 = -((\[Delta]2 (-j+l+s2 \[CapitalSigma]2))/(s1 (\[Delta]1-\[Delta]2)));
\[Beta]1 = -(\[Delta]2/(2 s1 (\[Delta]1-\[Delta]2)))(l^2 \[Delta]1+j^2 \[Delta]2+s1 (\[Delta]1-\[Delta]2) (s1+s2 \[CapitalSigma]1)+l s2 \[Delta]1 \[CapitalSigma]2+j(l (\[Delta]1+\[Delta]2)+s2 \[Delta]2 \[CapitalSigma]2));
\[Beta]2 = -(\[Delta]2/(2 s1 (\[Delta]1-\[Delta]2)))(l^2 \[Delta]1+j^2 \[Delta]2+s1 (\[Delta]1-\[Delta]2) (s1+s2 \[CapitalSigma]1)+l s2 \[Delta]1 \[CapitalSigma]2-j(l (\[Delta]1+\[Delta]2)+s2 \[Delta]2 \[CapitalSigma]2));
solPiece1 = (2/(A(x1-x2))^(1/2)) (( \[Beta]1  EllipticPi[(x2-x3)/(\[Alpha]1+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]1+x2)-( \[Beta]2 EllipticPi[(x2-x3)/(\[Alpha]2+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]2+x2)); 
solPiece2 = l\[Xi]10 -(solPiece1/.u->u0 );
l\[Xi]1sol = solPiece1 + solPiece2;
cos\[Kappa]2[t_]:= \[CapitalSigma]2 - \[Delta]1 s1/(\[Delta]2 s2) fcos\[Kappa]1[t];
lvecAzimuthAngle[t_]:= -\[Pi]/2+l\[Xi]1sol/.findu[t];
lvecPolarAngle[t_]:= ArcCos[(l^2+ l s1  fcos\[Kappa]1[t] + l s2  cos\[Kappa]2[t]  )/(l j)]; 
Lsol[t_]:= Module[{lx,ly,lz},
{lx,ly,lz} = vectorComponents[{l, lvecPolarAngle[t],lvecAzimuthAngle[t]}];
Return[ \[Mu] G M (Inverse[EulMat] . {lx,ly,lz})];]; (*Non reduced vector *)

(*S1 solution*)
s1vecJframe = EulMat . s1vec;
s1\[Xi]10 = sphericalAngles[s1vecJframe][[2]]+ \[Pi]/2;
\[Alpha]1s1 = \[Delta]2 (-j+s1+s2 \[CapitalSigma]1)/(l \[Delta]1);
\[Alpha]2s1 = \[Delta]2 ( j+s1+s2 \[CapitalSigma]1)/(l \[Delta]1);
\[Beta]1s1 = -1/(2 l \[Delta]1)\[Delta]2( l^2 \[Delta]1+( s1(\[Delta]1-\[Delta]2)+j \[Delta]2)(-j+s1+s2 \[CapitalSigma]1)+l s2 \[Delta]1 \[CapitalSigma]2);
\[Beta]2s1 =  1/(2 l \[Delta]1)\[Delta]2(-l^2 \[Delta]1+(-s1 \[Delta]1+(j+s1)\[Delta]2)( j+s1+s2 \[CapitalSigma]1)-l s2 \[Delta]1 \[CapitalSigma]2);
solPiece1s1 =(2/(A(x1-x2))^(1/2)) (( \[Beta]1s1  EllipticPi[(x2-x3)/(\[Alpha]1s1+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]1s1+x2)-(\[Beta]2s1 EllipticPi[(x2-x3)/(\[Alpha]2s1+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]2s1+x2))   + \[Epsilon] (v+eN Sin[v])/l^3 j \[Delta]2  ;
solPiece2s1 = s1\[Xi]10 -(solPiece1s1/.u->u0 );
s1\[Xi]1sol = solPiece1s1 + solPiece2s1 ;
cos\[Gamma][t_]:= \[CapitalSigma]1 +   l (\[Delta]1-\[Delta]2)/(\[Delta]2 s2) fcos\[Kappa]1[t];
s1vecAzimuthAngle[t_]:= -\[Pi]/2+s1\[Xi]1sol/.findu[t];
s1vecPolarAngle[t_]:= ArcCos[(s1^2+ l s1  fcos\[Kappa]1[t] + s1 s2  cos\[Gamma][t]  )/(s1 j)]; 
S1sol[t_]:= Module[{s1xx,s1yy,s1zz},
{s1xx,s1yy,s1zz} = vectorComponents[{s1, s1vecPolarAngle[t],s1vecAzimuthAngle[t]}];
Return[ \[Mu] G M (Inverse[EulMat] . {s1xx ,s1yy,s1zz})];];

(*S2 solution*)
\[CapitalSigma]\[CapitalSigma]1 = N[s2vec . s1vec/(s2 s1) + l/s1 (\[Delta]1-\[Delta]2)/\[Delta]1  s2vec . lvec/(l s2)];
\[CapitalSigma]\[CapitalSigma]2 = N[s1vec . lvec/(l s1) + (\[Delta]2 s2)/(\[Delta]1 s1)  s2vec . lvec/(l s2)];
s2vecJframe = EulMat . s2vec;
s2\[Xi]10 = sphericalAngles[s2vecJframe][[2]]+ \[Pi]/2;
\[Alpha]1s2 = (\[Delta]1 (-j+s2+s1 \[CapitalSigma]\[CapitalSigma]1)/(l \[Delta]2)+\[CapitalSigma]2)(-((\[Delta]2 s2)/(\[Delta]1 s1)));
\[Alpha]2s2 = (\[Delta]1 ( j+s2+s1 \[CapitalSigma]\[CapitalSigma]1)/(l \[Delta]2)+\[CapitalSigma]2)(-((\[Delta]2 s2)/(\[Delta]1 s1)));
\[Beta]1s2 = (-1/(2 l \[Delta]2)\[Delta]1( l^2 \[Delta]2+( s2(\[Delta]2-\[Delta]1)+j \[Delta]1)(-j+s2+s1 \[CapitalSigma]\[CapitalSigma]1)+l s1 \[Delta]2 \[CapitalSigma]\[CapitalSigma]2))(-((\[Delta]2 s2)/(\[Delta]1 s1)));
\[Beta]2s2 = ( 1/(2 l \[Delta]2)\[Delta]1(-l^2 \[Delta]2+(-s2 \[Delta]2+(j+s2)\[Delta]1)( j+s2+s1 \[CapitalSigma]\[CapitalSigma]1)-l s1 \[Delta]2 \[CapitalSigma]\[CapitalSigma]2))(-((\[Delta]2 s2)/(\[Delta]1 s1)));
solPiece1s2 =(2/(A(x1-x2))^(1/2)) (( \[Beta]1s2  EllipticPi[(x2-x3)/(\[Alpha]1s2+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]1s2+x2)-(\[Beta]2s2 EllipticPi[(x2-x3)/(\[Alpha]2s2+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]2s2+x2)) + \[Epsilon] (v+eN Sin[v])/l^3 j \[Delta]1  ;
solPiece2s2 = s2\[Xi]10 -(solPiece1s2/.u->u0 );
s2\[Xi]1sol = solPiece1s2 + solPiece2s2  ;
cos\[Kappa]2[t_]:= \[CapitalSigma]2 - \[Delta]1 s1/(\[Delta]2 s2) fcos\[Kappa]1[t];
cos\[Gamma][t_]:= \[CapitalSigma]1 +   l (\[Delta]1-\[Delta]2)/(\[Delta]2 s2) fcos\[Kappa]1[t];
s2vecAzimuthAngle[t_]:= -\[Pi]/2+s2\[Xi]1sol/.findu[t];
s2vecPolarAngle[t_]:= ArcCos[(s2^2+ l s2  cos\[Kappa]2[t] + s1 s2  cos\[Gamma][t]  )/(s2 j)]; 
S2sol[t_]:= Module[{s2xx,s2yy,s2zz},
{s2xx,s2yy,s2zz} = vectorComponents[{s2, s2vecPolarAngle[t],s2vecAzimuthAngle[t]}];
Return[ \[Mu] G M (Inverse[EulMat] . {s2xx ,s2yy,s2zz})];];

(*R solution*)
LsolJframe[t_]:=Module[{lx,ly,lz},
{lx,ly,lz} = vectorComponents[{l, lvecPolarAngle[t],lvecAzimuthAngle[t]}];
Return[(\[Mu] G M)({lx,ly,lz})];];

JtoLframeEulMat[t_]:=Module[{Lx,Ly,Lz, \[Xi]1L, \[Xi]2L},
{Lx,Ly,Lz} = LsolJframe[t];
{\[Xi]2L, \[Xi]1L}=sphericalAngles[{Lx,Ly,Lz}]+{0,\[Pi]/2};
EulMatJtoL ={ {Cos[\[Xi]1L], Sin[\[Xi]1L],0},{-Sin[\[Xi]1L]Cos[\[Xi]2L], Cos[\[Xi]1L]Cos[\[Xi]2L], Sin[\[Xi]2L]},{Sin[\[Xi]1L]Sin[\[Xi]2L], -Cos[\[Xi]1L]Sin[\[Xi]2L], Cos[\[Xi]2L]}}      ;
Return[EulMatJtoL];];

rLframe = JtoLframeEulMat[0] . EulMat . rvec ;
r\[Phi]0 = sphericalAngles[rLframe][[2]] ;

eN = Sqrt[1-2 l^2 \[CapitalEpsilon]];
\[Beta]eN = eN/(1+(1-eN^2)^(1/2));
vN = u + 2ArcTan[ (\[Beta]eN Sin[u])/(1-\[Beta]eN Cos[u])];

e\[Theta]prime = (1 - 2 l^2 \[CapitalEpsilon] + \[Epsilon](-((12 dotlseff \[CapitalEpsilon])/l^2) + 16 dotlseff \[CapitalEpsilon]^2) + \[Epsilon](-4 (-3+\[Delta]1+\[Delta]2)\[CapitalEpsilon] + l^2 \[CapitalEpsilon]^2 (-15+8(\[Delta]1+\[Delta]2)+\[Nu])))^(1/2);
\[Beta]e\[Theta]prime = e\[Theta]prime/(1+(1-e\[Theta]prime^2)^(1/2));
v\[Theta]prime = u + 2ArcTan[ (\[Beta]e\[Theta]prime Sin[u])/(1-\[Beta]e\[Theta]prime Cos[u])];

solPiece1r = (2/(A(x1-x2))^(1/2)) (( \[Beta]1 EllipticPi[(x2-x3)/(\[Alpha]1+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]1+x2)+( \[Beta]2 EllipticPi[(x2-x3)/(\[Alpha]2+x2),JacobiAmplitude[\[CapitalUpsilon],\[Beta]^2],\[Beta]^2])/(\[Alpha]2+x2));
solPiece2r = v\[Theta]prime - (vN \[Epsilon])/l^4 (3 dotlseff+l^2 (-3 + \[Delta]1 + \[Delta]2 + l^2 \[CapitalEpsilon] + h l^2 (1-3\[Nu]) - 3l^2 \[CapitalEpsilon] \[Nu]));
solPiece3r = solPiece1r + solPiece2r;
solPiece4r = r\[Phi]0 - (solPiece3r/.u->u0);
r\[Phi]solLframe = solPiece3r + solPiece4r;

rvecAzimuthAngleLframe[t_]:=r\[Phi]solLframe/.findu[t];
rMagnitude[t_] := ar (1-er Cos[u])/.findu[t];

Rsol[t_]:=Module[{rx,ry,rz},
{rx,ry,rz} = vectorComponents[{rMagnitude[t], \[Pi]/2,rvecAzimuthAngleLframe[t]}];
Return[ G M(Inverse[EulMat] . Inverse[JtoLframeEulMat[t]] . {rx,ry,rz}) ]; ];

(* P solution*)
pdn[t_] :=  ((-2\[CapitalEpsilon] + \[Epsilon](1- 3 \[Nu])\[CapitalEpsilon]^2) + (2(1 - \[Epsilon] (4-\[Nu])\[CapitalEpsilon]))/rMagnitude[t] + (-l^2 + \[Epsilon](6+\[Nu]))/rMagnitude[t]^2 + (-\[Epsilon] \[Nu] (l^2+ 2dotlseff/\[Nu]) )/rMagnitude[t]^3  )^(1/2) ;

pdninit = rvec . pvec/r; periodpdn = \[Pi]/nunscaled;
first0 = t0; (*Have to define first0 otherwise it does not work, i dont know why*)

If [first0 <=0, first0 = first0 + \[Pi]/nunscaled, first0=first0];
findpdnsign[t_]:=  QuotientRemainder[QuotientRemainder[t-first0,periodpdn][[1]],2][[2]];

pMagnitude[t_]:=( pdn[t]^2+l^2/rMagnitude[t]^2)^(1/2) ;
\[Phi]OffsetPvecNIFTemp[t_] :=  ArcSin[l/(rMagnitude[t]pMagnitude[t])] ;

\[Phi]OffsetPvecNIF[t_]:=Module[{tempVar} ,
If[(pdninit > 0 && findpdnsign[t] == 1) || (pdninit < 0 && findpdnsign[t] == 0),
    tempVar = \[Phi]OffsetPvecNIFTemp[t],
    tempVar = \[Phi]OffsetPvecNIFTemp[t] + 2 (\[Pi]/2 - \[Phi]OffsetPvecNIFTemp[t]);];
Return[tempVar];];

pvecAzimuthAngleLframe[t_]:= rvecAzimuthAngleLframe[t]+\[Phi]OffsetPvecNIF[t];
Psol[t_]:=Module[{px,py,pz},
{px,py,pz} = vectorComponents[{pMagnitude[t], \[Pi]/2,pvecAzimuthAngleLframe[t]}];
Return[ \[Mu](Inverse[EulMat] . Inverse[JtoLframeEulMat[t]] . {px,py,pz}) ]; ];
            
                finalvec=Re[{Rsol[\[Lambda]max-\[Lambda]0],Psol[\[Lambda]max-\[Lambda]0],
                S1sol[\[Lambda]max-\[Lambda]0],S2sol[\[Lambda]max-\[Lambda]0]}]//N;
                (*Print["The final state is"];*)
                Return[finalvec];

                 ]
  End[]

  EndPackage[]
  
  

