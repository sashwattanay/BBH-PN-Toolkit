(* ::Package:: *)

BeginPackage[ "BBHpnToolkit`Hfl`"]

     Hflow::usage = 
	"Hflow implements flow induced by H in phase space "

  Begin[ "`Private`"]

    Hflow15PN[G_,m1_, m2_, R0_,P0_,S10_, S20_,\[Lambda]max_,\[Epsilon]_]:=
               Module[{c,\[Lambda]0,M,\[Mu],\[Nu],rvec,r,pvec,p,lvec,l,s1vec,s2vec,s1,s2,jvec,j,\[Delta]1,\[Delta]2,dotlseff,\[CapitalEpsilon]0,h0,\[CapitalEpsilon],h,sphericalAngles,vectorComponents,
               \[Xi]1,\[Xi]2,EulMat,\[CapitalSigma]1,\[CapitalSigma]2,A,cubiceq,rtsol,x,x1,x2,x3,eN,\[Beta]eN,arN,nN,ar,et,er,n,nunscaled,drMagnitudedtInit,Absu0,Cosu,u0,t0,findu,vN,
               u,cos\[Kappa]10,\[Beta],sign,\[Alpha],\[CapitalUpsilon]N,cos\[Kappa]1,fcos\[Kappa]1,lvecJframe,l\[Xi]10,\[Alpha]1,\[Alpha]2,\[Beta]1,\[Beta]2,solPiece1,solPiece2,l\[Xi]1sol,fcos\[Kappa]2,lvecAzimuthAngle,lvecPolarAngle,
               Lsol,s1vecJframe,s1\[Xi]10,\[Alpha]1s1,\[Alpha]2s1,\[Beta]1s1,\[Beta]2s1,solPiece1s1,solPiece2s1,s1\[Xi]1sol,fcos\[Gamma],s1vecAzimuthAngle,s1vecPolarAngle,S1sol,\[CapitalSigma]\[CapitalSigma]1,\[CapitalSigma]\[CapitalSigma]2,
               s2vecJframe,s2\[Xi]10,\[Alpha]1s2,\[Alpha]2s2,\[Beta]1s2,\[Beta]2s2,solPiece1s2,solPiece2s2,s2\[Xi]1sol,s2vecAzimuthAngle,s2vecPolarAngle,S2sol,LsolJframe,JtoLframeEulMat,
               EulMatJtoL,rLframe,r\[Phi]0,ver,\[Beta]er,solPiece1r,solPiece2r,solPiece3r,solPiece4r,r\[Phi]solLframe,rvecAzimuthAngleLframe,rMagnitude,Rsol,pdninit,
               periodpdn,solr,A0,A1,B0,B1,C0,C1,D1,solrr,rext,rr,distances,closestIndices,solr1,solr2,artilde,ertilde,rMagnitudetilde,pdntilde,first0,
               findpdnsign,pMagnitudetilde,\[Phi]OffsetPvecNIFTemptilde,\[Phi]OffsetPvecNIFtilde,pvecAzimuthAngleLframetilde,Psol,finalvec},
 
(*Last Update 24/04/25*)

(*Set up all the basic parameters*)
c = 1/Sqrt[\[Epsilon]];   
\[Lambda]0=0; (*set initial time to 0 such that \[Lambda]max is the flow amount*) 
M = m1+m2 ;  \[Mu] = m1 m2/M ; \[Nu]=\[Mu]/M;
rvec = R0/(G M);r = Norm[rvec];pvec = P0/\[Mu];p = Norm[pvec];
lvec = Cross[rvec,pvec];l = Norm[lvec] ;
s1vec = S10/(\[Mu] G M); s2vec = S20/(\[Mu] G M);s1 = Norm[s1vec];s2 = Norm[s2vec];
jvec = lvec + s1vec + s2vec  ;  j = Norm[jvec] ;
\[Delta]1 = 2 m1 m2/M^2 (1+3/4 m2/m1);\[Delta]2 = 2 m2 m1/M^2 (1+3/4 m1/m2);dotlseff = lvec . (\[Delta]1  s1vec  +  \[Delta]2 s2vec);(*dotlseff is a constant*)
\[CapitalEpsilon]0 = -h0;h0= p^2/2 - 1/r ;
\[CapitalEpsilon] = -h;h = p^2/2 - 1/r + \[Epsilon]/r^3 dotlseff + \[Epsilon](1/8 (3\[Nu]-1)p^4+1/(2r^2)-1/(2r) ((3+\[Nu])p^2+\[Nu] Dot[rvec/r,pvec]^2));

(*Define functions from Cartesian to Spherical*)
sphericalAngles[V_]:=Module[{polarAngJ,azimuthAngJ},
	polarAngJ = ToSphericalCoordinates[Re[V]][[2]];azimuthAngJ =  ToSphericalCoordinates[Re[V]][[3]];
	Return[{polarAngJ,azimuthAngJ}] ; ];

vectorComponents[{V_, \[Theta]_, \[Phi]_}]:= Module[{Vx,Vy,Vz},
	Vz = V Cos[\[Theta]];Vx=V Sin[\[Theta]]Cos[\[Phi]];  Vy=V Sin[\[Theta]]Sin[\[Phi]];
	Return[{Vx, Vy, Vz}] ; ];

(*Define Euler Matrix from IF (J) to NIF (L)*)
{\[Xi]2, \[Xi]1}=sphericalAngles[jvec]+{0,\[Pi]/2};
EulMat ={{Cos[\[Xi]1], Sin[\[Xi]1],0},{-Sin[\[Xi]1]Cos[\[Xi]2], Cos[\[Xi]1]Cos[\[Xi]2], Sin[\[Xi]2]},{Sin[\[Xi]1]Sin[\[Xi]2], -Cos[\[Xi]1]Sin[\[Xi]2], Cos[\[Xi]2]}};

(*Define new constants via l . (s1 x s2)*)
\[CapitalSigma]1 = N[s1vec . s2vec/(s1 s2) - l/s2 (\[Delta]1-\[Delta]2)/\[Delta]2  s1vec . lvec/(l s1)];  
\[CapitalSigma]2 = N[s2vec . lvec/(l s2) + (\[Delta]1 s1)/(\[Delta]2 s2)  s1vec . lvec/(l s1)];     

(*Define Cubic polynome and its roots xi, x2<x3<x1*)
A = 2 l s1 \[Delta]1 (\[Delta]2-\[Delta]1);cubiceq = A x^3 -(l^2 (\[Delta]1-\[Delta]2)^2+2 \[Delta]2 l \[CapitalSigma]2 s2(\[Delta]2-\[Delta]1)+s1^2 \[Delta]1^2+2 \[Delta]1 \[Delta]2 \[CapitalSigma]1 s1 s2 +\[Delta]2^2  s2^2  )x^2+ (2 \[Delta]2 s2 (l \[CapitalSigma]1(\[Delta]2-\[Delta]1) +\[CapitalSigma]2(\[Delta]1 s1+ \[Delta]2 \[CapitalSigma]1 s2) )) x+-s2^2 \[Delta]2^2 (\[CapitalSigma]1^2+\[CapitalSigma]2^2-1);(*=A(x -x1)(x -x2)(x -x3)*)
rtsol = Solve[cubiceq==0,x];
x2 = x/.rtsol[[1]]; 
x3 = x/.rtsol[[2]];
x1 = x/.rtsol[[3]];
{x2, x3, x1} = Sort[{x2, x3, x1}];

(*qkp OPN, subscript N is for Newtonian order*)
eN = Sqrt[1-2 l^2 \[CapitalEpsilon]0];
\[Beta]eN = eN/(1+(1-eN^2)^(1/2));
arN = 1/(2 \[CapitalEpsilon]0);
nN = (2\[CapitalEpsilon]0)^(3/2);
(*qkp OPN + 1PN + 1.5PN : Samanta Tanay parametrization*)
ar = 1/(2\[CapitalEpsilon])+ (dotlseff \[Epsilon])/l^2+\[Epsilon]/4 (-7+\[Nu]);
er = (1 - 2 l^2 \[CapitalEpsilon] + \[Epsilon](-((8 dotlseff \[CapitalEpsilon])/l^2) + 8dotlseff \[CapitalEpsilon]^2) + \[Epsilon](12 \[CapitalEpsilon]-15 l^2 \[CapitalEpsilon]^2-2 \[CapitalEpsilon] \[Nu]+5 l^2 \[CapitalEpsilon]^2 \[Nu]))^(1/2);et = (1 - 2 l^2 \[CapitalEpsilon] - (4 dotlseff \[Epsilon] \[CapitalEpsilon])/l^2 + \[Epsilon](l^2 \[CapitalEpsilon]^2 (17 - 7\[Nu]) + 4\[CapitalEpsilon](\[Nu]-1)))^(1/2);n = ((2\[CapitalEpsilon])^(3/2)+ \[Epsilon] \[CapitalEpsilon]^(5/2) (-(15/Sqrt[2])+\[Nu]/Sqrt[2]));nunscaled = n/(G M);  

(*COS\[CapitalKappa]1 SOLUTION*)
drMagnitudedtInit = P0 . R0(2+\[Epsilon] (-1+3 \[Nu]) Norm[P0]^2-(2 \[Epsilon] (3+2 \[Nu]))/Norm[R0]);
Absu0 = ArcCos[Cosu/.(Solve[r == ar (1-er Cosu),Cosu]//N)[[1]]];
If[drMagnitudedtInit  > 0,u0 = Absu0;,u0 = -Absu0;];
t0 = -(t0/.Solve[ nunscaled t0 == u0 - et Sin[u0], t0] [[1]]);  (*This is t0 in n(t-t0) = u - et sin[u]*)
findu[t_]:=FindRoot[nunscaled (t-t0) == u - et Sin[u],{u, nunscaled (t-t0), nunscaled (t-t0) - 2 \[Pi], nunscaled (t-t0)+2 \[Pi]}, PrecisionGoal->50]; 
(*At t=0 we have u=u0     and      at t=t0 we have u=0*)

vN = u + 2ArcTan[ (\[Beta]eN Sin[u])/(1-\[Beta]eN Cos[u])];   (* PN extension of v= 2 ArcTan[Sqrt[(1+e)/(1-e)] Tan[u/2]] but without the ArcTan issues*)
cos\[Kappa]10 = Dot[lvec,s1vec]/(Norm[lvec]Norm[s1vec]) (*Initial value of cos\[Kappa]1 at t=0*);
\[Beta] = ((x3-x2)/(x1-x2))^(1/2);
If [ \[Delta]2/(c^2 r^3 s1 l) lvec . Cross[s1vec, s2vec] > 0, sign =1, sign =-1 ] ;
\[Alpha] =  sign 2/Sqrt[A(x1-x2)] EllipticF[ArcSin[Sqrt[(cos\[Kappa]10-x2)/(x3-x2)]],\[Beta]^2]- (+\[Epsilon] (vN+eN Sin[vN])/l^3/.findu[0]);(*Have to remove the value at t=0, so that the inputs match with cos\[Kappa]10 *)
\[CapitalUpsilon]N = 1/2 Sqrt[A(x1-x2)](\[Alpha]+\[Epsilon] (vN+eN Sin[vN])/l^3) ;
cos\[Kappa]1 = x2+(x3-x2)JacobiSN[\[CapitalUpsilon]N,\[Beta]^2]^2;
fcos\[Kappa]1[t_]:= cos\[Kappa]1/.findu[t];

(*L VECTOR SOLUTION*)
lvecJframe = EulMat . lvec;
l\[Xi]10 = sphericalAngles[lvecJframe][[2]]+ \[Pi]/2;

\[Alpha]1 = -((\[Delta]2 ( j+l+s2 \[CapitalSigma]2))/(s1 (\[Delta]1-\[Delta]2)));
\[Alpha]2 = -((\[Delta]2 (-j+l+s2 \[CapitalSigma]2))/(s1 (\[Delta]1-\[Delta]2)));
\[Beta]1 = -(\[Delta]2/(2 s1 (\[Delta]1-\[Delta]2)))(l^2 \[Delta]1+j^2 \[Delta]2+s1 (\[Delta]1-\[Delta]2) (s1+s2 \[CapitalSigma]1)+l s2 \[Delta]1 \[CapitalSigma]2+j(l (\[Delta]1+\[Delta]2)+s2 \[Delta]2 \[CapitalSigma]2));
\[Beta]2 = -(\[Delta]2/(2 s1 (\[Delta]1-\[Delta]2)))(l^2 \[Delta]1+j^2 \[Delta]2+s1 (\[Delta]1-\[Delta]2) (s1+s2 \[CapitalSigma]1)+l s2 \[Delta]1 \[CapitalSigma]2-j(l (\[Delta]1+\[Delta]2)+s2 \[Delta]2 \[CapitalSigma]2));

solPiece1 = (2/(A(x1-x2))^(1/2)) (( 
		\[Beta]1  EllipticPi[(x2-x3)/(\[Alpha]1+x2),JacobiAmplitude[\[CapitalUpsilon]N,\[Beta]^2],\[Beta]^2])/(\[Alpha]1+x2)
		-\[Beta]2 EllipticPi[(x2-x3)/(\[Alpha]2+x2),JacobiAmplitude[\[CapitalUpsilon]N,\[Beta]^2],\[Beta]^2]/(\[Alpha]2+x2)); 

solPiece2 = l\[Xi]10 -(solPiece1/.u->u0 );
l\[Xi]1sol = solPiece1 + solPiece2;
fcos\[Kappa]2[t_]:= \[CapitalSigma]2 - \[Delta]1 s1/(\[Delta]2 s2) fcos\[Kappa]1[t];
lvecAzimuthAngle[t_]:= -\[Pi]/2+l\[Xi]1sol/.findu[t];
lvecPolarAngle[t_]:= ArcCos[(l^2+ l s1  fcos\[Kappa]1[t] + l s2  fcos\[Kappa]2[t]  )/(l j)]; 
Lsol[t_]:= Module[{lx,ly,lz},
	{lx,ly,lz} = vectorComponents[{l, lvecPolarAngle[t],lvecAzimuthAngle[t]}];
	Return[ \[Mu] G M (Inverse[EulMat] . {lx,ly,lz})];]; (*Non reduced vector *)



(*S1 VECTOR SOLUTION*)
s1vecJframe = EulMat . s1vec;
s1\[Xi]10 = sphericalAngles[s1vecJframe][[2]]+ \[Pi]/2;

\[Alpha]1s1 = \[Delta]2 (-j+s1+s2 \[CapitalSigma]1)/(l \[Delta]1);
\[Alpha]2s1 = \[Delta]2 ( j+s1+s2 \[CapitalSigma]1)/(l \[Delta]1);
\[Beta]1s1 = -1/(2 l \[Delta]1)\[Delta]2( l^2 \[Delta]1+( s1(\[Delta]1-\[Delta]2)+j \[Delta]2)(-j+s1+s2 \[CapitalSigma]1)+l s2 \[Delta]1 \[CapitalSigma]2);
\[Beta]2s1 =  1/(2 l \[Delta]1)\[Delta]2(-l^2 \[Delta]1+(-s1 \[Delta]1+(j+s1)\[Delta]2)( j+s1+s2 \[CapitalSigma]1)-l s2 \[Delta]1 \[CapitalSigma]2);

solPiece1s1 =(2/(A(x1-x2))^(1/2)) (
     \[Beta]1s1 EllipticPi[(x2-x3)/(\[Alpha]1s1+x2),JacobiAmplitude[\[CapitalUpsilon]N,\[Beta]^2],\[Beta]^2]/(\[Alpha]1s1+x2)
	-\[Beta]2s1 EllipticPi[(x2-x3)/(\[Alpha]2s1+x2),JacobiAmplitude[\[CapitalUpsilon]N,\[Beta]^2],\[Beta]^2]/(\[Alpha]2s1+x2)
	+j \[Delta]2 \[CapitalUpsilon]N );
     
solPiece2s1 = s1\[Xi]10 -(solPiece1s1/.u->u0 );
s1\[Xi]1sol = solPiece1s1 + solPiece2s1 ;
fcos\[Gamma][t_]:= \[CapitalSigma]1 +   l (\[Delta]1-\[Delta]2)/(\[Delta]2 s2) fcos\[Kappa]1[t];
s1vecAzimuthAngle[t_]:= -\[Pi]/2+s1\[Xi]1sol/.findu[t];
s1vecPolarAngle[t_]:= ArcCos[(s1^2+ l s1  fcos\[Kappa]1[t] + s1 s2  fcos\[Gamma][t]  )/(s1 j)]; 
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

solPiece1s2 =(2/(A(x1-x2))^(1/2)) (
     \[Beta]1s2  EllipticPi[(x2-x3)/(\[Alpha]1s2+x2),JacobiAmplitude[\[CapitalUpsilon]N,\[Beta]^2],\[Beta]^2]/(\[Alpha]1s2+x2)
	-\[Beta]2s2 EllipticPi[(x2-x3)/(\[Alpha]2s2+x2),JacobiAmplitude[\[CapitalUpsilon]N,\[Beta]^2],\[Beta]^2]/(\[Alpha]2s2+x2)
	+j \[Delta]1 \[CapitalUpsilon]N);

solPiece2s2 = s2\[Xi]10 -(solPiece1s2/.u->u0 );
s2\[Xi]1sol = solPiece1s2 + solPiece2s2  ;
s2vecAzimuthAngle[t_]:= -\[Pi]/2+s2\[Xi]1sol/.findu[t];
s2vecPolarAngle[t_]:= ArcCos[(s2^2+ l s2  fcos\[Kappa]2[t] + s1 s2  fcos\[Gamma][t]  )/(s2 j)]; 
S2sol[t_]:= Module[{s2xx,s2yy,s2zz},
	{s2xx,s2yy,s2zz} = vectorComponents[{s2, s2vecPolarAngle[t],s2vecAzimuthAngle[t]}];
	Return[ \[Mu] G M (Inverse[EulMat] . {s2xx ,s2yy,s2zz})];];

(*R VECTOR SOLUTION*)
LsolJframe[t_]:=Module[{lx,ly,lz},
	{lx,ly,lz} = vectorComponents[{l, lvecPolarAngle[t],lvecAzimuthAngle[t]}];
	Return[(\[Mu] G M)({lx,ly,lz})];];
JtoLframeEulMat[t_]:=Module[{Lx,Ly,Lz, \[Xi]1L, \[Xi]2L},
	{Lx,Ly,Lz} = LsolJframe[t];{\[Xi]2L, \[Xi]1L}=sphericalAngles[{Lx,Ly,Lz}]+{0,\[Pi]/2};
	EulMatJtoL ={ {Cos[\[Xi]1L], Sin[\[Xi]1L],0},{-Sin[\[Xi]1L]Cos[\[Xi]2L], Cos[\[Xi]1L]Cos[\[Xi]2L], Sin[\[Xi]2L]},{Sin[\[Xi]1L]Sin[\[Xi]2L], -Cos[\[Xi]1L]Sin[\[Xi]2L], Cos[\[Xi]2L]}};
	Return[EulMatJtoL];];

rLframe = JtoLframeEulMat[0] . EulMat . rvec ;
r\[Phi]0 = sphericalAngles[rLframe][[2]] ;

(*For l/r^2 : need to go to 1.5PN qkp*)
ver = u + 2ArcTan[ (\[Beta]er Sin[u])/(1-\[Beta]er Cos[u])]; 
\[Beta]er = er/(1+(1-er^2)^(1/2));

solPiece1r = (2/(A(x1-x2))^(1/2)) (
       \[Beta]1 EllipticPi[(x2-x3)/(\[Alpha]1+x2),JacobiAmplitude[\[CapitalUpsilon]N,\[Beta]^2],\[Beta]^2]/(\[Alpha]1+x2)
	 + \[Beta]2 EllipticPi[(x2-x3)/(\[Alpha]2+x2),JacobiAmplitude[\[CapitalUpsilon]N,\[Beta]^2],\[Beta]^2]/(\[Alpha]2+x2));

solPiece2r = -\[Epsilon] l (4+\[Delta]1+\[Delta]2-2\[Nu]) (vN+eN Sin[vN])/l^3  -\[Epsilon] l \[CapitalEpsilon](\[Minus]1+3\[Nu]) vN/(arN^2 nN Sqrt[1-eN^2])+ (l (-(((-1+er et) ver)/(1-er^2)^(3/2))+((er-et) Sin[u])/((-1+er^2) (-1+er Cos[u]))))/(ar^2 n);
solPiece3r = solPiece1r + solPiece2r;
solPiece4r = r\[Phi]0 - (solPiece3r/.u->u0);
r\[Phi]solLframe = solPiece3r + solPiece4r;
rvecAzimuthAngleLframe[t_]:=r\[Phi]solLframe/.findu[t];
rMagnitude[t_] := ar (1-er Cos[u])/.findu[t];
Rsol[t_]:=Module[{rx,ry,rz},
	{rx,ry,rz} = vectorComponents[{rMagnitude[t], \[Pi]/2,rvecAzimuthAngleLframe[t]}];
	Return[ G M(Inverse[EulMat] . Inverse[JtoLframeEulMat[t]] . {rx,ry,rz}) ]; ];


(*P VECTOR SOLUTION*)
pdninit = rvec . pvec/r; 
periodpdn = \[Pi]/nunscaled;

solr = Solve[Sqrt[(A0+\[Epsilon] A1)rr^3+(B0+\[Epsilon] B1)rr^2+(C0+\[Epsilon] C1)rr+(\[Epsilon] D1)]==0,rr];
A0=-2 \[CapitalEpsilon];A1=\[CapitalEpsilon]^2-3 \[CapitalEpsilon]^2 \[Nu];B0=2;B1=-8 \[CapitalEpsilon]+2 \[CapitalEpsilon] \[Nu];C0=-l^2;C1=6+\[Nu];D1=-2 dotlseff-l^2 \[Nu];

solrr=rr/.solr//Re; rext = {ar(1-er),ar(1+er)};
distances=Table[EuclideanDistance[rext[[i]],solrr[[j]]],{i,Length[rext]},{j,Length[solrr]}];
closestIndices=Table[First@First@Position[distances[[i]],Min[distances[[i]]]],{i,Length[rext]}];solr1= solrr[[closestIndices[[1]]]]; 
solr2= solrr[[closestIndices[[2]]]];

artilde=(solr1+solr2)/2;
ertilde=(-solr1+solr2)/(solr1+solr2);

rMagnitudetilde[t_] := artilde (1-ertilde Cos[u])/.findu[t];
pdntilde[t_] :=  ((-2\[CapitalEpsilon] + \[Epsilon](1- 3 \[Nu])\[CapitalEpsilon]^2) + (2(1 - \[Epsilon] (4-\[Nu])\[CapitalEpsilon]))/rMagnitudetilde[t] + (-l^2 + \[Epsilon](6+\[Nu]))/rMagnitudetilde[t]^2 + (-\[Epsilon] \[Nu] (l^2+ 2dotlseff/\[Nu]) )/rMagnitudetilde[t]^3  )^(1/2) ; 

first0 = t0; (*Have to define first0, can not use only t0 otherwise it does not work, i dont know why*)
If [first0 <=0, first0 = first0 + \[Pi]/nunscaled, first0=first0]; (*take positive zero between [-\[Pi]/n , +\[Pi]/n]*)
findpdnsign[t_]:=  QuotientRemainder[QuotientRemainder[t-first0,periodpdn][[1]],2][[2]]; (*Calculate number of complete half periods*)

pMagnitudetilde[t_]:=( pdntilde[t]^2+l^2/rMagnitudetilde[t]^2)^(1/2) ;
\[Phi]OffsetPvecNIFTemptilde[t_] :=  ArcSin[l/(rMagnitudetilde[t]pMagnitudetilde[t])] ;
\[Phi]OffsetPvecNIFtilde[t_]:=Module[{tempVar} ,
	If[(pdninit > 0 && findpdnsign[t] == 1) || (pdninit < 0 && findpdnsign[t] == 0),tempVar = \[Phi]OffsetPvecNIFTemptilde[t],tempVar = \[Phi]OffsetPvecNIFTemptilde[t] + 2 (\[Pi]/2 - \[Phi]OffsetPvecNIFTemptilde[t]);];
	Return[tempVar];];

pvecAzimuthAngleLframetilde[t_]:= rvecAzimuthAngleLframe[t]+\[Phi]OffsetPvecNIFtilde[t];
Psol[t_]:=Module[{px,py,pz},{px,py,pz} = vectorComponents[{pMagnitudetilde[t], \[Pi]/2,pvecAzimuthAngleLframetilde[t]}];
	Return[\[Mu](Inverse[EulMat] . Inverse[JtoLframeEulMat[t]] . {px,py,pz}) ]; ];

(*Initial values of R, L, S1 and S2 will match the inputs but not P due to the adjustment of qkp*)
               
finalvec=Re[{Rsol[\[Lambda]max-\[Lambda]0],Psol[\[Lambda]max-\[Lambda]0],S1sol[\[Lambda]max-\[Lambda]0],S2sol[\[Lambda]max-\[Lambda]0]}]//N;
    Return[finalvec];]
  End[]

  EndPackage[]





