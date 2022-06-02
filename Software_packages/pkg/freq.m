(* ::Package:: *)

BeginPackage[ "pkg`freq`"]

     frequency::usage = 
	"frequency is a simple function."

  Begin[ "`Private`"]

    frequency[m1_, m2_, Rinit_,Pinit_,S1init_, S2init_,\[Lambda]_,\[Lambda]0_,\[Epsilon]_]:=Module[{G,c,M, \[Mu], \[Nu], Q1,Q2,Linit,Jinit,RN,PN,Rn,
    Pn,S1N,S2N,S1n,S2n,LN,Ln,JN,Jn,Jz,Seffinit,SeffLN,SeffL,En,H,\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]21,\[CapitalSigma]1,\[CapitalSigma]2,a3,a2,a1,a0,A,p,q,f1,f2,f3,k,B1,B2,
    D1,D2,\[Alpha]1sq,\[Alpha]2sq,B1s1,B2s1,D1s1,D2s1,\[Alpha]1s1sq,\[Alpha]2s1sq,B1s2,B2s2,D1s2,D2s2,\[Alpha]1s2sq,\[Alpha]2s2sq,\[CapitalDelta]\[Lambda]1,\[CapitalDelta]\[Lambda]2,\[CapitalDelta]\[Lambda]3,\[CapitalDelta]\[Lambda]4,\[CapitalDelta]\[Lambda]5,
    CC,J1,J2,J3,J4,J5,Jac,JCmat,JCmatN,Rx,Ry,Rz,Px,Py,Pz,S1x,S1y,S1z,S2x,S2y,S2z,CJmat,freqN},
G=1;     c = 1/Sqrt[\[Epsilon]]   ;
M= m1+m2;
\[Mu] =m1  m2 /(m1+m2);
 \[Nu]= \[Mu]/M;
Q1=(1+3 m2/(4 m1)); Q2=(1+3 m1/(4 m2));


Linit=Cross[Rinit, Pinit];
Jinit=Linit+S1init+S2init (*check Jinit is along z axis*);
RN=Norm[Rinit];
PN=Norm[Pinit];
S1N = Norm[S1init] ;    
S2N = Norm[S2init] ;
LN=Norm[Linit] ;
JN=Norm[Jinit];
Seffinit= Q1 S1init + Q2 S2init;
SeffLN= (Q1 S1init+ Q2 S2init) . Linit;
En=\[Mu] ((PN/\[Mu])^2/2-1/(RN/(G M))) +\[Mu]/c^2 (1/8 (3 \[Nu]-1)(PN/\[Mu])^4 +1/(2 (RN/(G M))^2)-1/(2 (RN/(G M))) ((3+\[Nu])(PN/\[Mu])^2 +\[Nu] ((Rinit/RN) . ( Pinit/\[Mu]))^2)) + ( (2 G)/(c^2 RN^3)  SeffLN)  ;
Print["The energy for the initial data is"];
Print[En];
\[CapitalDelta]1=(1/(Q1-Q2))((1/2)(Jn^2-Ln^2-S1n^2-S2n^2)-SeffL/Q2 );
\[CapitalDelta]2=(1/(Q1-Q2))((1/2)(Jn^2-Ln^2-S1n^2-S2n^2)-SeffL/Q1 );
\[CapitalDelta]21=SeffL/(Q1 Q2 );
\[CapitalSigma]1=((Q1-Q2) \[CapitalDelta]1)/(S1n * S2n);
\[CapitalSigma]2= SeffL/( Q2 * Ln* S2n );a3 =2 Q1 Q2(Q2-Q1);a2=2(\[CapitalDelta]1+\[CapitalDelta]2)(Q1-Q2)Q1 Q2- Ln^2 (Q1-Q2)^2 -Q1^2 S1n^2-Q2^2 S2n^2;
a1=2(Q1^2 S1n^2 \[CapitalDelta]2+Q2^2 S2n^2 \[CapitalDelta]1+Q1 Q2 \[CapitalDelta]1 \[CapitalDelta]2 (Q2-Q1));
a0=Ln^2 S1n^2 S2n^2-Q1^2 S1n^2 \[CapitalDelta]2^2-Q2^2 S2n^2 \[CapitalDelta]1^2;


A=2 Q1 Q2(Q2-Q1);
p=(3 a1 a3 -a2^2)/(3 a3^2); 
q=(2 a2^3-9 a1 a2 a3+27 a0 a3^2)/(27 a3^3);
f1=-a2/(3 a3)+2 Sqrt[-p/3]Cos[1/3 ArcCos[(3 q)/(2p) Sqrt[-3/p]]+(2 \[Pi] )/3];f2=-a2/(3 a3)+2 Sqrt[-p/3]Cos[1/3 ArcCos[(3 q)/(2p) Sqrt[-3/p]]+(2*2 \[Pi] )/3] ;f3=-a2/(3 a3)+2 Sqrt[-p/3]Cos[1/3 ArcCos[(3 q)/(2p) Sqrt[-3/p]]+(3*2 \[Pi] )/3];

k=(f2-f1)/(f3-f1);

B1=(1/2)((SeffL +Ln^2 (Q1+Q2))(Jn+Ln)+Ln(Q1  S1n^2+Q2  S2n^2+(Q1+Q2)(\[CapitalDelta]2 Q1 -\[CapitalDelta]1 Q2)));B2=(1/2)((SeffL +Ln^2 (Q1+Q2))(Jn-Ln)-Ln(Q1  S1n^2+Q2  S2n^2+(Q1+Q2)(\[CapitalDelta]2 Q1 -\[CapitalDelta]1 Q2)));D1=Ln(Ln+Jn)+(\[CapitalDelta]2 Q1 -\[CapitalDelta]1 Q2);D2=Ln(Ln-Jn)+(\[CapitalDelta]2 Q1 -\[CapitalDelta]1 Q2);
\[Alpha]1sq=((-f1+f2) (Q1-Q2))/(D1+f1 (-Q1+Q2)); \[Alpha]2sq=((-f1+f2) (Q1-Q2))/(D2+f1 (-Q1+Q2));



(*subspin space for spin S1*)
B1s1=1/2 (-S1n Q1 (Ln^2 -Jn S1n + S1n^2 +\[CapitalDelta]2 Q1)+(Jn-S1n)^2 S1n Q2-(Jn-2S1n)\[CapitalDelta]1 Q1 Q2 +(Jn-S1n)\[CapitalDelta]1 Q2^2);
B2s1=1/2 (S1n Q1 (Ln^2 +Jn S1n + S1n^2 +\[CapitalDelta]2 Q1)-(Jn+S1n)^2 S1n Q2-(Jn+2S1n)\[CapitalDelta]1 Q1 Q2 +(Jn+S1n)\[CapitalDelta]1 Q2^2);
D1s1=(S1n-Jn)S1n - \[CapitalDelta]1 Q2;
D2s1=(S1n+Jn)S1n - \[CapitalDelta]1 Q2;
\[Alpha]1s1sq=((-Q1)(f2-f1) )/(D1s1+f1 Q1); \[Alpha]2s1sq=((-Q1)(f2-f1))/(D2s1+f1 Q1);


(*subspin space for spin S2*)B1s2 =  1/2 (-Jn^2 S2n Q1+2 Jn S2n^2 Q1-S2n^3 Q1+Jn \[CapitalDelta]2 Q1^2-S2n \[CapitalDelta]2 Q1^2+Ln^2 S2n Q2-Jn S2n^2 Q2+S2n^3 Q2-Jn \[CapitalDelta]2 Q1 Q2+2 S2n \[CapitalDelta]2 Q1 Q2-S2n \[CapitalDelta]1 Q2^2) ;
B2s2=1/2 (Jn^2 S2n Q1+2 Jn S2n^2 Q1+S2n^3 Q1+Jn \[CapitalDelta]2 Q1^2+S2n  \[CapitalDelta]2 Q1^2-Ln^2 S2n Q2-Jn S2n^2 Q2-S2n^3 Q2-Jn \[CapitalDelta]2 Q1 Q2-2 S2n \[CapitalDelta]2 Q1 Q2+S2n \[CapitalDelta]1 Q2^2)  ;
D1s2= (Jn S2n-S2n^2-\[CapitalDelta]2 Q1);
D2s2= (-Jn S2n-S2n^2-\[CapitalDelta]2 Q1);
\[Alpha]1s2sq = ((-Q2)(f2-f1))/(D1s2+f1 Q2)  ;
\[Alpha]2s2sq = ((-Q2)(f2-f1))/(D2s2+f1 Q2)  ;


(*\[CapitalDelta]\[Lambda]s are the flow amount along each commuting constant *)

\[CapitalDelta]\[Lambda]1=(4  EllipticK[k])/Sqrt[A(f3-f1)]; (* Eq 66*)


\[CapitalDelta]\[Lambda]2=(-1/(2 Jn))(4/Sqrt[A (f3-f1)] ((B1 EllipticPi[\[Alpha]1sq,k])/(D1-f1 (Q1-Q2))+(B2 EllipticPi[\[Alpha]2sq,k])/(D2-f1 (Q1-Q2)))  )    ;  (* Eq 90*)
\[CapitalDelta]\[Lambda]3=(-1/(2 Ln)) (4 /Sqrt[A (f3-f1)] ((B1 EllipticPi[\[Alpha]1sq,k])/(D1-f1 (Q1-Q2))-(B2 EllipticPi[\[Alpha]2sq,k])/(D2-f1 (Q1-Q2))) -(Ln^2 (Q1+Q2)+SeffL+Q1 Q2 (\[CapitalDelta]1-\[CapitalDelta]2)) /Ln \[CapitalDelta]\[Lambda]1); (*Eq 100*)




\[CapitalDelta]\[Lambda]4=(-1/(2 S1n))(4/Sqrt[A (f3-f1)] (-((B1s1 EllipticPi[\[Alpha]1s1sq,k])/(D1s1+f1 Q1))+(B2s1 EllipticPi[\[Alpha]2s1sq,k])/(D2s1+f1 Q1))+S1n (Q2-Q1) \[CapitalDelta]\[Lambda]1 ) ;(*Eq116*)
\[CapitalDelta]\[Lambda]5=(-1/(2 S2n))(4/Sqrt[A (f3-f1)] (-((B1s2 EllipticPi[\[Alpha]1s1sq,k])/(D1s2+f1 Q2))+(B2s2 EllipticPi[\[Alpha]2s2sq,k])/(D2s2+f1 Q2))+S2n (Q1-Q2) \[CapitalDelta]\[Lambda]1 ) ;



CC={Jn,Jz,Ln,H, SeffL}; (* Commuting constants*)
J1=Jn; J2=Jz; J3=Ln;J4=-Ln +(G M \[Mu]^(3/2))/Sqrt[-2 H]+( G M)/c^2 ((3 G M \[Mu]^2)/Ln+ (Sqrt[-H] \[Mu]^(1/2) (\[Nu]-15))/Sqrt[32]-(2 G \[Mu]^3)/Ln^3 SeffL); J5=1/\[Pi] (SeffL  \[CapitalDelta]\[Lambda]1 + Jn^2 \[CapitalDelta]\[Lambda]2 + Ln^2  \[CapitalDelta]\[Lambda]3  +S1n^2 \[CapitalDelta]\[Lambda]4 +S2n^2 \[CapitalDelta]\[Lambda]5);
Jac={J1,J2,J3,J4,J5};
JCmat= Table[\!\(
\*SubscriptBox[\(\[PartialD]\), \(CC[\([j]\)]\)]\(Jac[\([i]\)]\)\),{i,1,Length[Jac]},{j,1,Length[CC]}];





 Rn=RN; Pn=PN;S1n=S1N;S2n=S2N; Ln=LN;Jn=JN;SeffL=SeffLN; H=En; 

JCmatN=JCmat;
CJmat=Inverse[JCmatN];

freqN= CJmat[[4]];
Print["The frequencies are"];
Print[freqN];
ListPlot[(Tooltip[{Re[#1],Im[#1]}]&)/@freqN,AspectRatio->1]


]

  End[]

  EndPackage[]
