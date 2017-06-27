(* ::Package:: *)

BeginPackage["StructuralColor3D`",{"StructuralColorFramework`"}]

Iteration3D::usage="Iteration3D[color, PhotonicSystem] returns color-matching results given the {X,Y,Z} coordinate of the desired color and a 3D PhotonicSystem."				

PhotonicSystemGlass::usage="PhotonicSystemGlass[n1, n2, sampleThickness] returns a PhotonicsSystem to use in an Iteration3D. n1 is for the particle material and n2 is for the matrix material. The thickness is set to a default of 15000.0 nm."

IterationReset3D::usage="IterationReset3D is used to reset after an iteration is complete."


Begin["Private`"]

defaultEps1=1*10^-3;
defaultEps2=1*10^-16;

fresnelReflection[n1_,n2_,incidentAngle_]:=Module[{theta,ms,rPerp,root,rPar,valTable},
theta=incidentAngle;
ms=(n2/n1);
valTable=Table[
Which[ ms<1&&thetaVal>=ArcSin[ms],
rPerp=1;rPar=1;,True,root=Sqrt[n2^2-(n1*Sin[thetaVal])^2]];rPar=(Abs[(n1*root-n2^2*Cos[thetaVal])/(n1*root+n2^2*Cos[thetaVal])])^2;rPerp=(Abs[(n1*Cos[thetaVal]-root)/(n1*Cos[thetaVal]+root)])^2;{rPar,rPerp},{thetaVal,theta}];
Transpose[valTable]
];

fresnelTransmission[n1_,n2_,incidentAngle_]:=Module[{rPar,rPerp},
{rPar,rPerp}=fresnelReflection[n1,n2,incidentAngle];
{1.0-rPar,1.0-rPerp}];

nStop[x_]:=Round[Abs[x+4.05*x^(1./3.)+2]];

dn1Down[z_,nmx_,nstop_,startVal_]:=Module[{dn},
dn=ConstantArray[0+0*I,nmx+1];
dn[[nmx]]=startVal;
Do[dn[[i]]=(i+1.)/z-1.0/(dn[[i+1]]+(i+1.)/z);, {i, nmx-1,0,-1}];
Table[dn[[i]], {i, 0, nstop-1}]];

lentzDn1[z_,n_,eps1_:defaultEps1,eps2_:defaultEps2]:=Module[{ai,numerator,denominator,product,ratio,ctr,aiVal,xi1,xi2},
ai[i_]:=(-1.)^(i+1)*2.*(n+i-0.5)/z;
numerator=ai[2]+1./ai[1];
denominator=ai[2];
product=ai[1]*numerator/denominator;
ratio =product;
ctr=3;
While[Abs[Re[product]-1]>eps2 \[Or]Abs[Im[product]]>eps2,
aiVal=ai[ctr];numerator=aiVal+1./numerator;denominator=aiVal+1./denominator;
If[Abs[numerator/aiVal]<eps1 \[Or]Abs[denominator/aiVal]<eps1,
xi1=1.+ai[ctr+1]*numerator;xi2=1.+ai[ctr+1]*denominator;ratio=ratio*xi1/xi2;numerator=ai[ctr+2]+numerator/xi1;denominator=ai[ctr+2]+denominator/xi2;ctr=ctr+2;,Nothing];
product=numerator/denominator;ratio=ratio*product;
ctr=ctr+1;];ratio-Floor[n/z]];

riccatiPsiXi[x_,nStop_]:=Module[{psiN,xIn},
psiN=Table[x*SphericalBesselJ[val,1.0*x]+0.0I, {val, 0,nStop}];
xIn=MapThread[#1+#2*I&,{psiN,Table[x*SphericalBesselY[val,1.0*x],{val,0,nStop}]}];
{psiN, xIn}];

scatCoeffs[m_,x_,nstopVal_,eps1_:defaultEps1,eps2_:defaultEps2]:=Module[{dnmx,n,psi,xi,psiShift,xiShift,an,bn},
dnmx=dn1Down[m*x,nstopVal+1,nstopVal+1,lentzDn1[m*x,nstopVal+1,eps1,eps2]];
n=Range[0,nstopVal];
{psi,xi}=riccatiPsiXi[x,nstopVal];
psiShift=Join[{0},psi][[1;;nstopVal+1]];
xiShift=Join[{0},xi][[1;;nstopVal+1]];
an=((dnmx/m+n/x)*psi-psiShift)/((dnmx/m+n/x)*xi-xiShift);
bn=((dnmx*m+n/x)*psi-psiShift)/((dnmx*m+n/x)*xi-xiShift);
{an[[2;;nstopVal+1]],bn[[2;;nstopVal+1]]}];

pisAndTaus[nstopVal_,theta_]:=Module[{mu,legendre0,val,muVal,pis,pishift,n,mus,taus},
mu=Cos[theta];
legendre0=Table[D[LegendreP[val,muVal],muVal]/.{val->value,muVal->mu},{value, 0, nstopVal}];
pis=legendre0[[1;;nstopVal+1]];
pishift=Join[{0},pis][[1;;nstopVal+1]];
n=Range[0,nstopVal];
mus=mu*ConstantArray[1,nstopVal+1];
taus=n*pis*mus-(n+1)*pishift;
{pis[[2;;nstopVal+1]],taus[[2;;nstopVal+1]]}];

amplitudeScatteringMatrix[nStopVal_,prefactor_,coeffs_,theta_]:=Module[{angfuncs,pis,taus,s1,s2},
angfuncs=pisAndTaus[nStopVal,theta];
pis=angfuncs[[1]];
taus=angfuncs[[2]];
s1=Total@(prefactor*(coeffs[[1]]*pis+coeffs[[2]]*taus));
s2=Total@(prefactor*(coeffs[[1]]*taus+coeffs[[2]]*pis));
{s2,s1}];

structureFactor[qd_,phi_]:=Module[{lambda1,lambda2,c},
lambda1=(1+2*phi)^2/(1-phi)^4;
lambda2=-(1+phi/2.)^2/(1-phi)^4;
c=-24*phi*(lambda1*(Sin[qd]-qd*Cos[qd])/qd^3-6*phi*lambda2*(qd^2*Cos[qd]-2*qd*Sin[qd]-2*Cos[qd]+2.0)/qd^4-(phi*lambda1/2.)*(qd^4*Cos[qd]-4*qd^3*Sin[qd]-12*qd^2*Cos[qd]+24*qd*Sin[qd]+24*Cos[qd]-24.0)/qd^6);
1.0/(1-c)];

calcAngDist[m_,x_,angles_]:=
Module[{ipar,iperp,nstopVal,coeffs,n,prefactor,asmat,par,perp},ipar={};iperp={};
nstopVal=nStop[x];
coeffs=scatCoeffs[m,x,nstopVal];
n=Range[0,nstopVal-1]+1.;
prefactor=(2*n+1.)/(n*(n+1.));
Table[asmat=amplitudeScatteringMatrix[nstopVal,prefactor,coeffs,angles[[i]]];
par=Abs[asmat[[1]]]^2;
ipar=Append[ipar,par];
perp=Abs[asmat[[2]]]^2;
iperp=Append[iperp,perp];,
{i, Range[Length[angles]]}];
{ipar,iperp}];

differentialCrossSection[m_,x_,angles_,volumeFraction_]:=Module[{formFactor,fPar,fPerp,qd,s,scatPar,scatPerp},
formFactor=calcAngDist[m,x,angles];
fPar=formFactor[[1]];fPerp=formFactor[[2]];qd=4*x*Sin[angles/2];
s=structureFactor[qd,volumeFraction];
scatPar=s*fPar;
scatPerp=s*fPerp;
{scatPar,scatPerp}];

integrateCrossSection[crossSection_,factor_,angles_]:=Module[{integrand},
integrand=crossSection*factor*Sin[angles];
2*Pi*Total@Table[(integrand[[num]]+integrand[[num+1]])/2*(angles[[num+1]]-angles[[num]]),{num, Range[Length[angles]-1]}]];

photonicGlassReflection[nParticle_,nMatrix_,nMedium_,wavelen_,radius_,volumeFraction_,thickness_:15000.0,thetaMin_:90,thetaMax_:179.9999,incidentAngle_:0,numAngles_:25,smallAngle_:5]:=Module[{nSample,m,x2,k,tMediumSample,rMediumSample,sinAlphaSample,thetaMinRefracted},nSample=nMedium*Sqrt[(2*nMedium^2+nParticle^2+2*volumeFraction*((nParticle^2)-(nMedium^2)))/(2*nMedium^2+nParticle^2-volumeFraction*((nParticle^2)-(nMedium^2)))];
m=nParticle/nSample;
x2=2*Pi*nSample/wavelen*radius;
k=2*Pi*nSample/wavelen;
tMediumSample=fresnelTransmission[nMedium,nSample,{incidentAngle*Pi/180}]//Flatten;
rMediumSample=fresnelReflection[nMedium,nSample,{incidentAngle*Pi/180}]//Flatten;
sinAlphaSample=Sin[Pi-thetaMin*Pi/180]*nMedium/nSample;
Which[sinAlphaSample>=1,thetaMinRefracted=Pi/2.0,True,thetaMinRefracted=Pi-ArcSin[sinAlphaSample]];

angles=Range[thetaMinRefracted,thetaMax*Pi/180,(thetaMax*Pi/180-thetaMinRefracted)/numAngles];
diffCs=differentialCrossSection[m,x2,angles,volumeFraction];
diffCs={ReplacePart[diffCs[[1]], -1->diffCs[[1,-2]]],ReplacePart[diffCs[[2]], -1->diffCs[[2,-2]]]};
transmission=fresnelTransmission[nSample,nMedium,Pi-angles];
transmission={ReplacePart[transmission[[1]], 1->0],ReplacePart[transmission[[2]], 1->0]};
sigmaDetectedPar=integrateCrossSection[diffCs[[1]],transmission[[1]]/k^2,angles];
sigmaDetectedPerp=integrateCrossSection[diffCs[[2]],transmission[[2]]/k^2,angles];
sigmaDetected=(sigmaDetectedPar+sigmaDetectedPerp)/2.0;

angles=Range[0.0+smallAngle*Pi/180,Pi-0.0001,(Pi-smallAngle*Pi/180)/numAngles];
diffCs=differentialCrossSection[m,x2,angles,volumeFraction];
sigmaTotalPar=integrateCrossSection[diffCs[[1]],1.0/k^2,angles];
sigmaTotalPerp=integrateCrossSection[diffCs[[2]],1.0/k^2,angles];
sigmaTotal=(sigmaTotalPar+sigmaTotalPerp)/2.0;

rho=3.0*volumeFraction/(4.0*Pi*radius^3);
factor=1.0-Exp[-rho*sigmaTotal*thickness];reflectedPar=tMediumSample[[1]]*sigmaDetectedPar/sigmaTotal*factor+rMediumSample[[1]];
reflectedPerp=tMediumSample[[2]]*sigmaDetectedPerp/sigmaTotal*factor+rMediumSample[[2]];
(reflectedPar+reflectedPerp)/2.0];

(*calculations*)

PhotonicSystemGlass[nOne_,nTwo_,thickness_:15000.0]:={nOne, nTwo,thickness};
photonicGlassSpectrum[radius_,volumeFraction_,{nOne_, nTwo_, thickness_}]:=Module[{wavelength},Interpolation[Table[{wavelength,photonicGlassReflection[nOne,nTwo,1.,wavelength,radius,volumeFraction,thickness]},{wavelength,400.,700.,10.}]][\[Lambda]]];
expandedMatrix3D[dimension_, resolution_, matrixDimension_:3]:=Module[{matrixRange,potentialMatrix},
matrixRange=Range[-Floor[matrixDimension/2], Floor[matrixDimension/2]];
potentialMatrix=Table[{dimension[[1]]+resolution[[1]]*a, dimension[[2]]+resolution[[2]]*b}, {a, matrixRange}, {b, matrixRange}];
Cases[Flatten[potentialMatrix,1], {a_, b_}/;100<=a<=1000&&0.01<=b<=0.7404]];
newDimensions3D[dimensions_, resolution_, matrixDimension_:3]:=Module[{matrices},matrices=expandedMatrix3D[#, resolution, matrixDimension]&/@dimensions;
Union@Flatten[matrices,1]];

discretizedSpectrum[desiredSpectrum_]:=Module[{discretizedGuess},
Table[Clip[desiredSpectrum, {0,1}], {\[Lambda], 400., 700., 1}];
];

color3D[dimension_, photonicSystem_]:=SpectrumToColor[photonicGlassSpectrum[dimension[[1]], dimension[[2]], photonicSystem]];

(*LAB Comparison*)

scoreLAB[specCoords_,guessSpectrum_]:=Module[{labCoords,discretizedGuess, zeroedGuess, xInt, yInt, zInt},
discretizedGuess=Table[guessSpectrum, {\[Lambda], 385., 745., 5}];
xInt=Total[discretizedGuess*xDiscrete]/17.172;
yInt=Total[discretizedGuess*yDiscrete]/17.172;
zInt=Total[discretizedGuess*zDiscrete]/17.172;
labCoords=List@@ColorConvert[{xInt,yInt,zInt},"XYZ"->"LAB"];
ColorDistance[LABColor@specCoords,LABColor@labCoords,DistanceFunction->"CIE2000"]
];

bestDimensionsLAB3D[colorCoords_, dimensionsVal_: dimensions3D, spectra3DVal_: spectra3D] := Quiet@Module[{specCoords,scores, rankings},
specCoords = List@@ColorConvert[colorCoords, "XYZ" -> "LAB"];
scores = Table[scoreLAB[specCoords, i], {i, spectra3DVal}];
rankings = Flatten[Position[scores, #] & /@ (Union@scores)[[;; 10]]];
dimensionsVal[[#]] & /@ rankings]; 

(*Finally, a 3D iterative search*)
    
dimensions3D=Tuples@{Range[25, 200, 10],Range[0.15,0.70,0.1]};
precision3D={5,0.05};

Iteration3D[color_, photonicSystem_, method_:"LABCompare"]:=Module[{ratioValues,bestDims, bestColors,colorDist,spectra3D,radiiCoords},
PrintTemporary[ProgressIndicator[Appearance -> "Indeterminate"]]; Pause[2];
spectra3D=photonicGlassSpectrum[#[[1]],#[[2]], photonicSystem]&/@dimensions3D;
bestDims=Which[method=="LABCompare",bestDimensionsLAB3D[color,dimensions3D,spectra3D],method=="XYZCompare", bestDimensionsXYZ3D[color,dimensions3D,spectra3D],method=="SubtractCompare", bestDimensionsSubtraction3D[color,dimensions3D,spectra3D],method=="RiemannCompare", bestDimensionsRiemann3D[color,dimensions3D,spectra3D], True, Print["Not a valid method, will default to LAB"];bestDimensionsLAB3D[color,dimensions3D,spectra3D]];
Print["These are the suggested dimensions:"];
Print[bestDims];
bestColors=color3D[#,photonicSystem]&/@bestDims;
Print["These are the corresponding colors:"];
Print[XYZColor/@bestColors];
colorDist=ColorDistance[#, XYZColor[color], DistanceFunction->"CIE2000"]&/@(XYZColor/@bestColors);
Print["These are the color distances in the CIELAB colorspace:"];
Print[colorDist];
Print["

The best match is spheres with a radius of "<>ToString[bestDims[[Position[colorDist, Min@colorDist][[1,1]]]][[1]]]<>" nm of Material 1 with a volume fraction of "<>ToString[bestDims[[Position[colorDist, Min@colorDist][[1,1]]]][[2]]]<>" within Material 2.
Here is the color of this match (right) compared to your desired color (left):"];
Print[Graphics[{XYZColor[color],Rectangle[],XYZColor@bestColors[[Position[colorDist, Min@colorDist][[1,1]]]],Rectangle[{1.2,0},{2.2,1}]},ImageSize->Small]];
Print[List@@XYZColor@bestColors[[Position[colorDist, Min@colorDist][[1,1]]]]];

precision3D=precision3D/2;
dimensions3D=newDimensions3D[bestDims,precision3D];

Print["This has a color distance of "<>ToString[Min@colorDist]<>" . 
If you want to keep iterating evaluate the Iteration3D function again"]
]

IterationReset3D:=Module[{},precision3D={5,0.05};dimensions3D=Tuples@{Range[25, 200, 10],Range[0.15,0.70,0.1]};]

End[ ]

EndPackage[ ]









