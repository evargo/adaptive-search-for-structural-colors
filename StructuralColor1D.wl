(* ::Package:: *)

BeginPackage["StructuralColor1D`",{"StructuralColorFramework`"}]

Spectrum1D::usage=
		"Spectrum1D[thickness1, thickness2, {n1, n2, layerNumber, theta}] returns a reflectance spectrum in terms of \[Lambda]. The default is a 5-layer crystal with material 1 as PMMA and material 2 as air, viewed from above"

Iteration1D::usage="Iteration1D[color, PhotonicSystem] returns color-matching results given the {X,Y,Z} coordinate of the desired color and a 1D PhotonicSystem."
				
PhotonicSystem1D::usage="PhotonicSystem1D[n1, n2, layerNumber, theta] returns a PhotonicsSystem to use in an Iteration1D."

IterationReset1D::usage="IterationReset[] is used to reset after an iteration is complete."
Begin["Private`"]

PhotonicSystem1D[nOne_,nTwo_,numberOfLayers_:5, incidentAngle_:0]:={nOne,nTwo,numberOfLayers, incidentAngle};

(*1D Spectrum calculation: originally published in http://demonstrations.wolfram.com/MultilayerPhotonicBandgap, 
full credit to the authors*)

kz[n_,\[Omega]_,c_,\[Beta]_]:=Sqrt[((n  \[Omega])/c)^2-\[Beta]^2]
k\[Theta][n_,\[Omega]_,c_,\[Theta]_]:=(n \[Omega])/c Cos[\[Theta]]
MTE[k1_,k2_,a_,b_] := Module[
{phase= Exp[I k1 a], invphase = Exp[-I k1 a],
sink2b = Sin[k2 b],cosk2b=Cos[k2 b], 
k1overk2 = k1/k2, k2overk1 = k2/k1,
sum,diff},
sum =I  (k1overk2 + k2overk1)/2;
diff = I  ( k1overk2 - k2overk1)/2;
{
{
phase (cosk2b + sum sink2b),-invphase diff sink2b
},
{
phase diff sink2b,invphase(cosk2b -  sum sink2b)
}
}
]
MTM[k1_,k2_,a_,b_, n1_, n2_] := Module[
{phase= Exp[I k1 a], invphase = Exp[-I k1 a],
sink2b = Sin[k2 b],cosk2b=Cos[k2 b], 
rat12sq , invrat12sq ,sum,diff},
rat12sq= (n2^2 k1)/(n1^2 k2);
invrat12sq= 1/rat12sq;

sum = I (rat12sq + invrat12sq)/2;
diff = I (rat12sq - invrat12sq)/2;
{
{
(cosk2b + sum sink2b)phase,(diff sink2b)invphase
},
{
(-diff sink2b)phase                 ,(cosk2b -  sum sink2b)invphase
}
}
]
K\[CapitalLambda][matrix_]:= ArcCos[Tr[matrix]/2]
 
CmagSquared[matrix_] := Abs[matrix[[2,1]]]^2

rE[ a_, b_,n1_,n2_, nn_, \[Omega]_, \[Theta]in_]:=
 Module[
{matrix,k\[Lambda], c= 3. 10^8,\[Theta]1,\[Theta]2},
\[Theta]1=ArcSin[Sin[\[Theta]in]/n1];
\[Theta]2= ArcSin[(n1 Sin[\[Theta]1])/n2];
matrix = MTE[k\[Theta][n1,\[Omega],c,\[Theta]1],k\[Theta][n2,\[Omega],c,\[Theta]2],a,b];

k\[Lambda] = K\[CapitalLambda][matrix];

1/(1 +Abs[ Sin[k\[Lambda]]/Sin[nn k\[Lambda]]]^2/CmagSquared[matrix])
]
rM[ a_, b_,n1_,n2_, nn_, \[Omega]_, \[Theta]in_]:=
 Module[
{matrix,k\[Lambda], c= 3. 10^8,\[Theta]1,\[Theta]2},
\[Theta]1=ArcSin[Sin[\[Theta]in]/n1];
\[Theta]2= ArcSin[(n1 Sin[\[Theta]1])/n2];
matrix = MTM[k\[Theta][n1,\[Omega],c,\[Theta]1],k\[Theta][n2,\[Omega],c,\[Theta]2],a,b,n1,n2];

k\[Lambda] = K\[CapitalLambda][matrix];

1/(1 +Abs[ Sin[k\[Lambda]]/Sin[nn k\[Lambda]]]^2/CmagSquared[matrix])
]

te=rE[lone/(10^9 2),ltwo/(10^9 2),n1,n2,numlay,((2 \[Pi]) 3 10^8)/(\[Lambda]/10^9),tin Pi/180];
tm =rM[lone/(10^9 2),ltwo/(10^9 2),n1,n2,numlay,((2 \[Pi]) 3 10^8)/(\[Lambda]/10^9),tin Pi/180];
nonPolarized = Mean[{te, tm}];

Spectrum1D[lOneThickness_, lTwoThickness_, {nOne_:1.4906, nTwo_:1, numberOfLayers_:5, incidentAngle_:0}]:=Module[{},nonPolarized/.{n1->nOne, n2->nTwo,numlay->numberOfLayers, tin->incidentAngle, lone->lOneThickness, ltwo->lTwoThickness}]

color1D[dimension_, photonicSystem_]:=SpectrumToColor[Spectrum1D[dimension[[1]], dimension[[2]], photonicSystem]]

(*LAB Comparison*)

bestDimensionsLAB1D[colorCoords_, dimensionsVal_: dimensions1D, spectra1DVal_: spectra1D] := Quiet@Module[{specCoords,scores, rankings},
specCoords = List@@ColorConvert[colorCoords, "XYZ" -> "LAB"];
scores = Table[scoreLAB[specCoords, i], {i, spectra1DVal}];
rankings = Flatten[Position[scores, #] & /@ (Union@scores)[[;; 10]]];
dimensionsVal[[#]] & /@ rankings]; 
    
dimensions1D=Tuples[Range[100, 1000, 50], 2];
precision1D=50;

Iteration1D[color_, photonicSystem_, method_:"LABCompare"]:=Module[{ratioValues,bestDims, bestColors,colorDist,spectra1D,radiiCoords},
PrintTemporary[ProgressIndicator[Appearance -> "Indeterminate"]]; Pause[2];
spectra1D=Table[Spectrum1D[val[[1]], val[[2]],photonicSystem], {val,dimensions1D}];
bestDims=bestDimensionsLAB1D[color,dimensions1D,spectra1D];
Print["These are the suggested dimensions:"];
Print[bestDims];
bestColors=color1D[#,photonicSystem]&/@bestDims;
Print["These are the corresponding colors:"];
Print[XYZColor/@bestColors];
colorDist=ColorDistance[#, XYZColor[color], DistanceFunction->"CIE2000"]&/@(XYZColor/@bestColors);
Print["These are the color distances in the CIELAB colorspace:"];
Print[colorDist];
Print["

The best match is "<>ToString[bestDims[[Position[colorDist, Min@colorDist][[1,1]]]][[1]]]<>" nm layers of Material 1 and "<>ToString[bestDims[[Position[colorDist, Min@colorDist][[1,1]]]][[2]]]<>" nm layers of Material 2.
Here is the color of this match (right) compared to your desired color (left):"];
Print[Graphics[{XYZColor[color],Rectangle[],XYZColor@bestColors[[Position[colorDist, Min@colorDist][[1,1]]]],Rectangle[{1.2,0},{2.2,1}]},ImageSize->Small]];

precision1D=precision1D/2;
dimensions1D=newDimensions[bestDims,precision1D];

Print["This has a color distance of "<>ToString[Min@colorDist]<>" . 
If you want to keep iterating evaluate the Iteration1D function again"]
];

IterationReset1D:=Module[{},precision1D=50;dimensions1D=Tuples[Range[100, 1000, 50], 2];]

End[ ]

EndPackage[ ]







