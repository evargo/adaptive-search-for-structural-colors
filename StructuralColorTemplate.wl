(* ::Package:: *)

BeginPackage["StructuralColorTemplate`",{"StructuralColorFramework`"}]


(* ::Text:: *)
(*Modify your usage statements depending on free and fixed parameters of your system :*)


SpectrumXX::usage=
		"SpectrumXX[freeVar1, freeVar2, {fixedParameter1,fixedParameter2}] returns a reflectance spectrum in terms of \[Lambda]."

IterationXX::usage="IterationXX[color, PhotonicSystem] returns color-matching results given the {X,Y,Z} coordinate of the desired color and an XX PhotonicSystem."
				
PhotonicSystemXX::usage="PhotonicSystemXX[fixedParameter1,fixedParameter2] returns a PhotonicsSystem to use in an IterationXX."

IterationResetXX::usage="IterationReset[] is used to reset after an iteration is complete."


Begin["Private`"]


(* ::Text:: *)
(*The PhotonicSystemXX function just served to put parameters into the order the spectrum function expects them.*)
(*  You can set default values, which will be used if no parameter is specified.*)


PhotonicSystemXX[freeParameter1_,freeParameter2_,freeParameter3_:defaultValue]:={freeParameter1,freeParameter2,freeParameter3};



(* ::Text:: *)
(*Here' s where you add your reflectance calculation. It can be broken into as many functions as needed, but at the end you should have a single function that takes in parameters and outputs a continuous spectrum.*)


(*code code code*)
spectrumXX[freeParameter1_, freeParameter2_, freeParameter3_:defaultValue,{fixedParameter1_,fixedParameter2_}]:=function[freeParameter1, freeParameter2,freeParameter3,fixedParameter1,fixedParameter2]

colorXX[dimension_, photonicSystem_]:=SpectrumToColor[SpectrumXX[dimension[[1]], dimension[[2]], photonicSystem]]

bestDimensionsLABXX[colorCoords_, dimensionsVal_: dimensionsXX, spectraXXVal_: spectraXX] := Quiet@Module[{specCoords,scores, rankings},
specCoords = List@@ColorConvert[colorCoords, "XYZ" -> "LAB"];
scores = Table[scoreLAB[specCoords, i], {i, spectraXXVal}];
rankings = Flatten[Position[scores, #] & /@ (Union@scores)[[;; 10]]];
dimensionsVal[[#]] & /@ rankings]; 
    


(* ::Text:: *)
(*Change the dimensions and precision to your initial set of dimensions (values of the free parameters) and precision. If you want different precisions for different parameters, look at the 3D package for an example.*)


dimensionsXX={{1,2,3},{1,2,3}};
precision1D=1;

IterationXX[color_, photonicSystem_]:=Module[{bestDims, bestColors,colorDist,spectraXX},
PrintTemporary[ProgressIndicator[Appearance -> "Indeterminate"]]; Pause[2];
spectraXX=Table[SpectrumXX[val[[1]], val[[2]],photonicSystem], {val,dimensionsXX}];
bestDims=bestDimensionsLABXX[color,dimensionsXX,spectraXX];
Print["These are the suggested dimensions:"];
Print[bestDims];
bestColors=colorXX[#,photonicSystem]&/@bestDims;
Print["These are the corresponding colors:"];
Print[XYZColor/@bestColors];
colorDist=ColorDistance[#, XYZColor[color], DistanceFunction->"CIE2000"]&/@(XYZColor/@bestColors);
Print["These are the color distances in the CIELAB colorspace:"];
Print[colorDist];
Print["

The best match is "<>ToString[bestDims[[Position[colorDist, Min@colorDist][[1,1]]]][[1]]]<>" (what the parameter corresponds to) and "<>ToString[bestDims[[Position[colorDist, Min@colorDist][[1,1]]]][[2]]]<>" (what the parameter correpsonds to).
Here is the color of this match (right) compared to your desired color (left):"];
Print[Graphics[{XYZColor[color],Rectangle[],XYZColor@bestColors[[Position[colorDist, Min@colorDist][[1,1]]]],Rectangle[{1.2,0},{2.2,1}]},ImageSize->Small]];

precisionXX=precisionXX/2;
dimensionsXX=newDimensions[bestDims,precisionXX];

Print["This has a color distance of "<>ToString[Min@colorDist]<>" . 
If you want to keep iterating evaluate the Iteration1D function again"]
];

IterationReset:=Module[{},precisionXX=1;dimensions1D={{1,2,3},{1,2,3}};]


End[ ]

EndPackage[ ]





(* ::Text:: *)
(*You did it! :)*)
