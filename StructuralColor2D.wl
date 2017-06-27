(* ::Package:: *)

BeginPackage["StructuralColor2D`",{"StructuralColorFramework`"}]

Iteration2D::usage="Iteration2D[color, PhotonicSystem] returns color-matching results given the {X,Y,Z} coordinate of the desired color and a 2D PhotonicSystem."

PhotonicSystem2D::usage="PhotonicSystem2D[n1, n2] returns a PhotonicsSystem to use in an Iteration2D. n1 is for the column material (usually air) and n2 is for the matrix material."

Begin["Private`"]

PhotonicSystem2D[nOne_,nTwo_]:={nOne^2, nTwo^2};

ratioList=Range[0.275, 0.5, 0.025];

writeInputFile[counter_,radius_:0.46,{eps1real_:9.0,eps2real_:1.0}]:=Module[{text},text=ToString[PaddedForm[eps1real,{5,4},NumberPadding->{""," "}]]<>"       epsrefreal\n"<>ToString[PaddedForm[eps1real,{5,4},NumberPadding->{""," "}]]<>"       eps1real\n"<>ToString[PaddedForm[eps2real,{3,4},NumberPadding->{""," "}]]<>"       eps2real\n"<>ToString[PaddedForm[radius,{3,4},NumberPadding->{""," "}]]<>"       radius";Export["/Users/emmavargo/inputfile"<>ToString[counter]<>".dat", text];]

convertTo\[Lambda][freq_, a_]:=a/(137*freq);

selectWavelengthData[a_, index_,orgData_]:=Module[{spectrum},
spectrum=Reverse@Cases[({convertTo\[Lambda][#[[1]], a], #[[2]]}&/@orgData[[index]]), {x_,y_}/;350<=x<=800];
#-{0,Min[spectrum[[All,2]]]}&/@spectrum];

addMoreRadii[listIn_]:=Module[{list,n,allOptions},
list=N@listIn;
n=1;
While[Length@Union@Flatten@list<10,list=Join[list,DeleteCases[0.0125/n+list,a_/;a>0.5],DeleteCases[-0.0125/n+list,a_/;a<=0.0]];n++];
Clip[#,{0,0.5}]&/@PadRight[DeleteDuplicates[Flatten[list],Abs[#1-#2]<=0.00001&],10,0.5][[1;;10]]]

color2D[{a_,r_},ratioListVal_:ratios,organizedData_]:=SpectrumToColor[Interpolation[selectWavelengthData[a,Position[ratioListVal,n_ /; n==r/a][[1,1]],organizedData]][Private`\[Lambda]]];

coords2D[{a_,r_},ratioListVal_:ratios]:={a,Position[ratioListVal,n_ /; n==r/a][[1,1]]};

bestDimensionsLAB2D[colorCoords_,ratioListVal_:ratioList,allSpectra2DVal_]:=Quiet@Module[{dimensionOptions,scores,options,accepted,accepted10,specCoords,spectrum,discretizedSpec,xSpec,ySpec,zSpec,rankings,bothMethodsPositions,favoredDimensionsBoth},
specCoords=List@@ColorConvert[colorCoords,"XYZ"->"LAB"];
scores=Table[scoreLAB[specCoords,i], {i,allSpectra2DVal}];
rankings=Flatten[Position[scores, #]&/@(Union@scores)[[;;30]]];
dimensionOptions=Flatten[Table[{a, ratioListVal[[i]]*a}, {a, 2000, 10000, 50}, {i, 1, 10}],1][[#]]&/@rankings;
(*"accepted" allows for a fabrication constraint or a set of fabrication constraints. {a_,b_}/;a>0 works as the null case*)
accepted=Cases[dimensionOptions,{a_,b_}/;b>250&&(a-2*b)>250];
If[Length@accepted<1, Print["you're being too picky"]];
PadRight[accepted, 10,accepted]]; 

Iteration2D[color_, photonicSystem_]:=Module[{ratioValues,bestDims, bestColors,colorDist,rawData,organizedData,locations2D,allSpectra2D,radiiCoords},
PrintTemporary[ProgressIndicator[Appearance -> "Indeterminate"]]; Pause[2];
ratioValues:=ratioList;
MapThread[writeInputFile[#1,#2,photonicSystem]&, {Range[20],N@PadRight[ratioValues,20,ratioValues]}];
Speak["If you wish to keep searching, please run the terminal file."];
Pause[30];
While[!FileExistsQ["/Users/emmavargo/markComplete"],Pause[30]];rawData=Partition[DeleteCases[Import["/Users/emmavargo/ref.dat", "Data"], a_/;Length[a]==3], 99];organizedData= MapThread[(rawData[[#1]]+rawData[[#2]])/2&, {Range[10], Range[11,20]}];

(*write functions*)
spectrum2D[a_,lineIndex_]:=Interpolation[selectWavelengthData[a,lineIndex,organizedData]];
locations2D=Flatten[Table[{a, ratioValues[[i]]*a}, {a, 2000, 10000,50}, {i, 1, 10}],1];
allSpectra2D=Flatten[Table[spectrum2D[a, i][\[Lambda]], {a, 2000, 10000, 50}, {i, 1, 10}],1];

(*find best colors*)
bestDims=bestDimensionsLAB2D[color,ratioList,allSpectra2D];
Print["These are the suggested dimensions:"];
Print[bestDims];
bestColors=color2D[#,ratioValues,organizedData]&/@bestDims;
Print["These are the corresponding colors:"];
Print[XYZColor/@bestColors];
(*Print[XYZColor[color]];*)
colorDist=ColorDistance[#, XYZColor[color], DistanceFunction->"CIE2000"]&/@(XYZColor/@bestColors);
Print["These are the color distances in the CIELAB colorspace:"];
Print[colorDist];
(*shows best guess*)
Print["

The best match has a lattice constant of "<>ToString[bestDims[[Position[colorDist, Min@colorDist][[1,1]]]][[1]]]<>" nm with a hole radius of "<>ToString[bestDims[[Position[colorDist, Min@colorDist][[1,1]]]][[2]]]<>" nm.
Here is the color of this match (right) compared to your desired color (left):"];
Print[Graphics[{XYZColor[color],Rectangle[],XYZColor@bestColors[[Position[colorDist, Min@colorDist][[1,1]]]],Rectangle[{1.2,0},{2.2,1}]},ImageSize->Small]];

(*set up for another iteration*)
radiiCoords=Union[(coords2D[#,ratioList]&/@bestDims[[Ordering[colorDist]]])[[All,2]]];
ratioList=addMoreRadii[ratioList[[radiiCoords]]];
Print["This has a color distance of "<>ToString[Min@colorDist]<>" . If you want to keep iterating evaluate the Iteration2D function again."];
]

End[ ]

EndPackage[ ]




