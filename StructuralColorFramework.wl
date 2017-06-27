(* ::Package:: *)

BeginPackage["StructuralColorFramework`"]

SpectrumToColor::usage =
        "SpectrumToColor[f[\[Lambda]]] uses the definition of the XYZ colorspace to convert a reflectance spectrum into an XYZ color"

ColorToSpectrum::usage= 
		"ColorToSpectrum[{X, Y, Z}] constructs a possible reflectance spectrum using the definition of the XYZ colorspace"

Begin["Private`"]

ChromaticityPlot["RGB"];

{x,y,z}=Interpolation[Thread[{Image`ColorOperationsDump`$wavelengths,#}]]&/@Transpose[Image`ColorOperationsDump`tris];

illuminant={{380,0.2661290322580645`},{385,0.32193548387096776`},{390,0.38225806451612904`},{395,0.4449193548387097`},{400,0.510483870967742`},{405,0.5791129032258064`},{410,0.6499999999999999`},{415,0.7220161290322581`},{420,0.7911290322580644`},{425,0.8532258064516128`},{430,0.9064516129032258`},{435,0.9495967741935484`},{440,0.9798387096774194`},{445,0.9955645161290323`},{450,1.`},{455,0.996774193548387`},{460,0.9927419354838709`},{465,0.9943548387096773`},{470,0.9983870967741935`},{475,1.000725806451613`},{480,0.9991935483870968`},{485,0.9912903225806452`},{490,0.9733870967741935`},{495,0.942741935483871`},{500,0.904032258064516`},{505,0.862741935483871`},{510,0.825`},{515,0.7968548387096774`},{520,0.7814516129032258`},{525,0.780483870967742`},{530,0.7903225806451613`},{535,0.8059677419354838`},{540,0.8233870967741935`},{545,0.8383064516129032`},{550,0.8483870967741935`},{555,0.8521774193548387`},{560,0.8491935483870967`},{565,0.8395967741935484`},{570,0.825`},{575,0.8076612903225807`},{580,0.7887096774193548`},{585,0.7695967741935484`},{590,0.7516129032258064`},{595,0.7356451612903225`},{600,0.7233870967741935`},{605,0.7163709677419354`},{610,0.7129032258064516`},{615,0.7112096774193548`},{620,0.7104838709677419`},{625,0.7101612903225807`},{630,0.7096774193548387`},{635,0.7085483870967741`},{640,0.7080645161290322`},{645,0.7095967741935483`},{650,0.7112903225806452`},{655,0.7112903225806452`},{660,0.7088709677419355`},{665,0.7033870967741935`},{670,0.6959677419354838`},{675,0.6879032258064516`},{680,0.6774193548387096`},{685,0.6629838709677419`},{690,0.646774193548387`},{695,0.6309677419354838`},{700,0.6153225806451612`},{705,0.5996774193548386`},{710,0.5838709677419355`},{715,0.567741935483871`},{720,0.5508064516129032`},{725,0.5346774193548387`},{730,0.5193548387096775`},{735,0.5064516129032258`},{740,0.4959677419354839`},{745,0.4854838709677419`},{750,0.4774193548387097`},{755,0.47177419354838707`},{760,0.4685483870967742`},{765,0.46774193548387094`},{770,0.4693548387096774`},{775,0.47177419354838707`},{780,0.47661290322580646`}};

illuminantInt=Interpolation[illuminant];

illuminantValues=illuminant[[All,2]];

xDiscrete=illuminantValues[[2;;74]]*Table[x[\[Lambda]],{\[Lambda],385,745,5}];
yDiscrete=illuminantValues[[2;;74]]*Table[y[\[Lambda]],{\[Lambda],385,745,5}];
zDiscrete=illuminantValues[[2;;74]]*Table[z[\[Lambda]],{\[Lambda],385,745,5}];

scalingConst=(NIntegrate[Through[{x,y,z}[\[Lambda]]]*illuminantInt[\[Lambda]],{\[Lambda],385,745}][[2]]);

SpectrumToColor[spectrum_]:=Re[NIntegrate[(Through[{x,y,z}[\[Lambda]]]*illuminantInt[\[Lambda]]*spectrum),{\[Lambda],385,745}, PrecisionGoal->2]]//#/scalingConst&

ColorToSpectrum[{xTristimulus_, yTristimulus_, zTristimulus_}]:=Module[{xWeight, yWeight, zWeight, function},
function=(xWeight*x[\[Lambda]]+yWeight*y[\[Lambda]]+zWeight*z[\[Lambda]])/.NSolve[{xWeight*NIntegrate[x[\[Lambda]]^2*illuminantInt[\[Lambda]],{\[Lambda],385,745}]+yWeight*NIntegrate[x[\[Lambda]]*y[\[Lambda]]*illuminantInt[\[Lambda]],{\[Lambda],385,745}]+
zWeight*NIntegrate[x[\[Lambda]]*z[\[Lambda]]*illuminantInt[\[Lambda]],{\[Lambda],385,745}]==xTristimulus*scalingConst,

yWeight*NIntegrate[y[\[Lambda]]^2*illuminantInt[\[Lambda]],{\[Lambda],385,745}]+xWeight*NIntegrate[y[\[Lambda]]*x[\[Lambda]]*illuminantInt[\[Lambda]],{\[Lambda],385,745}]+
zWeight*NIntegrate[y[\[Lambda]]*z[\[Lambda]]*illuminantInt[\[Lambda]],{\[Lambda],385,745}]==yTristimulus*scalingConst,

zWeight*NIntegrate[z[\[Lambda]]^2*illuminantInt[\[Lambda]],{\[Lambda],385,745}]+yWeight*NIntegrate[z[\[Lambda]]*y[\[Lambda]]*illuminantInt[\[Lambda]],{\[Lambda],385,745}]+
xWeight*NIntegrate[x[\[Lambda]]*z[\[Lambda]]*illuminantInt[\[Lambda]],{\[Lambda],385,745}]==zTristimulus*scalingConst}, {xWeight, yWeight, zWeight}(*\[Element]0\[LessEqual]#\[LessEqual]1&*)][[1]];
function]

expandedMatrix[dimension_, resolution_, matrixDimension_:3]:=Module[{matrixRange,potentialMatrix},
matrixRange=Range[-Floor[matrixDimension/2], Floor[matrixDimension/2]];
potentialMatrix=Table[{dimension[[1]]+resolution*a, dimension[[2]]+resolution*b}, {a, matrixRange}, {b, matrixRange}];
Cases[Flatten[potentialMatrix,1], {a_, b_}/;100<=a<=1000&&100<=b<=1000]];

newDimensions[dimensions_, resolution_, matrixDimension_:3]:=Module[{matrices},
matrices=expandedMatrix[#, resolution, matrixDimension]&/@dimensions;
Union@Flatten[matrices,1]];

discretizedSpectrum[desiredSpectrum_]:=Module[{discretizedGuess},
discretizedGuess=Table[Clip[desiredSpectrum, {0,1}], {\[Lambda], 400., 700., 1}];
#-Min@discretizedGuess&/@discretizedGuess];

scoreLAB[specCoords_,guessSpectrum_]:=Module[{labCoords,discretizedGuess, zeroedGuess, xInt, yInt, zInt},
discretizedGuess=Table[guessSpectrum, {\[Lambda], 385., 745., 5}];
xInt=Total[discretizedGuess*xDiscrete]/17.172;
yInt=Total[discretizedGuess*yDiscrete]/17.172;
zInt=Total[discretizedGuess*zDiscrete]/17.172;
labCoords=List@@ColorConvert[{xInt,yInt,zInt},"XYZ"->"LAB"];
ColorDistance[LABColor@specCoords,LABColor@labCoords,DistanceFunction->"CIE2000"]
];

End[ ]

EndPackage[ ]









