(* ::Package:: *)

If[$VersionNumber>=\[Infinity],
	BeginPackage["NVHamiltonian`", {"QuantumUtils`"}],
	BeginPackage["NVHamiltonian`", {"QuantumUtils`","VectorAnalysis`"}]
]


(* ::Section::Closed:: *)
(*Predicates*)


(* ::Subsection::Closed:: *)
(*Usage Declarations*)


CarbonShellPositionQ::usage = "CarbonShellPosition[pos] returns True iff pos is a {shell,index} pair.";
CarbonPhysicalPositionQ::usage = "CarbonPhysicalPosition[pos] returns True iff pos is a Cartesian vector.";
CarbonPositionQ::usage = "CarbonPositionQ[pos] returns True iff CarbonShellPositionQ or CarbonPhysicalPositionQ return True.";
CarbonPositionListQ::usage = "CarbonPositionList[carbonPositions] returns True iff carbonPositions is a list of carbon positions, i.e., each element satisfies CarbonPositionQ.";
DipolarHyperfineQ::usage = "DipolarHyperfineQ[x] returns True iff x is of the form {n,\"dipolar\"}.";
HyperfineTensorQ::usage = "HyperfineTensorQ[A] returns True iff A is a 3x3 matrix or a length 3 list";
NitrogenHyperfineSourceQ::usage = "NitrogenHyperfineSourceQ[str] returns True iff str is a valid nitrogen hyperfine source. See $nitrogenHyperfineTensorSources.";
NitrogenQuadrapolarSourceQ::usage = "NitrogenQuadrapolarSourceQ[str] returns True iff str is a valid nitrogen quadrapolar tensor source. See $nitrogenQuadrapolarTensorSources.";
CarbonHyperfineSourceQ::usage = "CarbonHyperfineSourceQ[str] returns True iff str is a valid carbon hyperfine source. See $carbonHyperfineTensorSources.";


(* ::Subsection::Closed:: *)
(*Implementations*)


Begin["`Private`"];


CarbonShellPositionQ[pos_]:=VectorQ[pos,IntegerQ]&&(Length[pos]===2)


CarbonPhysicalPositionQ[pos_]:=VectorQ[pos]&&(Length[pos]===3)


CarbonPositionQ[pos_]:=CarbonShellPositionQ[pos]||CarbonPhysicalPositionQ[pos]


CarbonPositionListQ[carbonPositions_]:=VectorQ[carbonPositions,CarbonPositionQ]


DipolarHyperfineQ[{n_,str_}]:=IntegerQ[n]&&(str==="dipolar")
DipolarHyperfineQ[str_]:=str==="dipolar"


HyperfineTensorQ[A_]:=(MatrixQ[A]||VectorQ[A])&&(Length[A]===3)


NitrogenHyperfineSourceQ[str_]:=MemberQ[NVHamiltonian`$nitrogenHyperfineTensorSources,str]


NitrogenQuadrapolarSourceQ[str_]:=MemberQ[NVHamiltonian`$nitrogenQuadrapolarTensorSources,str]


CarbonHyperfineSourceQ[str_]:=MemberQ[NVHamiltonian`$carbonHyperfineTensorSources,str]
CarbonHyperfineSourceQ[{n_,str_}]:=IntegerQ[n]&&MemberQ[NVHamiltonian`$carbonHyperfineTensorSources,str]


FalseQ[x_]:=TrueQ[Not[x]]


End[];


(* ::Section::Closed:: *)
(*Geometric/Rotational Definitions*)


(* ::Subsection:: *)
(*Usage Declarations*)


$e0::usage = "The crystal origin.";
$ez::usage = "The z unit direction in a tetrahedron centred around the origin";
$e1::usage = "The first bottom vertex of tetrahedron centred around the origin";
$e2::usage = "The second bottom vertex of tetrahedron centred around the origin";
$e3::usage = "The third bottom vertex of tetrahedron centred around the origin";
RotMat::usage = "RotMat[n,\[Theta]] returns the 3x3 rotation matrix by \[Theta] about the vector n using the right-hand-rule";
RotMatX::usage = "RotMatX[\[Theta]] returns the 3x3 rotation matrix by \[Theta] about the x direction using the right hand rule";
RotMatY::usage = "RotMatY[\[Theta]] returns the 3x3 rotation matrix by \[Theta] about the y direction using the right hand rule";
RotMatZ::usage = "RotMatZ[\[Theta]] returns the 3x3 rotation matrix by \[Theta] about the z direction using the right hand rule";
RotMatZYZ::usage = "RotMatZYZ[\[Theta]z1,\[Theta]y,\[Theta]z2] performs three rotations sequentially. These angles are what is meant by \"Extrinsic ZYZ Euler angles\".";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*First thing's first, get rid of this annoyance:*)


Unprotect[ArcTan];
ArcTan[0,0]=0;
Protect[ArcTan];


(* ::Text:: *)
(*Define various rotation matrices in O(3). All rotations of solid bodies will be done with extrinsic Euler angles using the ZYZ convention; hence we define the ZYZ rotation matrix.*)


RotMat[n_,\[Theta]_]:=With[{m=n/Sqrt[n.n]},MatrixExp[\[Theta]({
 {0, -m[[3]], m[[2]]},
 {m[[3]], 0, -m[[1]]},
 {-m[[2]], m[[1]], 0}
})]]
RotMatX[\[Theta]_]=RotMat[{1,0,0},\[Theta]];
RotMatY[\[Theta]_]=RotMat[{0,1,0},\[Theta]];
RotMatZ[\[Theta]_]=RotMat[{0,0,1},\[Theta]];
RotMatZYZ[\[Theta]z1_,\[Theta]y_,\[Theta]z2_]=RotMatZ[\[Theta]z2].RotMatY[\[Theta]y].RotMatZ[\[Theta]z1];


(* ::Text:: *)
(*Define the unit cell directions. If O is the centroid of a regular tetrahedron, then ArcCos(-1/3) is the angle AOB, where A and B are any two of the four vertices.*)
(*$ez lies in the z axis, and the other three are below the xy-plane forming a triangle. They are normalized so that the distance between any two vertices is unity (i.e. dist($ez,$e1)=dist($ez,$e2)=dist($ez,$e3)=dist($e1,$e2)=dist($e2,$e3)=dist($e3,$e1)=1).*)


$e0={0,0,0};
$ez={0,0,1};
$e1=RotMatY[ArcCos[-1/3]].$ez;
$e2=RotMatZ[2Pi/3].$e1;
$e3=RotMatZ[4Pi/3].$e1;


End[];


(* ::Section::Closed:: *)
(*NV Geometry*)


(* ::Subsection:: *)
(*Usage Declarations*)


$nitrogenPosition::usage = "The dimensionless location of the nitrogen in the NV PAS. Multiply this vector by \[Lambda] to get physical units.";
$vacancyPosition::usage = "The dimensionless location of the vacancy in the NV PAS. Multiply this vector by \[Lambda] to get physical units.";
CarbonDirections::usage = "CarbonDirections[{shell,index}] returns a path of unit vectors specifying how to get to the carbon at the specified shell and index from the origin in the NV PAS. Multiply this vector by \[Lambda] to get physical units.";
CarbonPositions::usage = "CarbonDirections[{shell,index}] returns the unitless position of the carbon at the specified shell and index in the NV PAS. Multiply this vector by \[Lambda] to get physical units.";
ShellSize::usage = "ShellSize[shell] returns the number of indeces in a given shell.";
NVAngles::usage = "NVAngles[nvOrientation] returns the Euler angles specifying how to rotate from the Crystal Frame to the given position.";
CarbonCoordinateList::usage = "CarbonCoordinateList[radius,latticeParameter:0.357] returns a list of of the positions of all carbons in Cartesian coordinates within radius of the NV at the center. The lattice is oriented such that the nitrogen is at (0,0,\[Lambda]). Default units are nm, but you can change this with latticeParameter.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Define where the nitrogen and vacancy are with respect to the labeled carbon positions.*)


$nitrogenPosition=$ez;
$vacancyPosition=$e0;


(* ::Text:: *)
(*Define the carbon position vectors. We do this by defining how to get to each Carbon from the origin using a sequence of the unit vectors above. This makes adding new shells trivial. A "shell" is defined as a collection of carbons with equivalent hyperfine tensors.*)
(*This is the correspondence between "shell" notation and Felton09's notation:*)
(*	shell=1  \[DoubleLongLeftRightArrow]  Subscript[C, a]*)
(*	shell=2  \[DoubleLongLeftRightArrow]  Subscript[C, b]*)
(*	shell=3  \[DoubleLongLeftRightArrow]  Subscript[C, d]*)
(*	shell=4  \[DoubleLongLeftRightArrow]  Subscript[C, f]*)
(*	shell=5  \[DoubleLongLeftRightArrow]  Subscript[C, g]*)
(*	shell=6  \[DoubleLongLeftRightArrow]  Subscript[C, h]*)


CarbonDirections[pos_?CarbonShellPositionQ]:=With[{shell=pos[[1]],index=pos[[2]]},Piecewise[{
{{$vacancyPosition,$e1},shell==1&&index==1},
{{$vacancyPosition,$e2},shell==1&&index==2},
{{$vacancyPosition,$e3},shell==1&&index==3},

{{$nitrogenPosition,-$e1},shell==2&&index==1},
{{$nitrogenPosition,-$e2},shell==2&&index==2},
{{$nitrogenPosition,-$e3},shell==2&&index==3},

{{$vacancyPosition,$e1,-$e2},shell==3&&index==1},
{{$vacancyPosition,$e1,-$e3},shell==3&&index==2},
{{$vacancyPosition,$e2,-$e3},shell==3&&index==3},
{{$vacancyPosition,$e2,-$e1},shell==3&&index==4},
{{$vacancyPosition,$e3,-$e1},shell==3&&index==5},
{{$vacancyPosition,$e3,-$e2},shell==3&&index==6},

{{$nitrogenPosition,-$e1,$e2},shell==4&&index==1},
{{$nitrogenPosition,-$e1,$e3},shell==4&&index==2},
{{$nitrogenPosition,-$e2,$e3},shell==4&&index==3},
{{$nitrogenPosition,-$e2,$e1},shell==4&&index==4},
{{$nitrogenPosition,-$e3,$e1},shell==4&&index==5},
{{$nitrogenPosition,-$e3,$e2},shell==4&&index==6},

{{$vacancyPosition,$e1,-$e2,$e1},shell==5&&index==1},
{{$vacancyPosition,$e1,-$e3,$e1},shell==5&&index==2},
{{$vacancyPosition,$e2,-$e3,$e2},shell==5&&index==3},
{{$vacancyPosition,$e2,-$e1,$e2},shell==5&&index==4},
{{$vacancyPosition,$e3,-$e1,$e3},shell==5&&index==5},
{{$vacancyPosition,$e3,-$e2,$e3},shell==5&&index==6},

{{$vacancyPosition,$e1,-$e3,$e2},shell==6&&index==1},
{{$vacancyPosition,$e2,-$e1,$e3},shell==6&&index==2},
{{$vacancyPosition,$e3,-$e2,$e1},shell==6&&index==3}
},{}]];
CarbonDirections[pos_?CarbonPhysicalPositionQ]:={$vacancyPosition,pos}


CarbonPositions[pos_?CarbonShellPositionQ]:=Total@CarbonDirections[pos]
CarbonPositions[pos_?CarbonPhysicalPositionQ]:=pos


(* ::Text:: *)
(*Return the number of indices in a given shell.*)


ShellSize[shell_]:=Piecewise[{
{3,shell==1},
{3,shell==2},
{6,shell==3},
{6,shell==4},
{6,shell==5},
{3,shell==6}
},{}];


(* ::Text:: *)
(*Define the 8 sets of Euler angles which define the 8 possible NV orientations within the lattice. The convention is to have the carbon-1 (i.e. shell-1, index-1 carbon) in the $ez direction whenever the PAS is not in the $ez direction, and in the $e1 direction whenever the PAS is in the $ez direction.*)


NVAngles[nvOrientation_]:=Piecewise[{
{{0,0,0},nvOrientation==1},
{{\[Pi],ArcCos[-1/3],0},nvOrientation==2},
{{\[Pi],ArcCos[-1/3],2 \[Pi]/3},nvOrientation==3},
{{\[Pi],ArcCos[-1/3],4\[Pi]/3},nvOrientation==4},
{{0,\[Pi],0},nvOrientation==5},
{{\[Pi],ArcCos[-1/3]-\[Pi],0},nvOrientation==6},
{{\[Pi],ArcCos[-1/3]-\[Pi],2\[Pi]/3},nvOrientation==7},
{{\[Pi],ArcCos[-1/3]-\[Pi],4\[Pi]/3},nvOrientation==8}
},{}];


CarbonCoordinateList[Radius_,LatticeParameter_:.357]:=
	Module[
	{
		cellcoords=0.25*{
			{0,0,0},{4,0,0},{0,4,0},{4,4,0},{2,2,0},
			{1,1,1},{3,3,1},
			{2,0,2},{0,2,2},{2,4,2},{4,2,2},
			{3,1,3},{1,3,3},
			{0,0,4},{4,0,4},{0,4,4},{4,4,4},{2,2,4}
		},
		N=Round[Radius/LatticeParameter]+1,
		output
	},
		output=Union[Flatten[Table[TranslationTransform[{j,k,l}][cellcoords],{j,-N,N,1},{k,-N,N,1},{l,-N,N,1}],{1,2,3,4}]];
		output=Select[output,0<LatticeParameter*Norm[#]<Radius&];
		output=(RotMatZYZ[-\[Pi]/4,-ArcCos[-1/3]/2,0].#&/@output)//Chop;
		output=Select[output,Not[#[[1]]==0&&#[[2]]==0&&0.433<#[[3]]<.434]&];
		output*LatticeParameter
	]


End[];


(* ::Section::Closed:: *)
(*Constants and Units*)


(* ::Subsection:: *)
(*Usage Declarations*)


$splitConst::usage = "$splitConst is the Zeeman splitting of an NV one would see at 1 Gauss when aligned.";
\[Gamma]e;\[Gamma]c;\[Gamma]n14;\[Gamma]n15;\[Mu]0;\[Lambda];\[HBar];
$s;$ms;$\[Mu]s;$ns;$Hz;$kHz;$MHz;$GHz;$Tesla;$Gauss;$Henry;$Joule;$m;$mm;$\[Mu]m;$nm;$rad;$rot;
$physicalConstants::usage = "A set of rules which define physical constants.";
$units::usage = "A set of rules defining common unit relationships.";
\[Gamma]n::usage = "\[Gamma]n[n] returns \[Gamma]n15 if n==15 and \[Gamma]n14 if n==14; all other inputs are left undefined.";


FillConstants::usage = "FillConstants[consts] takes the list of rules constants, and sets each Key to its Value. The default value for consts is $physicalConstants.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Various unit conversions in SI *)


$splitConst=5.604991461855897;


$units={
$s -> 1,
$ms -> 10^-3,
$\[Mu]s -> 10^-6,
$ns -> 10^-9,

$Hz -> 1,
$kHz -> 10^3,
$MHz -> 10^6,
$GHz -> 10^9,

$Tesla -> 1,
$Gauss -> 1/10000,
$Henry -> 1,
$Joule -> 1,

$m -> 1,
$mm -> 10^-3,
$\[Mu]m -> 10^-6,
$nm -> 10^-9,

$rad -> 1,
$rot -> 2\[Pi]
}


$physicalConstants={
\[Gamma]e -> 1.76086*10^5*$rad*$MHz/$Tesla,  (*gyromagnetic ratio of an electron*)
\[Gamma]c -> 67.262*$rad*$MHz/$Tesla,        (*gyromagnetic ratio of 13C*)
\[Gamma]n14 -> 19.331*$rad*$MHz/$Tesla,      (*gyromagnetic ratio of 14N*)
\[Gamma]n15 -> -27.116*$rad*$MHz/$Tesla,     (*gyromagnetic ratio of an 15N*)
\[Mu]0 -> 4\[Pi]*10^-7*$Henry/$m,             (*magnetic constant*)
\[Lambda] -> 0.154*$nm,                       (*distance between two sites in a diamond lattice*)
\[HBar] -> 1.054571726*10^\[Minus]34*$Joule*$s     (*reduced Planck constant*)
}/.$units


(* ::Text:: *)
(*A function to choose between nitrogen isotope symbols*)


\[Gamma]n[14]:=\[Gamma]n14;
\[Gamma]n[15]:=\[Gamma]n15;


FillConstants[consts_:$physicalConstants]:=If[NumberQ[#],#,#=N[#/.consts]]&/@consts[[All,1]];


End[];


(* ::Section::Closed:: *)
(*NV Parameters from Literature*)


(* ::Subsection:: *)
(*Usage Declarations*)


NitrogenDim::usage = "NitrogenDim[nitrogenIsotope] returns the dimension of the Hilbert space for a nitrogen atom with spin nitrogenIsotope.";
NitrogenSpin::usage = "NitrogenDim[nitrogenIsotope] returns the spin value for a nitrogen atom with isotope nitrogenIsotope.";
$typicalZFSConstant::usage = "A typical ZFS in Hz";
$nitrogenHyperfineTensorSources::usage = "A list of the valid nitrogen hyperfine source strings, e.g. \"Felton09\".";
$nitrogenQuadrapolarTensorSources::usage = "A list of the valid nitrogen quadrapolar tensor source strings, e.g. \"Felton09\".";
$carbonHyperfineTensorSources::usage = "A list of the valid carbon hyperfine source strings, e.g. \"Felton09\".";
$defaultNitrogenHyperfineTensorSource::usage = "The default nitrogen hyperfine source";
$defaultNitrogenQuadrapolarTensorSource::usage = "The default nitrogen quadrapolar source";
$defaultCarbonHyperfineTensorSource::usage = "The default carbon hyperfine source";
NitrogenHyperfineTensor::usage = "NitrogenHyperfineTensor[nitrogenIsotope,source] returns nitrogen hyperfine tensors from the literature. See code for source strings.";
NitrogenHyperfineTensorInfo::usage = "NitrogenHyperfineTensor[nitrogenIsotope,source] returns a description of where in the literature the given data was found. See code for source strings.";
NitrogenQuadrapolarTensor::usage = "NitrogenQuadrapolarTensor[nitrogenIsotope,source] returns nitrogen quadrapolar tensors from the literature. See code for source strings.";
NitrogenQuadrapolarTensorInfo::usage = "NitrogenQuadrapolarTensor[nitrogenIsotope,source] returns a description of where in the literature the given data was found. See code for source strings.";
CarbonHyperfineTensor::usage = "CarbonHyperfineTensor[nitrogenIsotope,source] returns carbon hyperfine tensors from the literature. See code for source strings.";
CarbonHyperfineTensorInfo::usage = "CarbonHyperfineTensor[nitrogenIsotope,source] returns a description of where in the literature the given data was found. See code for source strings.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Nitrogen dimension and spin*)


NitrogenDim[nitrogenIsotope_]:=Piecewise[{{3,nitrogenIsotope==14},{2,nitrogenIsotope==15}},1];
NitrogenSpin[nitrogenIsotope_]:=Piecewise[{{1,nitrogenIsotope==14},{1/2,nitrogenIsotope==15}},0];


(* ::Text:: *)
(*Near the standard value for the ZFS splitting frequency (Hz)*)


$typicalZFS=2.88*$GHz;


(* ::Text:: *)
(*Create the Nitrogen hyperfine tensors.*)


$defaultNitrogenHyperfineTensorSource = "Felton09";
$nitrogenHyperfineTensorSources = {"Felton09"};
NitrogenHyperfineTensor[nitrogenIsotope_,source_]:=Piecewise[{
{10^6*({
 {-2.70, 0, 0},
 {0, -2.70, 0},
 {0, 0, -2.14}
}),nitrogenIsotope==14&&source=="Felton09"},
{10^6*({
 {3.65, 0, 0},
 {0, 3.65, 0},
 {0, 0, 3.3}
}),nitrogenIsotope==15&&source=="Felton09"}
},({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
})];
NitrogenHyperfineTensor[nitrogenIsotope_]:=NitrogenHyperfineTensor[nitrogenIsotope,$defaultNitrogenHyperfineTensorSource];


(* ::Text:: *)
(*Define a function to keep track of where the numbers are coming from, etc.*)


NitrogenHyperfineTensorInfo[nitrogenIsotope_,source_]:=Piecewise[{
{"Felton et al., PRB 79, 075203 (2009). Table II. Room Temperature",nitrogenIsotope==14&&source=="Felton09"},
{"Felton et al., PRB 79, 075203 (2009). Table II. Room Temperature",nitrogenIsotope==15&&source=="Felton09"}
},0];
NitrogenHyperfineTensorInfo[nitrogenIsotope_]:=NitrogenHyperfineTensorInfo[nitrogenIsotope,$defaultNitrogenHyperfineTensorSource];


(* ::Text:: *)
(*Create the Nitrogen quadrapolar tensors*)


$defaultNitrogenQuadrapolarTensorSource = "Felton09";
$nitrogenQuadrapolarTensorSources = {"Felton09"};
NitrogenQuadrapolarTensor[nitrogenIsotope_,source_]:=Piecewise[{
{10^6*({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, -5.01}
}),nitrogenIsotope==14&&source=="Felton09"},
{({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
}),nitrogenIsotope==15&&source=="Felton09"}
},({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
})];
NitrogenQuadrapolarTensor[nitrogenIsotope_]:=NitrogenQuadrapolarTensor[nitrogenIsotope,$defaultNitrogenQuadrapolarTensorSource];


(* ::Text:: *)
(*Define a function to keep track of where the numbers are coming from, etc.*)


NitrogenQuadrapolarTensorInfo[nitrogenIsotope_,source_]:=Piecewise[{
{"Felton et al., PRB 79, 075203 (2009). Table II. Room Temperature",nitrogenIsotope==14&&source=="Felton09"},
{"Felton et al., PRB 79, 075203 (2009). Table II. Room Temperature",nitrogenIsotope==15&&source=="Felton09"}
},0];
NitrogenQuadrapolarTensorInfo[nitrogenIsotope_]:=NitrogenQuadrapolarTensorInfo[nitrogenIsotope,$defaultNitrogenQuadrapolarTensorSource];


(* ::Text:: *)
(*Define the carbon hyperfine tensors. Note that they depend on the nitrogen isotope, and hence it is important that the case where nitrogenIsotope=0 is taken into account, because that means we are ignoring nitrogen.*)


$defaultCarbonHyperfineTensorSource = "Felton09";
$carbonHyperfineTensorSources = {"Felton09"};
CarbonHyperfineTensor[shell_,nitrogenIsotope_,source_]:=Piecewise[{
{10^6*({
 {120.3, 0, 0},
 {0, 120.3, 0},
 {0, 0, 199.7}
}),shell==1&&nitrogenIsotope<=14&&source=="Felton09"},
{10^6*({
 {120.8, 0, 0},
 {0, 120.8, 0},
 {0, 0, 198.2}
}),shell==1&&nitrogenIsotope==15&&source=="Felton09"},
{10^6*({
 {13.26, 0, 0},
 {0, 13.26, 0},
 {0, 0, 18.49}
}),shell==5&&nitrogenIsotope<=15&&source=="Felton09"}
},({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
})];
CarbonHyperfineTensor[shell_,nitrogenIsotope_]:=CarbonHyperfineTensor[shell,nitrogenIsotope,$defaultCarbonHyperfineTensorSource];


(* ::Text:: *)
(*Define a function to keep track of where the numbers are coming from, etc.*)


CarbonHyperfineTensorInfo[shell_,nitrogenIsotope_,source_]:=Piecewise[{
{"Felton et al., PRB 79, 075203 (2009). Table III. At 10K.",shell==1&&nitrogenIsotope<=14&&source=="Felton09"},
{"Felton et al., PRB 79, 075203 (2009). Table III. Room temperature.",shell==1&&nitrogenIsotope==15&&source=="Felton09"},
{"Felton et al., PRB 79, 075203 (2009). Table III. At 10K.",shell==5&&nitrogenIsotope<=15&&source=="Felton09"}
},0];
CarbonHyperfineTensorInfo[shell_,nitrogenIsotope_]:=CarbonHyperfineTensorInfo[shell,nitrogenIsotope,$defaultCarbonHyperfineTensorSource];


End[];


(* ::Section::Closed:: *)
(*NV Frame change functions*)


(* ::Subsection:: *)
(*Usage Declarations*)


RotateBtoNVPAS::usage = "RotateBtoNVPAS[BLab,crystalOrientation,nvOrientation,BCoords] takes a vector, Blab, in the lab frame in BCoords coordinates and outputs the vector in the NV PAS with the same coordinate system.";
RotateBfromNVPAS::usage = "RotateBfromNVPAS[Bnvcoords,crystalOrientation,nvOrientation,BCoords] takes a vector, Bnvpas, in the NVPAS in BCoords coordinates and outputs the vector in the lab frame with the same coordinate system.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Define the function that rotates the B field in the lab axis into the NV PAS axis.*)
(*	BLab			The magnetic field in the lab frame in spherical coordinates ({B, \[Theta], \[Phi]}, *)
(*				where \[Theta] is the azimuthal angle, and \[Phi] is the polar angle)*)
(*	crystalOrientation 	The Euler angles describing the orientation of the crystal lattice with *)
(*				respect to the lab frame (ZYZ extrinsic convention).*)
(*	nvOrientation		An integer between 1 and 8 inclusive specifying the orientation of *)
(*				the NV centre with respect to the lattice. See NVAngles*)


RotateBtoNVPAS[BLab_,crystalOrientation_,nvOrientation_,BCoords_]:=
CoordinatesFromCartesian[
RotMatZYZ[Sequence@@Reverse@-NVAngles[nvOrientation]].RotMatZYZ[Sequence@@Reverse@-crystalOrientation].CoordinatesToCartesian[BLab,BCoords],
BCoords
]


RotateBfromNVPAS[Bnvpas_,crystalOrientation_,nvOrientation_,BCoords_]:=
CoordinatesFromCartesian[
RotMatZYZ[Sequence@@crystalOrientation].RotMatZYZ[Sequence@@NVAngles[nvOrientation]].CoordinatesToCartesian[Bnvpas,BCoords],
BCoords
]


(* ::Text:: *)
(*Suppose A is the hyperfine tensor between a spin at location $ez+$vacancyPosition and the NV. This function rotates the tensor so that the spin now sits at position vec+$vacancyPosition.*)
(*	A	The hyperfine tensor, a 3x3 matrix*)
(*	vec	Where to rotate*)


RotateTensorToNVPAS[A_, vec_]:=With[
{angles=CoordinatesFromCartesian[vec-$vacancyPosition,Spherical][[{2,3}]]},
With[{R=RotMatZ[angles[[2]]].RotMatY[angles[[1]]]},
R.A.R\[Transpose]
]]


End[];


(* ::Section::Closed:: *)
(*Pieces of the Total NV Hamiltonian*)


(* ::Subsection:: *)
(*Usage Declarations*)


ZFSHamiltonian::usage = "ZFSHamiltonian[\[CapitalDelta]] returns the zero-field-splitting Hamiltonian of the NV Centre. \[CapitalDelta] can be a number, the diagonal elements of the tensor, or the whole tensor.";
NVZeemanHamiltonian::usage = "NVZeemanHamiltonian[B] returns the Zeeman Hamiltonian for the NV Centre. B is in cartesian coordinates, and in units of Gauss.";
NitrogenZeemanHamiltonian::usage = "NitrogenZeemanHamiltonian[B,nitrogenIsotope] returns the Zeeman Hamiltonian for a nitrogen with specified isotope. B is in cartesian coordinates, and in units of Gauss.";
CarbonZeemanHamiltonian::usage = "CarbonZeemanHamiltonian[B,numCarbon] returns the Zeeman Hamiltonian for a specified number of carbon 13s. B is in cartesian coordinates, and in units of Gauss.";
NitrogenHyperfineHamiltonian::usage = "NitrogenHyperfineHamiltonian[P, nitrogenIsotope] returns the hyperfine Hamiltonian for a Nitrogen coupled to an electron for a given tensor A. A can be the 3x3 tensor, or just the diagonal elements, and is in units of Hz.";
NitrogenQuadrapolarHamiltonian::usage = "NitrogenQuadrapolarHamiltonian[A, nitrogenIsotope] returns the quadrapolar Hamiltonian for a Nitrogen coupled to an electron for a given tensor P. P can be the 3x3 tensor, or just the diagonal elements, and is in units of Hz.";
CarbonHyperfineHamiltonian::usage = "CarbonHyperfineHamiltonian[As, nitrogenIsotope] returns the hyperfine Hamiltonian for a bunch  of carbon 13s coupled to an electron for a given list of tensors A. The tensors in A can be the 3x3 tensor, or just the diagonal elements, and is in units of Hz.";
NitrogenCarbonDipolarHamiltonian::usage = "NitrogenCarbonDipolarHamiltonian[nitrogenIsotope,carbonPositions] returns the dipolar Hamiltonian of a bunch of carbons in the position list with a nitrogen atom. See CarbonPositions.";
CarbonCarbonDipolarHamiltonian::usage = "CarbonCarbonDipolarHamiltonian[carbonPositions] returns the dipolar Hamiltonian between all of the Carbons in the position list. See CarbonPositions.";
DipolarHyperfineTensor::usage = "DipolarHyperfineTensor[pos] computes the hyperfine tensor due to a dipole-dipole coupling between a spin at cartesian coordinate pos (in PAS frame) and an electron at $vacancyPosition.";
DipolarVector::usage = "DipolarVector[carbonposition,inputunit:$nm,outputunit:$kHz] returns the Cartesian dipole coupling vector between a spin at cartesian coordinate carbonposition (in PAS frame) and an electron at the origin. Default input units are nm, and output units are kHz.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*The Zero-Field-Splitting Hamiltonian. Notice there are three patterns for this function: one where you enter just a number for the Sz^2 coefficient, one where you enter a lits of the diagonal coefficients, and one where you enter the entire 3x3 tensor.*)


ZFSHamiltonian[\[CapitalDelta]_]:=2\[Pi]*\[CapitalDelta]*SpinZ[1].SpinZ[1];
ZFSHamiltonian[\[CapitalDelta]_?ListQ]:=2\[Pi]*(\[CapitalDelta][[1]]SpinX[1].SpinX[1]+\[CapitalDelta][[2]]SpinY[1].SpinY[1]+\[CapitalDelta][[3]]SpinZ[1].SpinZ[1]);
ZFSHamiltonian[\[CapitalDelta]_?MatrixQ]:=2\[Pi]*(Flatten[Outer[Dot[#1,#2]&,{x,y,z},{x,y,z}]].Flatten[\[CapitalDelta]]/.{x->SpinX[1],y->SpinY[1],z->SpinZ[1]});


(* ::Text:: *)
(*The Zeeman Hamiltonians:*)


NVZeemanHamiltonian[B_]:=ZeemanHamiltonian[\[Gamma]e,B $Gauss,1];


NitrogenZeemanHamiltonian[B_,nitrogenIsotope_]:=ZeemanHamiltonian[\[Gamma]n[nitrogenIsotope],B $Gauss,NitrogenSpin[nitrogenIsotope]]


CarbonZeemanHamiltonian[B_,numCarbon_]:=
With[{zeemanHam=ZeemanHamiltonian[\[Gamma]c,B $Gauss,1/2]},
Total@Table[Subscript[\[DoubleStruckOne],2^(n-1)]\[CircleTimes]zeemanHam\[CircleTimes]Subscript[\[DoubleStruckOne],2^(numCarbon-n)],{n,numCarbon}]]


(* ::Text:: *)
(*The Nitrogen quadrapolar Hamiltonian. The tensor is expected in units of Hz (and NOT rad/s).*)


NitrogenQuadrapolarHamiltonian[P_?MatrixQ,nitrogenIsotope_]:=
With[{S=Spin[NitrogenSpin[nitrogenIsotope]]},
2\[Pi]*(Flatten[Outer[Dot[#1,#2]&,{x,y,z},{x,y,z}]].Flatten[P]/.{x->S[[1]],y->S[[2]],z->S[[3]]})
]
NitrogenQuadrapolarHamiltonian[P_?ListQ,nitrogenIsotope_]:=
With[{S=Spin[NitrogenSpin[nitrogenIsotope]]},
2\[Pi]*(P[[1]]S[[1]].S[[1]]+P[[2]]S[[2]].S[[2]]+P[[3]]S[[3]].S[[3]])
]
NitrogenQuadrapolarHamiltonian[P_,0]:=0


(* ::Text:: *)
(*The hyperfine Hamiltonians. The tensor is expected in units of Hz (and NOT rad/s).*)


NitrogenHyperfineHamiltonian[nitrogenHyperfineTensor_,nitrogenIsotope_]:=
HyperfineHamiltonian[2\[Pi]*nitrogenHyperfineTensor,1,NitrogenSpin[nitrogenIsotope]]
NitrogenHyperfineHamiltonian[2\[Pi]*nitrogenHyperfineTensor_,0]:=0


CarbonHyperfineHamiltonian[carbonHyperfineTensors_,nitrogenIsotope_]:=
With[
{numCarbon=Length[carbonHyperfineTensors],
carbonHyperfineTensorList=Sort[carbonHyperfineTensors, #1[[1]] < #2[[1]] &][[All, 2]]},
Total@Table[
HyperfineHamiltonian[2\[Pi]*carbonHyperfineTensorList[[n]],1,1/2,0,NitrogenDim[nitrogenIsotope]*2^(n-1),2^(numCarbon-n)],{n,numCarbon}]
]


(* ::Text:: *)
(*The dipolar Hamiltonians:*)


NitrogenCarbonDipolarHamiltonian[nitrogenIsotope_,carbonPositions_]:=
With[{numCarbon=Length[carbonPositions]},
Total@Table[
DipolarHamiltonian[
\[Lambda]*($nitrogenPosition-CarbonPositions[carbonPositions[[n]]]),
\[Mu]0*\[Gamma]n[nitrogenIsotope]*\[Gamma]c*\[HBar],
NitrogenSpin[nitrogenIsotope],
1/2,
0,
2^(n-1),
2^(numCarbon-n)],
{n,numCarbon}]]


CarbonCarbonDipolarHamiltonian[carbonPositions_]:=
With[{numCarbon=Length[carbonPositions]},
Total@Flatten[Table[
DipolarHamiltonian[
\[Lambda]*(CarbonPositions[carbonPositions[[n]]]-CarbonPositions[carbonPositions[[m]]]),
\[Mu]0*\[Gamma]c*\[Gamma]c*\[HBar],
1/2,
1/2,
2^(n-1),
2^(m-n-1),
2^(numCarbon-m)],
{n,1,numCarbon},{m,n+1,numCarbon}],1]]


(* ::Text:: *)
(*Compute the hyperfine tensor due to a dipole-dipole coupling between a carbon at pos and an electron at $vacancyPosition. We divide by 2\[Pi] because of the 2\[Pi] in the CarbonHyperfineHamilonian function.*)


DipolarHyperfineTensor[pos_]:=Module[{R,r},
r=\[Lambda]*(pos-$vacancyPosition);R=Sqrt[r.r];
(((-\[Mu]0 \[Gamma]c \[Gamma]e \[HBar])/(4\[Pi] R^3))*(Subscript[\[DoubleStruckOne], 3]-3*{r}\[Transpose].{r}/R^2))/(2\[Pi])
]


DipolarVector[carbonposition_,inputunit_:$nm,outputunit_:$kHz]:=
	Module[
		{R,r,e},
		r=carbonposition*inputunit;
		R=Norm[r];
		e=r/R;
		(\[Mu]0 \[Gamma]e \[Gamma]c \[HBar])/(8 \[Pi]^2 R^3) {3e[[3]]e[[1]],3e[[3]]e[[2]],3e[[3]]^2-1}/outputunit
	]


End[];


(* ::Section::Closed:: *)
(*NV Parameter Object Functions*)


(* ::Subsection:: *)
(*Usage Declarations*)


FillInNVParameters::usage = "FillInNVParameters[params] takes a (possibly) incomplete list of NV Parameters and inserts default values for all those parameters missing. Keywords for tensors (such as \"Felton09\" or \"dipolar\") are also replaced by the correct tensor matrices.";
GrabCarbonHyperfineTensor::usage = "GrabCarbonHyperfineTensorQ[carbonPosition,input,nitrogenIsotope] returns a 3x3 hyperfine tensor, where input can be (1) a 3x3 hyperfine tensor, (2) the empty list {}, (3) the string \"dipolar\", (4) a valid carbon hyperfine source.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*GrabNitrogenqQuadrapolarTensor[A.isotope] returns a 3x3 tensor where A is a valid source string, or the tensor itself*)


GrabNitrogenQuadrapolarTensor[source_?NitrogenQuadrapolarSourceQ,nitrogenIsotope_]:=NitrogenQuadrapolarTensor[nitrogenIsotope,source];
GrabNitrogenQuadrapolarTensor[A_?HyperfineTensorQ,nitrogenIsotope_]:=A;


(* ::Text:: *)
(*GrabNitrogenHyperfineTensor[A.isotope] returns a 3x3 tensor where A is a valid source string, or the tensor itself*)


GrabNitrogenHyperfineTensor[source_?NitrogenHyperfineSourceQ,nitrogenIsotope_]:=NitrogenHyperfineTensor[nitrogenIsotope,source];
GrabNitrogenHyperfineTensor[A_?HyperfineTensorQ,nitrogenIsotope_]:=A;


(* ::Text:: *)
(*The following set the default behaviour in all four situations of carbonposition type and analytic form*)


GrabCarbonHyperfineTensor[n_,carbonPosition_?CarbonShellPositionQ,{},nitrogenIsotope_,analyticForm_?FalseQ]:=RotateTensorToNVPAS[CarbonHyperfineTensor[carbonPosition[[1]],nitrogenIsotope,$defaultCarbonHyperfineTensorSource],CarbonPositions[carbonPosition]];
GrabCarbonHyperfineTensor[n_,carbonPosition_?CarbonPositionQ,{},nitrogenIsotope_,analyticForm_?TrueQ]:=Symbol["Global`"<>#<>ToString[n]]&/@{"A","B","C"};
GrabCarbonHyperfineTensor[n_,carbonPosition_?CarbonPhysicalPositionQ,{},nitrogenIsotope_,analyticForm_?FalseQ]:={0,0,0};
GrabCarbonHyperfineTensor[n_,carbonPosition_?CarbonPhysicalPositionQ,{},nitrogenIsotope_,analyticForm_?TrueQ]:=Symbol["Global`"<>#<>ToString[n]]&/@{"A","B","C"};


(* ::Text:: *)
(*The following handle when the user has actually input something*)


GrabCarbonHyperfineTensor[n_,carbonPosition_?CarbonPositionQ,source_?CarbonHyperfineSourceQ,nitrogenIsotope_,analyticForm_]:=RotateTensorToNVPAS[CarbonHyperfineTensor[carbonPosition[[1]],nitrogenIsotope,source[[2]]],CarbonPositions[carbonPosition]];
GrabCarbonHyperfineTensor[n_,carbonPosition_?CarbonPositionQ,dipolarString_?DipolarHyperfineQ,nitrogenIsotope_,analyticForm_]:=DipolarHyperfineTensor[CarbonPositions[carbonPosition]];
GrabCarbonHyperfineTensor[n_,carbonPosition_?CarbonPositionQ,hyperfineTensor_?HyperfineTensorQ,nitrogenIsotope_,analyticForm_]:=hyperfineTensor;


(* ::Text:: *)
(*Takes whatever the user has put in (or not put in) to the "carbonHyperfineTensors" option and puts in/calculates/gets the default 3x3 tensors*)


SelectCarbonTensor[hyperfineTensors_,n_]:=With[{s=Flatten[Select[hyperfineTensors,(#[[1]]===n)&,1],1]},If[s==={},{},s[[2]]]];
FillInCarbonHyperfineTensors[carbonPositions_,hyperfineTensors_,nitrogenIsotope_,analyticForm_]:=
With[{numCarbon=Length[carbonPositions]},
Table[{n,GrabCarbonHyperfineTensor[n,carbonPositions[[n]],SelectCarbonTensor[hyperfineTensors,n],nitrogenIsotope,analyticForm]},{n,numCarbon}]
]


(* ::Text:: *)
(*This function fills in any missing parameters with the default ones.*)


FillInNVParameters[params_]:=
Module[{totalParams,
nitrogenHyperfineTensor,
nitrogenQuadrapolarTensor,
carbonHyperfineTensors,
analyticForm="analyticForm"/.Join[params,{"analyticForm"->False}],
nitrogenIsotope="nitrogenIsotope"/.Join[params,{"nitrogenIsotope"->0}]},
totalParams=Join[params,{
"\[CapitalDelta]"->Global`\[CapitalDelta],
"B"->{0,0,0},
"crystalOrientation"->{0,0,0},
"nvOrientation"->1,
"nitrogenIsotope"->nitrogenIsotope,
"carbonPositions"->{},
"nitrogenHyperfineTensor"->If[analyticForm,{Global`AN,Global`BN,Global`CN},$defaultNitrogenHyperfineTensorSource],
"nitrogenQuadrapolarTensor"->If[analyticForm,If[nitrogenIsotope==14,{Global`DN,Global`EN,Global`FN},{0,0,0}],$defaultNitrogenQuadrapolarTensorSource],
"carbonHyperfineTensors"->{},
"zfsActivated"->True,
"nvZeemanActivated"->True,
"nitrogenZeemanActivated"->True,
"carbonZeemanActivated"->True,
"nitrogenQuadrapolarActivated"->True,
"nitrogenHyperfineActivated"->True,
"carbonHyperfineActivated"->True,
"nitrogenCarbonDipolarActivated"->True,
"carbonCarbonDipolarActivated"->True,
"analyticForm"->analyticForm,
"BCoords"->Spherical}];
nitrogenHyperfineTensor=GrabNitrogenHyperfineTensor["nitrogenHyperfineTensor"/.totalParams,nitrogenIsotope];
nitrogenQuadrapolarTensor=GrabNitrogenQuadrapolarTensor["nitrogenQuadrapolarTensor"/.totalParams,nitrogenIsotope];
carbonHyperfineTensors=FillInCarbonHyperfineTensors[Sequence@@({"carbonPositions","carbonHyperfineTensors",nitrogenIsotope,analyticForm}/.totalParams)];
totalParams=DeleteCases[totalParams,_?(#[[1]]==="nitrogenHyperfineTensor"&)];
totalParams=DeleteCases[totalParams,_?(#[[1]]==="nitrogenQuadrapolarTensor"&)];
totalParams=DeleteCases[totalParams,_?(#[[1]]==="carbonHyperfineTensors"&)];
totalParams=Join[totalParams,{
"nitrogenHyperfineTensor"->nitrogenHyperfineTensor,
"nitrogenQuadrapolarTensor"->nitrogenQuadrapolarTensor,
"carbonHyperfineTensors"->carbonHyperfineTensors}];
GatherBy[totalParams, #[[1]] &][[All, 1]] (* this gets rid of duplicate Rules - purely asthetic *)
]


End[];


(* ::Section::Closed:: *)
(*The NV Hamiltonian*)


(* ::Subsection:: *)
(*Usage Declarations*)


NVHamiltonian::usage = "NVHamiltonian[\[CapitalDelta],B,crystalOrientation,nvOrientation,nitrogenIsotope,nitrogenHyperfineTensor,nitrogenQuadrapolarTensor,carbonPositions,carbonHyperfineTensors] returns the Hamiltonian with the given parameters. See code for details.";
NVHamiltonian::usage = "NVHamiltonian[params] returns the Hamiltonian corresponding to whatever is set in the params object. See code for details.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*The NVHamiltonianExplicitForm with explicit inputs. It is expected that the average user will never call this pattern; all of the hyperfines need to be explicitly entered.*)
(*Use the NVHamiltonian definition that comes right after this one.*)
(*A breakdown of the inputs:*)
(*	\[CapitalDelta]					The ZFS tensor. Can be a number, a list of the diagonal elements, or the full 3x3 tensor.*)
(*	B					The magnetic field in the lab frame in spherical coordinates ({B, \[Theta], \[Phi]}, where \[Theta] is the azimuthal angle, and \[Phi] is the polar angle). Units of Gauss.*)
(*	crystalOrientation 			The Euler angles describing the orientation of the crystal lattice with respect to the lab frame (ZYZ extrinsic convention).*)
(*	nvOrientation				An integer between 1 and 8 inclusive specifying the orientation of the NV centre with respect to the lattice. See NVAngles*)
(*	nitrogenIsotope			The spin of the nitrogen (so either 14 or 15). Optionally, enter 0 to ignore the nitrogen.*)
(*	nitrogenHyperfineTensor		The 3x3 nitrogen hyperfine tensor, or a list of the diagonal elements.*)
(*	nitrogenQuadrapolarTensor		The 3x3 nitrogen quadrapolar tensor, or a list of the diagonal elements.*)
(*	carbonPositions 			A list of shell-index pairs, specifying which carbons have spin, e.g. {{1,1},{1,2},{1,3}} for all three closest carbons. This is the order they will take in the tensor product.*)
(*	carbonHyperfineTensors		A list of the 3x3 carbon hyperfine tensors (or diagonal elements), one for each of the entries in carbonPositions.*)
(*	zfsActivated				True or False; whether the ZFS Hamiltonian is to be included.*)
(*	nvZeemanActivated			True or False; whether the NV Zeeman Hamiltonian is to be included.*)
(*	nitrogenZeemanActivated		True or False; whether the nitrogen Zeeman Hamiltonian is to be included.*)
(*	carbonZeemanActivated		True or False; whether the carbon Zeeman Hamiltonian is to be included.*)
(*	nitogenQuadrapolarActivated	True or False; whether the nitrogen quadrapolar Hamiltonian is to be included.*)
(*	nitrogenHyperfineActivated		True or False; whether the nitrogen hyperfine Hamiltonian is to be included.*)
(*	carbonHyperfineActivated		True or False; whether the carbon hyperfine Hamiltonian is to be included.*)
(*	nitrogenCarbonDipolarActivated	True or False; whether the nitrogen-carbon dipolar Hamiltonian is to be included.*)
(*	carbonCarbonDipolarActivated	True or False; whether the carbon-carbon dipolar Hamiltonian is to be included.*)
(*	analyticForm				True or False; true to avoid numerical values*)
(*	BCoords				A coordinate system from the VectorAnalysis package, such as Spherical or Cartesian*)
(*The order of the tensor product structure is Subscript[\[ScriptCapitalH], NV]\[CircleTimes]Subscript[\[ScriptCapitalH], Nitrogen]\[CircleTimes]Subscript[\[ScriptCapitalH], Carbon 1]\[CircleTimes]Subscript[\[ScriptCapitalH], Carbon 2]\[CircleTimes]...\[CircleTimes]Subscript[\[ScriptCapitalH], Carbon m]where Length[carbonPositions]=m.*)


NVHamiltonianExplicitForm[\[CapitalDelta]_,B_,crystalOrientation_,nvOrientation_,
nitrogenIsotope_,nitrogenHyperfineTensor_,nitrogenQuadrapolarTensor_,
carbonPositions_,carbonHyperfineTensors_,
zfsActivated_,
nvZeemanActivated_,nitrogenZeemanActivated_,carbonZeemanActivated_,
nitrogenQuadrapolarActivated_,nitrogenHyperfineActivated_,carbonHyperfineActivated_,
nitrogenCarbonDipolarActivated_,carbonCarbonDipolarActivated_,analyticForm_,BCoords_]:=
With[{
Beff=CoordinatesToCartesian[RotateBtoNVPAS[B,crystalOrientation,nvOrientation,BCoords],BCoords],
numCarbon=Length[carbonPositions],
nvDim=3,
nitrogenDim=NitrogenDim[nitrogenIsotope],
totalCarbonDim=2^Length[carbonPositions],
numericReplacements=If[analyticForm,{},Join[$physicalConstants,$units,{\[CapitalDelta]->($typicalZFS)/.$units}]]
},
(
If[zfsActivated,ZFSHamiltonian[\[CapitalDelta]]\[CircleTimes]Subscript[\[DoubleStruckOne], nitrogenDim]\[CircleTimes]Subscript[\[DoubleStruckOne], totalCarbonDim],0]+
If[nvZeemanActivated,NVZeemanHamiltonian[Beff]\[CircleTimes]Subscript[\[DoubleStruckOne], nitrogenDim]\[CircleTimes]Subscript[\[DoubleStruckOne], totalCarbonDim],0]+
If[nitrogenZeemanActivated,Subscript[\[DoubleStruckOne], nvDim]\[CircleTimes]NitrogenZeemanHamiltonian[Beff,nitrogenIsotope]\[CircleTimes]Subscript[\[DoubleStruckOne], totalCarbonDim],0]+
If[carbonZeemanActivated,Subscript[\[DoubleStruckOne], nvDim]\[CircleTimes]Subscript[\[DoubleStruckOne], nitrogenDim]\[CircleTimes]CarbonZeemanHamiltonian[Beff,numCarbon],0]+

If[nitrogenHyperfineActivated,NitrogenHyperfineHamiltonian[nitrogenHyperfineTensor,nitrogenIsotope]\[CircleTimes]Subscript[\[DoubleStruckOne], totalCarbonDim],0]+
If[nitrogenQuadrapolarActivated,Subscript[\[DoubleStruckOne], nvDim]\[CircleTimes]NitrogenQuadrapolarHamiltonian[nitrogenQuadrapolarTensor,nitrogenIsotope]\[CircleTimes]Subscript[\[DoubleStruckOne], totalCarbonDim],0]+
If[nitrogenCarbonDipolarActivated,Subscript[\[DoubleStruckOne], nvDim]\[CircleTimes]NitrogenCarbonDipolarHamiltonian[nitrogenIsotope,carbonPositions],0]+

If[carbonHyperfineActivated,CarbonHyperfineHamiltonian[carbonHyperfineTensors,nitrogenIsotope],0]+
If[carbonCarbonDipolarActivated,Subscript[\[DoubleStruckOne], nvDim]\[CircleTimes]Subscript[\[DoubleStruckOne], nitrogenDim]\[CircleTimes]CarbonCarbonDipolarHamiltonian[carbonPositions],0]
)/.numericReplacements
]


(* ::Text:: *)
(*The NVHamiltonian function that the average user should call. There is only one argument, params, which should be a list of Rules, e.g., params={"B"->{1000,0,0},"nitrogenIsotope"->15}.*)
(*The rules that can be assigned are as follows:*)
(*	"\[CapitalDelta]"					The ZFS tensor. Can be a number, a list of the diagonal elements, or the full 3x3 tensor.*)
(*	"B"					The magnetic field in the lab frame in spherical coordinates ({B, \[Theta], \[Phi]}, where \[Theta] is the azimuthal angle, and \[Phi] is the polar angle). Units of Gauss.*)
(*	"crystalOrientation "			The Euler angles describing the orientation of the crystal lattice with respect to the lab frame (ZYZ extrinsic convention).*)
(*	"nvOrientation"			An integer between 1 and 8 inclusive specifying the orientation of the NV centre with respect to the lattice. See NVAngles*)
(*	"nitrogenIsotope"			The spin of the nitrogen (so either 14 or 15). Optionally, enter 0 to ignore the nitrogen.*)
(*	"nitrogenHyperfineTensor"		The 3x3 nitrogen hyperfine tensor, or a list of the diagonal elements, or a string specifying the source.*)
(*	"nitrogenQuadrapolarTensor"	The 3x3 nitrogen quadrapolar tensor, or a list of the diagonal elements, or a string specifying the source.*)
(*	"carbonPositions "			A list of shell-index pairs, specifying which carbons have spin, e.g. {{1,1},{1,2},{1,3}} for all three closest carbons. *)
(*						The order of this list is important: it decides the order of the tensor product structure, and also the n's in the "carbonHyperfineTensors" parameter use this order.*)
(*	"carbonHyperfineTensors"		A list of the form {{Subscript[n, 1],Subscript[T, 1]},{Subscript[n, 2],Subscript[T, 2]},{Subscript[n, 3],Subscript[T, 3]},...} where the n's are integers refering to the *)
(*						carbon index number (see "carbonPositions"), and the T's are one of four things:*)
(*							- A 3x3 tensor matrix*)
(*							- A length 3 list of the diagonal elements of the tensor matrix*)
(*							- A string specifying the source to pull the tensor from*)
(*							- The string "dipolar", meaning the dipole-dipole tensor for the carbon-electron coupling is calculated for you based on the carbon position*)
(*	"zfsActivated"				True or False; whether the ZFS Hamiltonian is to be included.*)
(*	"nvZeemanActivated"		True or False; whether the NV Zeeman Hamiltonian is to be included.*)
(*	"nitrogenZeemanActivated"		True or False; whether the nitrogen Zeeman Hamiltonian is to be included.*)
(*	"carbonZeemanActivated"		True or False; whether the carbon Zeeman Hamiltonian is to be included.*)
(*	"nitogenQuadrapolarActivated"	True or False; whether the nitrogen quadrapolar Hamiltonian is to be included.*)
(*	"nitrogenHyperfineActivated"	True or False; whether the nitrogen hyperfine Hamiltonian is to be included.*)
(*	"carbonHyperfineActivated"		True or False; whether the carbon hyperfine Hamiltonian is to be included.*)
(*	"nitrogenCarbonDipolarActivated"	True or False; whether the nitrogen-carbon dipolar Hamiltonian is to be included.*)
(*	"carbonCarbonDipolarActivated"	True or False; whether the carbon-carbon dipolar Hamiltonian is to be included.*)
(*	"analyticForm" 			True or False; true to avoid numerical values*)
(*	"BCoords"				A coordinate system from the VectorAnalysis package, such as Spherical or Cartesian*)
(*ALL of the above parameters are optional; it will pick default values for them if you don't specify. *)
(*Note: If you happen to both give a tensor and the source, it will choose to use the tensor you gave it instead of the tensor corresponding to the source.*)


NVHamiltonian[params_?ListQ]:=
NVHamiltonianExplicitForm[Sequence@@({"\[CapitalDelta]","B","crystalOrientation","nvOrientation",
"nitrogenIsotope","nitrogenHyperfineTensor","nitrogenQuadrapolarTensor",
"carbonPositions","carbonHyperfineTensors","zfsActivated",
"nvZeemanActivated","nitrogenZeemanActivated","carbonZeemanActivated",
"nitrogenQuadrapolarActivated","nitrogenHyperfineActivated","carbonHyperfineActivated",
"nitrogenCarbonDipolarActivated","carbonCarbonDipolarActivated",
"analyticForm","BCoords"}/.FillInNVParameters[params])];


(* ::Text:: *)
(*Below we make sequence inputs work for NVHamiltonian, e.g., NVHamiltonian["B"->{1,0,0}], will work, and the order of inputs doesn't matter.*)


NVHamiltonian[params___]:=NVHamiltonian[{params}];


End[];


(* ::Section::Closed:: *)
(*Other*)


(* ::Subsubsection:: *)
(*Usage Declarations*)


Secularize::usage = "Secularize[H] secularizes the Hamiltonian H with respect to the ZFS term";


$NVBasis::usage = "An orthogonal basis of operators for u(3), which includes the operators commonly used in NV Math, and then some extras: Sx, Sy, Sz, Sz.Sz, Syp, Sxp, projection onto ms=0, and two double quantum operators.";
$NVBasisLabels::usage = "A list of strings to go along with NVBasis.";


NVForm::usage = "NVForm[M] prints the matrix M in the basis you care about, with the basis elements written as Strings.";
NVMatForm::usage = "NVForm[M] prints the matrix M in the basis you care about, with the basis elements written with MatrixForm.";


(* ::Subsubsection:: *)
(*Implementations*)


Begin["`Private`"];


$NVBasis={Sx,Sy,Sz,Sz.Sz,Syp,Sxp,S0,Sxx,Sxy};
$NVBasisLabels={"Sx","Sy","Sz","\!\(\*SuperscriptBox[\(Sz\), \(2\)]\)","Sy'","Sx'","S0","Sxx","Sxy"};


NVForm[M_]:=
	Block[
		{n=Length[M],nc,basis,labels,printer},
		If[
			Mod[n,9]==0,
				nc=Log[2,n/9];
				basis=MultipartiteBasis[Join[{$NVBasis,$NVBasis},Array[{Subscript[\[DoubleStruckOne], 2],X,Y,Z}&,nc]]];
				labels=MultipartiteBasis[Join[{$NVBasisLabels,$NVBasisLabels},Array[{"I","X","Y","Z"}&,nc]]];,
				nc=Log[2,n/3];
				basis=MultipartiteBasis[Join[{$NVBasis},Array[{Subscript[\[DoubleStruckOne], 2],X,Y,Z}&,nc]]];
				labels=MultipartiteBasis[Join[{$NVBasisLabels},Array[{"I","X","Y","Z"}&,nc]]];
		];
		(Simplify@Chop@ExpandInOrthogonalBasis[M,basis]).(TableForm[{{Style[#,Bold]}}]&/@labels)
	]


NVMatForm[M_]:=
	Block[
		{n=Length[M],nc,basis},
		If[
			Mod[n,9]==0,
				nc=Log[2,n/9];
				basis=MultipartiteBasis[Join[{$NVBasis,$NVBasis},Array[{Subscript[\[DoubleStruckOne], 2],X,Y,Z}&,nc]]];,
				nc=Log[2,n/3];
				basis=MultipartiteBasis[Join[{$NVBasis},Array[{Subscript[\[DoubleStruckOne], 2],X,Y,Z}&,nc]]];
		];
		(Simplify@Chop@ExpandInOrthogonalBasis[M,basis]).(MatrixForm/@basis)
	]


Secularize[H_]:=
Block[{Ui,Uic,t,dim,test},
dim=Length[H]/3;
Ui=DiagonalMatrix[Flatten[Outer[Times,{t,1,t} ,Table[1,{dim}]]]]; (*so t=Exp[-it\[CapitalDelta]], but there is no point in actually doing the exponential*)
Uic=DiagonalMatrix[Flatten[Outer[Times,{1/t,1,1/t} ,Table[1,{dim}]]]];
test=FreeQ[#,t]&;
Map[ (*Now we search through every element of H and remove any expression that contains the symbol t*)
If[Head[#]===Plus,
Select[#,test],
If[test[#],#,0]
]&,
Expand[Ui.H.Uic],{2}
]
]


End[];


(* ::Section::Closed:: *)
(*Parameter Visualization*)


(* ::Subsection:: *)
(*Usage Declarations*)


VisualizeNVParameters::usage = "VisualizeNVParameters[params] outputs a visual display of the NV parameters. params is the same list of rules that you plug into NVHamiltonian[].";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


cylinderRadius=0.01;
sphereRadius=0.1;
activeOutline[scale_,pos_]:={Glow[Red], Black, FaceForm[Opacity[0.2]], Sphere[pos,1.2*scale*sphereRadius]};


CarbonTooltipString[shell_,index_]:="Shell:  "<>ToString[shell]<>"\nIndex: "<>ToString[index];


(* ::Text:: *)
(*Draws the nitrogen and vacancy spheres, as well as the line between them. Highlight the nitrogen if isotope is 15 or 14.*)


NVGraphics[rotation_,scale_,nitrogenIsotope_]:=
With[{R=scale*rotation,nitrogenLabel="Nitrogen "<>Piecewise[{{"(off)",nitrogenIsotope==0}}, ToString[nitrogenIsotope]]},{
Tooltip[{Blue,Sphere[R.$nitrogenPosition,scale*sphereRadius]},nitrogenLabel],
If[nitrogenIsotope==0,{},Tooltip[activeOutline[scale,R.$nitrogenPosition],nitrogenLabel]],
Tooltip[{White,Sphere[R.$vacancyPosition,scale*sphereRadius]},"Vacancy"],
{Green,Cylinder[{R.$vacancyPosition,R.$nitrogenPosition},scale*cylinderRadius]}
}]


(* ::Text:: *)
(*Draws all of the carbons up to maxshell, and the lines between them. Highlights the ones that are in the "positions" input.*)


CarbonGraphics[rotation_,scale_,maxshell_,positions_]:=
With[{R=scale*rotation},
Table[{
Tooltip[{GrayLevel[0.3],Sphere[R.CarbonPositions[{shell,index}],scale*sphereRadius]},CarbonTooltipString[shell,index]],
If[MemberQ[positions,{shell,index}],Tooltip[activeOutline[scale,R.CarbonPositions[{shell,index}]],CarbonTooltipString[shell,index]],{}],
{Green,Cylinder[{R.CarbonPositions[{shell,index}],R.(CarbonPositions[{shell,index}]-Last@CarbonDirections[{shell,index}])},scale*cylinderRadius]}
},{shell,maxshell},{index,ShellSize[shell]}]
]


(* ::Text:: *)
(*Draws a red wire frame.*)


CrystalBoxGraphics[rotation_,scale_]:=With[{
top1={2,4,1},top2={-2,4,1},top3={-2,-4,1},top4={2,-4,1},
bot1={2,4,-1},bot2={-2,4,-1},bot3={-2,-4,-1},bot4={2,-4,-1},
R=0.3*scale*rotation,
c1=Orange,c2=Lighter[Orange,0.6],c3=Red,c4=Lighter[Red,0.6]
},
{Thick,
c1,Line[{R.top1,R.top2}],
c2,Line[{R.top2,R.top3}],
c3,Line[{R.top3,R.top4}],
c4,Line[{R.top4,R.top1}],

c1,Line[{R.bot1,R.bot2}],
c2,Line[{R.bot2,R.bot3}],
c3,Line[{R.bot3,R.bot4}],
c4,Line[{R.bot4,R.bot1}],

c1,Line[{R.bot1,R.top1}],
c2,Line[{R.bot2,R.top2}],
c3,Line[{R.bot3,R.top3}],
c4,Line[{R.bot4,R.top4}]}
]


(* ::Text:: *)
(*Draws a red arrow.*)


BFieldGraphics[B_,BCoords_,scale_]:={Red,Arrowheads[0.03],Arrow[Tube[{{0,0,0},CoordinatesToCartesian[{If[B[[1]]==0,0,scale],B[[2]],B[[3]]},BCoords]},0.01]]}


(* ::Text:: *)
(*Concatinate all of the above functions using inputs from "params".*)


VisualizeNVParameters[params_?ListQ]:=Module[{p=FillInNVParameters[params],crystalRot,totalRot,maxShell},
crystalRot=RotMatZYZ[Sequence@@("crystalOrientation"/.p)];
totalRot=crystalRot.RotMatZYZ[Sequence@@NVAngles["nvOrientation"/.p]];
maxShell=Max[Part["carbonPositions"/.p,All,1],1];
Graphics3D[{
NVGraphics[totalRot,1,"nitrogenIsotope"/.p],
CarbonGraphics[totalRot,1,maxShell,"carbonPositions"/.p],
CrystalBoxGraphics[crystalRot,1],
BFieldGraphics["B"/.p,"BCoords"/.p,1]
},
PlotRange->1.1*Max[1.2,Sqrt[Total[(1.0*#)^2]]&/@Flatten[Table[CarbonPositions[{sh,ind}],{sh,maxShell}, {ind, ShellSize[sh]}], 1]],
Axes->True,
Ticks->None,
AxesLabel->{"x","y","z"},
ImageSize->{500,500},
PlotLabel->"NV Paramater Visualization\n The black frame is the lab axis, and the red frame is the crystal axis."
]]


(* ::Text:: *)
(*Below we make sequence inputs to VisualizeNVParameters work, as we did for NVHamiltonian.*)


VisualizeNVParameters[params___]:=VisualizeNVParameters[{params}];


End[];


EndPackage[];
