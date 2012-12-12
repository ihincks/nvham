(* ::Package:: *)

BeginPackage["NVCalibrationTools`",{"VectorAnalysis`","NVHamiltonian`"}]


(* ::Section:: *)
(*Predicates*)


(* ::Section:: *)
(*Four NV Magnetometry*)


(* ::Subsection:: *)
(*Usage Declarations*)


FourNVZeemanSplittings::usage = "FourNVZeemanSplittings[B,BCoords,crystalOrientation] returns the splitting one would see in ESR due the magnetic field B with coordinate system BCoords. The parameters have the same input format as NVHamiltonian. The last two parameters are defaulted to Spherical and {0,0,0} if omitted.}";
FourNVTransformation::usage = "FourNVTransformation[crystalOrientation] returns the 4x3 matrix which converts the actual magnetic field (in Cartesian coordinates) to ESR Zeeman splittings four each of the four NV orientations.";
FourNVTransformationInv::usage = "FourNVTransformation[crystalOrientation] returns the 3x4 matrix which is the Moore Penrose Pseudoinverse of FourNVTransformation[crystalOrientation].";
FourNVAllPossibilities::usage = "FourNVAllPossibilities[data,BCoords,crystialOrientation] returns all 16 of the magnetic fields that spectralDifferences could come from given that we don't know the sign of any element of spectralDifferences.
- data is a length four list of the Zeeman splittings in Gauss taken from ESR spectra, in the order $ez, $e1, $e2, $e3 
- BCoords is the coordinate system you want the output to be in, default is Spherical
- crystalOrientation is the orientaion of the crystal with respect to the lab frame, default is {0,0,0}. See NVHamiltonian paramaters for details.
- returns a 16x3 matrix, where each row is a possibility
";
FourNVBestChoices::usage = "FourNVBestChoices[data,BCoords,crystalOrientation] returns the B fields which are most consistent with the data.";
FourNVBestChoice::usage = "FourNVBestChoice[data,BCoords,crystalOrientation] returns the first element of FourNVBestChoices[data,BCoords,crystalOrientation]";
FourNVRandomData::usage = "FourNVRandomData[B,\[Sigma],N,BCoords,crystalOrientation] returns N data points with noise taken from a Gaussian centered at B with standard devation \[Sigma]. \[Sigma] can either be a single number, or four different deviations, one for each direction.";
FourNVErrorBounds::usage = "FourNVErrorBounds[data,\[Sigma],N,BCoords,crystalOrientation] estimates the standard deviation of each component of the B field (in whatever coordinate system you specify) given your four data points. N is the number of Monte Carlo points to take.";


(* ::Subsection:: *)
(*Implementations*)


Begin["Private`"];


FourNVZeemanSplittings[B_,BCoords_:Spherical,crystalOrientation_:{0,0,0},method_:"secular"]:=
If[method=="secular",
Table[With[
{H=NVHamiltonian["B"->B,"BCoords"->BCoords,"crystalOrientation"->crystalOrientation,"nvOrientation"->n]},
H[[1,1]]-H[[3,3]]],{n,4}]*((1/$Gauss/(2*\[Gamma]e))/.Join[$units,$physicalConstants])
,
"Exact splitting method not implemented yet; the secular approximation is pretty good anyways."
]


FourNVTransformation[crystalOrientation_:{0,0,0}]:=(#/Sqrt[#.#])&/@((RotMatZYZ[Sequence@@crystalOrientation].#)&/@{$ez,$e1,$e2,$e3})


FourNVTransformationInv[crystalOrientation_:{0,0,0}]:=PseudoInverse[FourNVTransformation[crystalOrientation]]


FourNVAllPossibilities[data_,BCoords_:Spherical,crystalOrientation_:{0,0,0}]:=
(CoordinatesFromCartesian[FourNVTransformationInv[crystalOrientation].#,BCoords])&/@((Abs[spectralDifferences]*#)&/@Tuples[{1,-1},4])


(* ::Text:: *)
(*A helper function for FourNVBestChoices: searches for the minimum element of outputs and returnns the corresponding element of inputs*)


ArgMinimum[outputs_, inputs_] := Pick[inputs, #==Min[outputs] & /@outputs]


FourNVBestChoices[data_,BCoords_:Spherical,crystalOrientation_:{0,0,0}]:=
With[{allpos=(Abs[data]*#)&/@Tuples[{1,-1},4]},
ArgMinimum[(Norm[(FourNVTransformation[crystalOrientation].FourNVTransformationInv[crystalOrientation]-IdentityMatrix[4]).#])&/@allpos,
(CoordinatesFromCartesian[FourNVTransformationInv[crystalOrientation].#,BCoords])&/@allpos]]


FourNVBestChoice[data_,BCoords_:Spherical,crystalOrientation_:{0,0,0}]:=FourNVBestChoices[data,BCoords,crystalOrientation][[1]]


FourNVRandomData[B_,\[Sigma]_?NumberQ,N_,BCoords_:Spherical,crystalOrientation_:{0,0,0},method_:"secular"]:=
RandomVariate[NormalDistribution[Re[#],Re[\[Sigma]]],N]&/@FourNVZeemanSplittings[B,BCoords,crystalOrientation,method]


FourNVRandomData[B_,\[Sigma]_?ListQ,N_,BCoords_:Spherical,crystalOrientation_:{0,0,0},method_:"secular"]:=
Inner[RandomVariate[NormalDistribution[Re[#1],Re[#2]],N]&,FourNVZeemanSplittings[B,BCoords,crystalOrientation,method],\[Sigma],List]


FourNVErrorBounds[data_,\[Sigma]_,N_:10000,BCoords_:Spherical,crystalOrientation_:{0,0,0}]:=
Module[{bestchoice,monteCarloData,monteCarloBFields},
bestchoice=FourNVBestChoice[data,BCoords,crystalOrientation];
monteCarloData=FourNVRandomData[bestchoice,\[Sigma],N,BCoords,crystalOrientation];
monteCarloBFields=CoordinatesFromCartesian[#,BCoords]&/@Transpose[FourNVTransformationInv[crystalOrientation].monteCarloData];
StandardDeviation/@Transpose[monteCarloBFields]
]



End[];


(* ::Section::Closed:: *)
(*Magnetic Field From A Disk*)


(* ::Subsubsection:: *)
(*Usage Declarations*)


DiskComputeFieldAtPoint::usage = "ComputeFieldPoint[{r,\[Phi],z},m:1000,d:7.5,h:2.5] computes the magnetic field at cylindrical coordinates {r,\[Phi],z} due to a solid cylindrical magnet at the origin of diameter d and height h, with a magnetization factor m. The return value is also in cylindrical coordinates. {0,0,0} is returned for locations within the magnet.";
ComputeDiskFieldAndSave::usage = "ComputeDiskFieldAndSave[filename,{rmin,rmax,dr},{zmin,zmax,dz},m:1000,d:7.5,h:2.5,parallel:False] computes the field using DiskComputeFieldAtPoint on each point of the rectangular array defined by {rmin,rmax,dr},{zmin,zmax,dz}. The results are stored to file using Save2DField. Calculation will be done in parallel if you set the parallel flag to True.";
Save2DField::usage = "Save2DField[filename,{rmin,rmax},{zmin,zmax},magmat] saves the magnetic field data stored in magmat to disk, where rmin, rmax, zmin, and zmax indicate the physical extend of the data in millimetres.";
Load2DField::usage = "Load2DField[filename,M:1000] is the inverse of Save2DField, with the output format {{rmin,rmax},{zmin,zmax},M*magmat}.";
Interpolate2DField::usage = "Interpolate2DField[data] accepts the output of Load2DField, and interpolates the matrix to third order, returning an InterpolationFunciton.";
CylindricalFieldExtrapolation::usage = "FieldInterpolation[{x,y,z},fun2d] returns the field at position {x,y,z} of a cylindrically symmetric, radial field defined by the 2D vector field fun2d, which is a slice of the field through the z-x plane. You almost certainly want to set fun2d to the output of Interpolate2DField.";
Plot2DField::usage = "Plot2DField[fun2d] (OR Plot2DField[{fun2d1,fun2d2,...}] (doesn't work)) takes a InterpolatingFunction which is a 2D vector field and plots it. You probably want fun2d to be the output of Interpolate2DField.";
CompareSolutionFields::usage = "CompareSolutionFields[{z1,z2,z3},fun2d,{zz1,zz2,zz3},funn2d] takes the normalized inner product between solution translations {z1,z2,z3} and {zz1,zz2,zz3} over the possible volume swept by the magnet motors.";


(* ::Subsubsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*ff and gg are intermediate functions for DiskComputeFieldAtPoint*)


ff[z_,ZZ_,R_,r_,\[Phi]_]:=-((2 (z-ZZ) (R^2+(z-ZZ)^2-r R Cos[\[Phi]]))/(Sqrt[r^2+R^2+z^2-2 z ZZ+ZZ^2-2 r R Cos[\[Phi]]] (-R^2-2 (z-ZZ)^2+R^2 Cos[2 \[Phi]])))
gg[z_,ZZ_,R_,r_,\[Phi]_]:=( r Sqrt[r^2+R^2-2 r R Cos[\[Phi]]] Cos[ArcTan[R-r Cos[\[Phi]],r Sin[\[Phi]]]])/(r^2+R^2+z^2-2 z ZZ+ZZ^2-2 r R Cos[\[Phi]])^(3/2)


(* ::Text:: *)
(*The following is equivalent to doing the 3D integral over a solid cylinder of the dipole formula in cylindrical coordinates (without the Subscript[\[Mu], 0]/4\[Pi] factor). The functions ff and gg are used because some of the integrals are analytically solvable, so there is no need to waste time numerically doing them. ff and gg contain the jacobian, by the way. m=1000 is a decent fit to the magnet we have.*)


DiskComputeFieldAtPoint[{R_?NumberQ,\[CapitalPhi]_?NumberQ,ZZ_?NumberQ},m_:1000,d_:7.5,h_:2.5]:=
	Piecewise[{{{0,0,0},-d/2<=R<=d/2 && -h/2<=ZZ<=h/2}},
		Block[{rm,zm,\[Phi]m,r,\[Phi],z,v,\[Theta],Bdz,Bdr},
			rm=(r^2+R^2-2r R Cos[\[Phi]])^(1/2);
			zm=ZZ-z;
			\[Phi]m=ArcTan[R-r Cos[\[Phi]],r Sin[\[Phi]]];
			v=(rm^2+zm^2)^(1/2);
			\[Theta]=ArcTan[rm/zm];
			Bdz=2m NIntegrate[
				ff[h/2,ZZ,R,d/2,\[Phi]]-ff[h/2,ZZ,R,0,\[Phi]]
				-ff[-h/2,ZZ,R,d/2,\[Phi]]+ff[-h/2,ZZ,R,0,\[Phi]],
				{\[Phi],0,\[Pi]},AccuracyGoal->5];
			Bdr=2m NIntegrate[
				gg[h/2,ZZ,R,r,\[Phi]]-gg[-h/2,ZZ,R,r,\[Phi]],
				{\[Phi],0,\[Pi]},{r,0,d/2},AccuracyGoal->5];
			{Bdr ,\[CapitalPhi],Bdz}
		]
	]


ComputeDiskFieldAndSave[filename_,{rmin_,rmax_,dr_},{zmin_,zmax_,dz_},m_:1000,d_:7.5,h_:2.5,parallel_:False]:=
	Module[{data,tablefcn},
		tablefcn=If[parallel,ParallelTable,Table];
		data=Map[
				{#[[1]],#[[3]]}&,
				tablefcn[
					DiskComputeFieldAtPoint[{r,0,z},m,d,h],
					{z,zmin,zmax,dz},
					{r,rmin,rmax,dr}
				],
				{2}
			];
		Save2DField[
			filename,
			{rmin,rmax},
			{zmin,zmax},
			data
		]
	]


Save2DField[filename_,rbounds_,zbounds_,magmat_]:=
	Block[{rmin,rmax,zmin,zmax,mat},
		rmin=rbounds[[1]];
		rmax=rbounds[[2]];
		zmin=zbounds[[1]];
		zmax=zbounds[[2]];
		mat=magmat;
		Save[filename,{rmin,rmax,zmin,zmax,mat}];
		Print["Save Complete."];
	]


Load2DField[filename_String,M_]:=
	Block[{rmin,rmax,zmin,zmax,mat},
		Get[filename];
		{{rmin,rmax},{zmin,zmax},(M/1000)*mat}
	]


Interpolate2DField[data_]:=
	Block[{rmin,rmax,zmin,zmax,mat,W,H},
		rmin = data[[1,1]];
		rmax = data[[1,2]];
		zmin = data[[2,1]];
		zmax = data[[2,2]];
		H = Dimensions[data[[3]]][[1]];
		W = Dimensions[data[[3]]][[2]];
		Interpolation[
			Flatten[
				Table[
					{
						rmin+(rmax-rmin)(j-1)/(W-1), 
						zmin+(zmax-zmin)(i-1)/(H-1), 
						data[[3,i,j]]
					},
					{i,H},{j,W}
				],
			1]
		]
	]


(* ::Text:: *)
(*Note that we assume B[{x,y,-z}]=B[{x,y,z}]*{-1,-1,1} in the following function:*)


CylindricalFieldExtrapolation[{x_,y_,z_},intfun2d_]:=
	Block[{r,\[Phi],val,out},
		r=Sqrt[x^2+y^2];
		\[Phi]=If[x==0&&y==0,0,ArcTan[x,y]];
		val=intfun2d[r,Abs[z]];
		out=CoordinatesToCartesian[{val[[1]],\[Phi],val[[2]]},Cylindrical];
		If[z>=0,out,out*{-1,-1,1}]
	]


(* ::Text:: *)
(*Plot2DField with the list isn't working...*)


Plot2DField[fun2d_]:=
	With[{minr=fun2d[[1,1,1]],maxr=fun2d[[1,1,2]],minz=fun2d[[1,2,1]],maxz=fun2d[[1,2,2]]},
		StreamPlot[fun2d[r,z],{r,minr,maxr},{z,minz,maxz}]
	]
Plot2DField[fun2dlist_List]:=
	Block[{minr,maxr,minz,maxz,functions},
		minr=#[[1,1,1]]&/@fun2dlist;
		maxr=#[[1,1,2]]&/@fun2dlist;
		minz=#[[1,2,1]]&/@fun2dlist;
		maxz=#[[1,2,2]]&/@fun2dlist;
		functions=
			Function[{r,z},
				Evaluate[Piecewise[
					{{{0,0},(r<minr[[#]])||(r>maxr[[#]])||(z<minz[[#]])||(z>maxz[[#]])}},
					(fun2dlist[[#]])[r,z]
				]]
			]&/@{1,2};
		StreamPlot[
			(#[r,z])&/@functions,
			{r,Min[minr],Max[maxr]},
			{z,Min[minz],Max[maxz]}
		]
	]


CompareSolutionFields[{z1_,z2_,z3_},fun2d_,{zz1_,zz2_,zz3_},funn2d_]:=
	With[{
		a=NIntegrate[
			Hold[
				CylindricalFieldExtrapolation[MagToLab[{z1,z2,z3}-{x,y,z}],fun2d].
				CylindricalFieldExtrapolation[MagToLab[{zz1,zz2,zz3}-{x,y,z}],funn2d]
			],
			{x,0,25},
			{y,0,25},
			{z,10,25}
		],
		b=NIntegrate[
			Hold[Norm[CylindricalFieldExtrapolation[MagToLab[{z1,z2,z3}-{x,y,z}],fun2d]]^2],
			{x,0,25},
			{y,0,25},
			{z,10,25}
		],
		c=NIntegrate[
			Hold[Norm[CylindricalFieldExtrapolation[MagToLab[{zz1,zz2,zz3}-{x,y,z}],funn2d]]^2],
			{x,0,25},
			{y,0,25},
			{z,10,25}
		]
		},
			a/Sqrt[b*c]
	]
	


End[];


(* ::Section::Closed:: *)
(*Initializing Magnetic Fields*)


(* ::Subsubsection:: *)
(*Usage Declarations*)


Field::usage = "Field[{x,y,z},m:1] returns the magnetic field at the Cartesian point {x,y,z} in Gauss, and in Cartesian Coordinates, due to a fixed permanent magnet. The result is multiplied by m. To initialize the function, call LoadField.";
LoadField::usage = "LoadField[data,method] populates the Field[] function with a field. Method can be one of the following strings, which require different data formats:
\"CylindricalFromFile\" - data must be of the form {filename,m} where filename is the name of a file as stored by Save2DFile, and m is the corresponding magnetization";


(* ::Subsubsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Initialize the field to return an error if you try to call it.*)


Field[x___]:=(Print["An attempt to sample the magnetic field was made, but no magnetic field has been loaded. Call LoadField first."];Abort[];);


LoadField[data_,method_]:=
	Catch[
		Clear[Field];
		Which[
			method=="CylindricalFromFile",
				If[Not[MatchQ[data,{_String,_?NumericQ}]],Throw["Invalid input; data of the form {\"filename\",m} is expected."];];
				If[FindFile[data[[1]]]==$Failed,Throw["File "<>data[[1]]<>" not found."];];
				With[{fun2d=Interpolate2DField[Load2DField[Sequence@@data]]},
					Field[{x_,y_,z_},m_:1]:=m*CylindricalFieldExtrapolation[{x,y,z},fun2d];
				],
			True,
				Throw["method string not understood by LoadField."];
		];
		Print["Field Loaded Successfully. Access it using the Field[{x,y,z},m:1] function."];
	]


End[];


(* ::Section::Closed:: *)
(*Magnet Calibration From Position Array of Splitting Data	*)


(* ::Subsubsection:: *)
(*Usage Declarations*)


ImportSplittingData::usage = "ImportSplittingData[filename,scale:{1,1,1,1,1}] uses Get (so the file can be any number of formats, including .mat) to import data files where each row is a datapoint and each column is of the form {x,y,z,Splitting,nvOrientation}. Each datapoint is element-wise multiplied with scale.
ImportSplittingData[filename,nvOrientation,scale:{1,1,1,1,1}] uses Get (so the file can be any number of formats, including .mat) to import data files where each row is a datapoint and each column is of the form {x,y,z,Splitting} and sets each datapoint to have the orientation nvOrientation. Each datapoint is element-wise multiplied with scale.";
ExtractPositionArray::usage = "ExtractPositionArray[data] extracts just the positions from the data.";
CompareDataWithFourOrientations::usage = "CompareDataWithFourOrientations[data] calls PlotSplittingData on data, as well as on FakeData at the same positions, and one for each orientation."


MagToLab::usage = "MagToLab[{x,y,z}] converts the numbers on the magnet motor controller to labframe oriented coordinates (Note: without the translation, which is what we are trying to find).";
LabToMag::usage = "LabToMag[{x,y,z}] converts labframe oriented coordinates to the weird coordinates of the magnet motor (Note: without the translation, which is what we are trying to find).";
NVZeemanSplitting::usage = "NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},nvOrientation] computes the Zeeman splitting of the NV electron given the magnetic field, the crystal orientation, and the nv orientation. Uses the secular approximation..";


FakeData::usage = "FakeData[\[Sigma],points,nvOrientation_Integer,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m]
generates fake data of the form {{x1,y1,z1,splitting1,nvorientation1},...} by sampling from a Gaussian of variance \[Sigma] at each point of the list points={{x1,y1,z1},{x2,y2,z3},...}, with the given frame transformation z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3 and magnitization m.
FakeData[\[Sigma],points,nvOrientations_List,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m]
generates fake data of the form {{x1,y1,z1,splitting1,nvorientation1},...} by sampling from a Gaussian of variance \[Sigma] at each point of the list points={{x1,y1,z1},{x2,y2,z3},...} with respective nv orientations nvOrientations={nv1,nv2,...}, with the given frame transformation z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3 and magnetization m.
FakeData[\[Sigma],points,nvOrientations_List,solution] calls FakeData[\[Sigma],points,nvOrientations_List,Sequence@@solution]
FakeData, in all cases, returns an output which is compatible with the data input of both PlotSplittingData and CalibrateField.";
PlotSplittingData::usage = "PlotSplittingData[data,minsplit:0] 
shows you the splitting at each point in space in your data set. 
All data points which have a splitting less than minsplit will be Black, otherwise, each NV orientation gets a different colour.";


MostInformativeSpot::usage = "MostInformativeSpot[{z11,z12,z13},{\[Alpha]11,\[Alpha]12,\[Alpha]13},m1,{z21,z22,z23},{\[Alpha]21,\[Alpha]22,\[Alpha]23},m2,nvOrientation,minsplitting:15,maxsplitting:280,minz:6]
computes where you should do your next measurement to maximally distinguish (which means to maximize the difference of the splittings) the two solutions given by the two z's, \[Alpha]'s, and m's. You can specifiy the maximum and minimum splitting you can tolerate, as well as the lowest z3 value
MostInformativeSpot[solution1,solution2,nvOrientation,minsplitting:15,maxsplitting:280,minz:6] calls MostInformativeSpot[Sequence@@solution1,Sequence@@solution2,nvOrientation,minsplitting:15,maxsplitting:280,minz:6]";
MagnetFittingCostFunction::usage = "MagnetFittingCostFunction[{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m,data]
runs the cost function used by CalibrateField on the potential solution {z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3}, where m is the magnetization factor, and data is of the same form as accepted by CalibrateField.
MagnetFittingCostFunction[solution,data] returns MagnetFittingCostFunction[Sequence@@solution,data]";
CalibrateField::usage = "CalibrateField[data,m,constraints:{{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]}},method:\"NelderMead\",minsplit:13]
 does a numerical minimization to find the best frame transformation between the magnet origin and the NV PAS given the data.
The data should be a list of the form {{x1,y1,z1,S1,nv1},{x2,y2,z2,S2,nv2},...} where the x,y,z are the positions of the magnet motor, the S's are the splittings in MHz, and the nv's are the NV orientations. m is the magenetization factor, and minsplit specifies which data points to throw out.
The output of the function is of the form {mincostfunctionvalue,{{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m}} where mincostfunctionvalue is in units \!\(\*SuperscriptBox[\(MHz\), \(2\)]\), {z1,z2,z3} is the translation of the stage, {\[Alpha]1,\[Alpha]2,\[Alpha]3} are the Euler angles, and m is the magnetization (the same number you input).";
CalibrateFieldForBulk::usage = "CalibrateFieldForBulk is the same as CalibrateField, except that the cost function is also averaged over nv orientations.";


$typicalSolution::usage = "$typicalSolution is the sort of output you expect to get from CalibrateField. This is useful for feeding into functions like FakeData and PredictSplitting.";


PositionGivenField::usage = "PositionGivenField[{z1,z2,z3},m,B,BCoords:Cartesian,minz:6]
finds the position you should set the motors to in order to get the field B in coordinate system BCoords. You can demand that the z-motor be greater than 6.
The output is in the form {costfunctionvalue,{{x,y,z},m}}, where costfunctionvalue has units of MHz, {x,y,z} is where you should put the motors, and m is the magnetization ratio (the same number you input).
PositionGivenField[solution,B,BCoords:Cartesian,minz:6] returns PositionGivenField[solution[[1]],solution[[3]],B,BCoords:Cartesian,minz:6].";
AlignField::usage = "AlignField[{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m,nvOrientation,absB,minz:6] invokes PositionGivenField to find the motor position at which the magnetic field is aligned with the given nvOrienation, and with field strength Babs in Gauss.
AlignField[solution,nvOrientation,absB,minz:6] returns AlignField[Sequence@@solution,nvOrientation,absB,minz:6]";


PredictSplitting::usage = "PredictSplitting[{r1,r2,r3},nvOrientation,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m] 
predicts what the splitting will be at the magnet position {r1,r2,r3} given the the solution specified by z, \[Alpha] and m.
PredictSplitting[{r1,r2,r3},nvOrientation,solution] returns PredictSplitting[{r1,r2,r3},nvOrientation,Sequence@@solution].";
PredictNVPASVector::usage = "PredictNVPASVector[{r1,r2,r3},nvOrientation,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m] 
predicts what the magnetic field vector will be at the magnet position {r1,r2,r3} in the NV PAS given the the solution specified by z, \[Alpha] and m.
PredictNVPASVector[{r1,r2,r3},nvOrientation,solution] returns PredictNVPASVector[{r1,r2,r3},nvOrientation,Sequence@@solution].";


(* ::Subsubsection:: *)
(*Implementations*)


Begin["`Private`"];


$typicalSolution = {{1.34,23.5,-5.3},{0,-0.977,-0.04},.71};


ImportSplittingData[filename_String,scale_List:{1,1,1,1,1}]:=
	Catch[
		If[FindFile[filename]==$Failed,Throw["File "<>filename<>" not found."];];
		Block[{x=Flatten[Import[filename],1]},
			If[Dimensions[x][[2]]!=5,Throw["Your data doesn't appear to include NV orientations. Call ImportSplittingData[filename,nvOrientation] to specify all your orientations as equal."]];
			x=(#*scale)&/@x;
			ReplacePart[x\[Transpose],5->Round[x[[All,5]]]]\[Transpose]
		]
	]
ImportSplittingData[filename_String,nvOrientation_Integer,scale_List:{1,1,1,1,1}]:=
	Catch[
		If[FindFile[filename]==$Failed,Throw["File "<>filename<>" not found."];];
		Block[{x=Flatten[Import[filename],1]},
			x=Join[x[[All,{1,2,3,4}]]\[Transpose],{Table[nvOrientation,{k,Length[x]}]}]\[Transpose];
			(#*scale)&/@x
		]
	]


ExtractPositionArray[data_]:={#[[1]],#[[2]],#[[3]]}&/@data


CompareDataWithFourOrientations[data_]:=
	GraphicsRow[
		Join[
			{PlotSplittingData[data]},
			PlotSplittingData[#]&/@Table[FakeData[$MachineEpsilon,ExtractPositionArray[data],n,{1,23,-5},{0,-ArcCos[-1/3]/2,0},1],{n,4}]
		],
		ImageSize->1500
	]


(* ::Text:: *)
(*The following functions can be derived by going to the lab and looking at which ways the motors go. When facing the setup with your back to the computer, we want the lab frame to be the same as NVHamiltonian: z is up, x is to the right, and y is back towards the wall (right handed coordinates system).*)


MagToLab[{x_,y_,z_}]:={y,-x,-z};
LabToMag[{x_,y_,z_}]:={-y,x,-z};


(* ::Text:: *)
(*The following were found by running NVHamiltonian, taking the secular approximation, and finding the difference between the first and third diagonal component (so the following are a secular approximation to the splitting).*)


NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},1]=2*1000*Abs[With[{cz2=Cos[\[Theta]z2],sz2=Sin[\[Theta]z2],cy=Cos[\[Theta]y],sy=Sin[\[Theta]y]},0.0028025 Bz cy+0.0028025 Bx cz2 sy+0.0028025 By sy sz2]];
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},2]=2*1000*Abs[With[{cz1=Cos[\[Theta]z1],sz1=Sin[\[Theta]z1],cz2=Cos[\[Theta]z2],sz2=Sin[\[Theta]z2],cy=Cos[\[Theta]y],sy=Sin[\[Theta]y]},-0.00264222 Bz cz1 sy-0.000934165 Bx cz2 sy+0.00264222 By cz2 sz1-0.000934165 By sy sz2-0.00264222 Bx sz1 sz2+cy (-0.000934165 Bz+cz1 (0.00264222 Bx cz2+0.00264222 By sz2))]];
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},3]=2*1000*Abs[With[{cz1=Cos[\[Theta]z1],sz1=Sin[\[Theta]z1],cz2=Cos[\[Theta]z2],sz2=Sin[\[Theta]z2],cy=Cos[\[Theta]y],sy=Sin[\[Theta]y]},-0.000934165 Bx cz2 sy-0.00132111 By cz2 sz1+0.00228823 Bz sy sz1-0.000934165 By sy sz2+0.00132111 Bx sz1 sz2+cz1 (0.00228823 By cz2+0.00132111 Bz sy-0.00228823 Bx sz2)+cy (-0.000934165 Bz-0.00228823 Bx cz2 sz1-0.00228823 By sz1 sz2+cz1 (-0.00132111 Bx cz2-0.00132111 By sz2))]];
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},4]=2*1000*Abs[With[{cz1=Cos[\[Theta]z1],sz1=Sin[\[Theta]z1],cz2=Cos[\[Theta]z2],sz2=Sin[\[Theta]z2],cy=Cos[\[Theta]y],sy=Sin[\[Theta]y]},-0.000934165 Bx cz2 sy-0.00132111 By cz2 sz1-0.00228823 Bz sy sz1-0.000934165 By sy sz2+0.00132111 Bx sz1 sz2+cz1 (-0.00228823 By cz2+0.00132111 Bz sy+0.00228823 Bx sz2)+cy (-0.000934165 Bz+0.00228823 Bx cz2 sz1+0.00228823 By sz1 sz2+cz1 (-0.00132111 Bx cz2-0.00132111 By sz2))]];


NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},5]=NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},1];
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},6]=NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},2];
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},7]=NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},3];
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},8]=NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},4];


FakeData[\[Sigma]_,points_,nvOrientation_Integer,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},m_]:=
	With[{x=#[[1]],y=#[[2]],z=#[[3]]},
		{
			x,y,z,
			RandomVariate[NormalDistribution[
				NVZeemanSplitting[Field[MagToLab[{z1,z2,z3}-{x,y,z}],m],{\[Alpha]1,\[Alpha]2,\[Alpha]3},nvOrientation],\[Sigma]]],
			nvOrientation
		}
	]&/@points
FakeData[\[Sigma]_,points_,nvOrientations_List,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},m_]:=
	With[{x=#[[1]],y=#[[2]],z=#[[3]],nv=#[[4]]},
		{
			x,y,z,
			RandomVariate[NormalDistribution[
				NVZeemanSplitting[Field[MagToLab[{z1,z2,z3}-{x,y,z}],m],{\[Alpha]1,\[Alpha]2,\[Alpha]3},nv],\[Sigma]]],
			nv
		}
	]&/@(Join[points\[Transpose],{nvOrientations}]\[Transpose])
FakeData[\[Sigma]_,points_,nvOrientation_,solution_]:=FakeData[\[Sigma],points,nvOrientation,Sequence@@solution]


PlotSplittingData[data_,minsplit_:0]:=
	Block[{max,colours},
		max = Max[data[[All,4]]];
		colours[nn_]:=(Mod[nn-1,4]+1)/.{1->Blue,2->Green,3->Red,4->Cyan};
		Graphics3D[
			{
				If[#[[4]]>minsplit,colours[#[[5]]],Black],
				Sphere[{#[[1]],#[[2]],#[[3]]},#[[4]]/max]
			}&/@data,
			Axes->True,
			AxesLabel->{"x (mm)","y (mm)","z (mm)"},
			PlotRange->{0,25},
			PlotLabel->"Max Splitting: "<>ToString[max]<>"MHz",
			Epilog-> Inset[Framed[Column[Table[Style["Orientation "<>ToString[n], 12,colours[n]],{n,4}]], Background->LightYellow], {Right, Bottom}, {Right, Bottom}]]
	]


(* ::Text:: *)
(*Here is the cost function for the fitting. WARNING: it is duplicated in the CalibrationField for efficiency reasons. If you change one, be sure to change the other.*)


MagnetFittingCostFunction[{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},m_,data_]:=
	Sum[
			(NVZeemanSplitting[
				Field[MagToLab[{z1,z2,z3}-data[[k,{1,2,3}]]],m],
				{\[Alpha]1,\[Alpha]2,\[Alpha]3},
				data[[k,5]]
			]-data[[k,4]])^2
		,
		{k,1,Length[data]}
	]/Length[data];
MagnetFittingCostFunction[solution_,data_]:=MagnetFittingCostFunction[Sequence@@solution,data]


MostInformativeSpot[{z11_,z12_,z13_},{\[Alpha]11_,\[Alpha]12_,\[Alpha]13_},m1_,{z21_,z22_,z23_},{\[Alpha]21_,\[Alpha]22_,\[Alpha]23_},m2_,nvOrientation_,minsplitting_:15,maxsplitting_:280,minz_:6]:=
	Block[{r1,r2,r3,S1,S2,cost,result,values},
		S1[rr1_Real,rr2_Real,rr3_Real]:=NVZeemanSplitting[Field[MagToLab[{z11-rr1,z12-rr2,z13-rr3}],m1],{\[Alpha]11,\[Alpha]12,\[Alpha]13},nvOrientation];
		S2[rr1_Real,rr2_Real,rr3_Real]:=NVZeemanSplitting[Field[MagToLab[{z21-rr1,z22-rr2,z23-rr3}],m2],{\[Alpha]21,\[Alpha]22,\[Alpha]23},nvOrientation];
		cost[rr1_Real,rr2_Real,rr3_Real]:=(S1[rr1,rr2,rr3]-S2[rr1,rr2,rr3])^2;
		{result,values}=NMaximize[
			{cost[r1,r2,r3],0<=r1<=25,0<=r2<=25,minz<=r3<=25,minsplitting<S1[r1,r2,r3]<maxsplitting,minsplitting<S2[r1,r2,r3]<maxsplitting},
			{r1,r2,r3},
			MaxIterations->3000,
			Method->"NelderMead"
		];
		{r1,r2,r3}={r1,r2,r3}/.values;
		{Sqrt[result],{r1,r2,r3},S1[r1,r2,r3],S2[r1,r2,r3]}
	]
MostInformativeSpot[solution1_,solution2_,nvOrientation_,minsplitting_:15,maxsplitting_:280,minz_:6]:=MostInformativeSpot[Sequence@@solution1,Sequence@@solution2,nvOrientation,minsplitting,maxsplitting,minz]


CalibrateField[data_,m_,constraints_:{{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]}},method_:"NelderMead",minsplit_:13]:=
	Block[{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3,cost,sdata,cons,result,values},
		sdata=Select[data,(#[[4]]>=minsplit)&];
		cost[zz1_,zz2_,zz3_,\[Alpha]\[Alpha]1_,\[Alpha]\[Alpha]2_,\[Alpha]\[Alpha]3_]:=
			Sum[
				(
					NVZeemanSplitting[
						Field[MagToLab[{z1,z2,z3}-sdata[[k,{1,2,3}]]],m],
						{\[Alpha]1,\[Alpha]2,\[Alpha]3},sdata[[k,5]]]-sdata[[k,4]]
				)^2,
				{k,1,Length[sdata]}
			]/Length[sdata];
		{result,values}=NMinimize[
			{
				Hold[cost[z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3]],
				constraints[[1,1]]<=z1<=constraints[[1,2]],
				constraints[[2,1]]<=z2<=constraints[[2,2]],
				constraints[[3,1]]<=z3<=constraints[[3,2]],
				constraints[[4,1]]<=\[Alpha]1<=constraints[[4,2]],
				constraints[[5,1]]<=\[Alpha]2<=constraints[[5,2]],
				constraints[[6,1]]<=\[Alpha]3<=constraints[[6,2]]
			},
			{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3},
			Method->method,
			MaxIterations->3000
		];
		{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3}={z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3}/.values;
		{result,{{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m}}
	]


(* ::Text:: *)
(*This function assumes you have scanned your data along an axis where all four nv orientations give the same splitting.*)


CalibrateFieldForBulk[data_,m_,constraints_:{{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]}},method_:"NelderMead",minsplit_:13]:=
	Block[{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3,cost,sdata,cons,result,values},
		sdata=Select[data,(#[[4]]>=minsplit)&];
		cost[zz1_,zz2_,zz3_,\[Alpha]\[Alpha]1_,\[Alpha]\[Alpha]2_,\[Alpha]\[Alpha]3_]:=
			Sum[
				Sum[
					(
						NVZeemanSplitting[
							Field[MagToLab[{z1,z2,z3}-sdata[[k,{1,2,3}]]],m],
							{\[Alpha]1,\[Alpha]2,\[Alpha]3},nvOrientation]-sdata[[k,4]]
					)^2,
					{nvOrientation,4}
				]/4,
				{k,1,Length[sdata]}
			]/Length[sdata];
		{result,values}=NMinimize[
			{
				Hold[cost[z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3]],
				constraints[[1,1]]<=z1<=constraints[[1,2]],
				constraints[[2,1]]<=z2<=constraints[[2,2]],
				constraints[[3,1]]<=z3<=constraints[[3,2]],
				constraints[[4,1]]<=\[Alpha]1<=constraints[[4,2]],
				constraints[[5,1]]<=\[Alpha]2<=constraints[[5,2]],
				constraints[[6,1]]<=\[Alpha]3<=constraints[[6,2]]
			},
			{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3},
			Method->method,
			MaxIterations->3000
		];
		{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3}={z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3}/.values;
		{result,{{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},m}}
	]


PositionGivenField[{z1_,z2_,z3_},m_,B_,BCoords_:Cartesian,minz_:6]:=
	Block[{cost,Bcart,result,values,r1,r2,r3},
		Bcart=CoordinatesToCartesian[B,BCoords];
		cost[rr1_Real,rr2_Real,rr3_Real]:=Norm[Bcart-Field[MagToLab[{z1,z2,z3}-{rr1,rr2,rr3}],m]]^2;
		{result,values}=NMinimize[
			{
				cost[r1,r2,r3],
				0<=r1<=25,0<=r2<=25,minz<=r3<=25
			},
			{r1,r2,r3},
			MaxIterations->3000,
			Method->"NelderMead"
		];
		{Sqrt[result],{{r1,r2,r3}/.values,m}}
	]
PositionGivenField[solution_,B_,BCoords_:Cartesian,minz_:6]:=PositionGivenField[solution[[1]],solution[[3]],B,BCoords,minz]


AlignField[{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},m_,nvOrientation_,absB_,minz_:6]:=
	PositionGivenField[
		{z1,z2,z3},
		m,
		RotateBfromNVPAS[{0,0,-absB},{\[Alpha]1,\[Alpha]2,\[Alpha]3},nvOrientation,Cartesian],
		Cartesian,
		minz
	];
AlignField[solution_,nvOrientation_,absB_,minz_:6]:=AlignField[Sequence@@solution,nvOrientation,absB,minz]


PredictSplitting[{r1_,r2_,r3_},nvOrientation_,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},m_]:=
	NVZeemanSplitting[Field[MagToLab[{z1,z2,z3}-{r1,r2,r3}],m],{\[Alpha]1,\[Alpha]2,\[Alpha]3},nvOrientation]
PredictSplitting[{r1_,r2_,r3_},nvOrientation_,solution_]:=PredictSplitting[{r1,r2,r3},nvOrientation,Sequence@@solution]


PredictNVPASVector[{r1_,r2_,r3_},nvOrientation_,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},m_]:=
	RotateBtoNVPAS[CoordinatesFromCartesian[Field[MagToLab[{z1,z2,z3}-{r1,r2,r3}],m],Spherical],{\[Alpha]1,\[Alpha]2,\[Alpha]3},8,Spherical]
PredictNVPASVector[{r1_,r2_,r3_},nvOrientation_,solution_]:=PredictNVPASVector[{r1,r2,r3},nvOrientation,Sequence@@solution]


End[];


(* ::Section:: *)
(*Microwave Wire Functions*)


(* ::Subsubsection:: *)
(*Usage Declarations*)


FieldFromWire::usage="FieldFromWire[pos,end1,end2,current,distanceunit:$mm] calculates the magnetic field from a straight wire passing thruogh the coordinates end1 and end2 at the position pos in Gauss. Input the current in amps.";


DistanceFromWire::usage="DistanceFromWire[pos,end1,end2] calculates the closest distance between pos and the line in 3D determined by the points end1 and end 2.";
ClosestWirePoint::usage="ClosestWirePoint[pon,end1,end2] returns the coordinates of the closest point on the line defined by end1 and end2 to the position pos.";


(* ::Subsubsection:: *)
(*Implementations*)


Begin["`Private`"];


FieldFromWire[pos_,end1_,end2_,current_,distanceunit_:$mm]:=
	With[{dist=DistanceFromWire[pos,end1,end2]*distanceunit},
		(\[Mu]0/(2\[Pi]) current/dist 1/$Gauss)/.$units
	]


ClosestWirePoint[pos_,end1_,end2_]:=
	With[{diff=end2-end1},
		end1+(diff.(pos-end1))*diff/(diff.diff)
	]


DistanceFromWire[pos_,end1_,end2_]:=Norm[ClosestWirePoint[pos,end1,end2]-pos]


End[];


(* ::Section:: *)
(*Epilog*)


EndPackage[];
