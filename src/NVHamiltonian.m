(* ::Package:: *)

BeginPackage["NVHamiltonian`"];


(* ::Text:: *)
(*Todo:*)
(*-Change the static field frame if necessary*)
(*-Move the SpectrumData function here, or into QuantumUtils` from NVUtils`*)


(* ::Section::Closed:: *)
(*Ensure dependencies are met*)


(* ::Text:: *)
(*QuantumUtils` is a private work-a-day Mathematica package maintained by CoryLab for internal use. Many functions required by the NVSim code package are implemented in QuantumUtils`. Since many users of NVSim` at CoryLab will also want to use QuantumUtils` in the same Kernel session, the following logic is necessary to avoid writing multiple definitions of the same symbol in different contexts. Note that this logic is self contained, and will not appear in other sections of this package.*)


$hasQuantumUtils = Length[Position[$Packages,"QuantumUtils`"]]>0;
$hasNVUtils = Length[Position[$Packages,"NVUtils`"]]>0;


(* ::Text:: *)
(*Automatically import NVUtils` if neither of the *Utils` are around.*)


If[Not[$hasQuantumUtils||$hasNVUtils],
	Print["Warning: Neither QuantumUtils` nor NVUtils` was detected. Automatically importing NVUtils`."];
	Needs["NVUtils`"];
];


(* ::Text:: *)
(*Since the BeginPackage command makes all contexts inactive besides System` and NVHamiltonian`, we need to call Needs between the BeginPackage and EndPackage command. This is usually acheived using the second argument to BeginPackage, but we have special circumstances.*)


If[$hasQuantumUtils,Needs["QuantumUtils`"]];
If[$hasNVUtils,Needs["NVUtils`"]];


(* ::Text:: *)
(*The following unlucky functions are defined QuantumUtils`, but have different definitions in NVHamiltonian`. Therefore, it makes the most sense just to remove them here.*)


If[$hasQuantumUtils,
	With[{removeMe={QuantumUtils`ZeemanHamiltonian,QuantumUtils`HyperfineHamiltonian}},
		ClearAll/@removeMe;
		Remove/@removeMe;
	]
];


(* ::Section::Closed:: *)
(*Error Messages*)


(* ::Subsection:: *)
(*Frames and Vectors*)


Frame::notorthonormal = "The first three arguments of Frame must form an orthonormal basis for R^3.";


(* ::Subsection:: *)
(*Options*)


Numerical::badinput = "The option Numerical must be set to one of True, False, or Automatic.";


AngularUnits::badinput = "The option AngularUnits must be set to one of True, False, or Automatic.";


ZeroFieldSplitting::badinput = "The ZeroFieldSplitting must be a single number/symbol, or a list of length three (Each one will be multiplied by Sx.Sx, Sy.Sy, and Sz.Sz respectively.), or a list of two, {Dpar, Dperp}, Dpar Sz.Sz and Dperp (Sx.Sx-Sy.Sy).";


OutputFrame::badframe = "The output frame must be one of LabFrame, CrystalFrame, ZFSFrame, or ZeemanFrame.";


StaticField::badframe = "The static field's frame must be one of LabFrame, CrystalFrame, or ZFSFrame.";
StaticField::notvector = "The static field must be input as a Vector, e.g., Vector[{0,0,10}, Cartesian].";


(* ::Subsection:: *)
(*Hamiltonians*)


DipoleDipoleHamiltonian::equallocations = "The dipole-dipole Hamilotian was requested for two spins at the same physical locatian; this would result in division by 0. Instead, the 0 Hamiltonian has been return as the dipole-dipole Hamiltonian.";


(* ::Section::Closed:: *)
(*Predicates*)


(* ::Subsection:: *)
(*Usage Declarations*)


CarbonQ::usage = "CarbonQ[c] returns True iff the Head of c is Carbon";
NitrogenQ::usage = "NitrogenQ[c] returns True iff the Head of c is Carbon";
NucleusQ::usage = "CarbonQ[c] returns True iff one of CarbonQ or NitrogenQ holds.";


Vector3Q::usage = "Vector3Q[v] returns True iff v is a vector of length 3.";
Vector2Q::usage = "Vector2Q[v] returns True iff v is a vector of length 2.";
Matrix3Q::usage = "Matrix3Q[m] returns True iff m is a matrix of height 3.";


ValidReferenceFrameQ::usage = "ValidReferenceFrameQ[input] returns True iff input is one of LabFrame, CrystalFrame, ZFSFrame, or ZeemanFrame."


(* ::Text:: *)
(*We need to initialize the following things now so that they fall in the right context.*)


LabFrame::usage = "";ZFSFrame::usage = "";CrystalFrame::usage =  "";ZeemanFrame::usage = "";
Carbon::usage = ""; Nitrogen::usage = "";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


CarbonQ[c_]:=Head[c]===Carbon
NitrogenQ[c_]:=Head[c]===Nitrogen
NucleusQ[c_]:=CarbonQ[c]||NitrogenQ[c]


Vector2Q[v_]:=VectorQ[v]&&Length[v]==2
Vector3Q[v_]:=VectorQ[v]&&Length[v]==3
Matrix3Q[m_]:=MatrixQ[m]&&Length[m]==3


ValidReferenceFrameQ[input_]:=MemberQ[{LabFrame,ZFSFrame,CrystalFrame,ZeemanFrame},input]


End[];


(* ::Section:: *)
(*Physical Quantities and Spin*)


(* ::Subsection:: *)
(*Usage Declarations*)


\[Gamma]e::usage = "The gyromagnetic ratio of an electron. Units of Hz/G.";
\[Gamma]c::usage = "The gyromagnetic ratio of 13-Carbon. Units of Hz/G.";
\[Gamma]n15::usage = "The gyromagnetic ratio of 15-Nitrogen. Units of Hz/G."
\[Gamma]n14::usage = "The gyromagnetic ratio of 14-Nitrogen. Units of Hz/G."


\[Mu]0::usage = "The magnetic constant. Synonyms: the vacuum permeability, the permeability of free space. Units of 10^8*H/m. The multiplication by 10^8 is because we prefer to write gyromagnetic ratios in Hz/G instead of Hz/T.";
\[HBar]::usage = "Planck's reduced constant. Untis of J*s.";
\[Lambda]::usage = "The bond length in a diamond crystal. Units in metres.";


\[CapitalDelta]::usage = "The Zero Field Splitting (ZFS) of an NV center. Units of Hz.";


$constants::usage = "A list of replacement rules for the numerical values of physical constants.";
InsertConstants::usage = "InsertConstants[expr] replaces all physical constants appearing in the expression expr with their numerical value. Check the usage text of a given physical constant for the units used.";


SpinDim::usage = "SpinDim[halfInteger] returns the Hilbert space dimension of the given spin number";
IdentityInsert::usage = "IdentityInsert[C,dimA,dimB,n1,n2,n3] assumes C is a square matrix acting on a bipartite system with dimensions dimA and dimB and proceeds to add identity operations before, in-between, and after the bipartite system, with dimensions n1, n2, and n3, respectively.";


X::usage = "The 2x2 Pauli X operator.";
Y::usage = "The 2x2 Pauli Y operator.";
Z::usage = "The 2x2 Pauli Z operator.";


Si::usage="Spin-1 Identity Matrix";
Sx::usage="Spin-1 X Matrix.";
Sy::usage="Spin-1 Y Matrix.";
Sz::usage="Spin-1 Z Matrix.";
Syp::uage="Spin-1 Y Matrix in the interaction frame of Sz.Sz.";
Sxp::usage="Spin-1 X Matrix in the interaction frame of Sz.Sz.";
Sxx::usage="Spin-1 basis filler: gets the -1->1 transition";
Syy::usage="Spin-1 basis filler: gets the -1->1 transition";
S0::usage="Spin-1 basis filler: the projection onto ms=0";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


$constants = 
{
	\[CapitalDelta]     -> 2.87*10^9,
	\[Gamma]e    -> 2.802495266*10^6,
	\[Gamma]c    -> 1.0705*10^3,
	\[Gamma]n14  -> 0.3077*10^3,
	\[Gamma]n15  -> -.4316*10^3,
	\[Mu]0    -> 4\[Pi]*10^-7*10^8,
	\[Lambda]     -> 0.154*10^-9,
	\[HBar]     -> 1.054571726*10^\[Minus]34
};


SpinDim[s_?NumericQ]:=2*s+1


IdentityInsert[C_,dimA_,dimB_,n_]:=
	If[Length[C]!=dimA*dimB,
		Message[IdentityInsert::baddimensions];Abort;,
		ArrayFlatten[Map[IdentityMatrix[n]\[CircleTimes]#&,Partition[C,{dimB,dimB}],{2}]]
	];
IdentityInsert[C_,dimA_,dimB_,n1_,n2_,n3_]:=IdentityMatrix[n1]\[CircleTimes]IdentityInsert[C,dimA,dimB,n2]\[CircleTimes]IdentityMatrix[n3]


InsertConstants[expr_]:=expr/.$constants


X={{0,1},{1,0}};
Y={{0,-I},{I,0}};
Z={{1,0},{0,-1}};


Sx={{0,1,0},{1,0,1},{0,1,0}}/Sqrt[2];
Sy={{0,-I,0},{I,0,-I},{0,I,0}}/Sqrt[2];
Sz={{1,0,0},{0,0,0},{0,0,-1}};

Si=IdentityMatrix[3];
Syp={{0,-I,0},{I,0,I},{0,-I,0}}/Sqrt[2];
Sxp={{0,-1,0},{-1,0,1},{0,1,0}}/Sqrt[2];
Sxx={{0,0,1},{0,0,0},{1,0,0}};
Syy={{0,0,-I},{0,0,0},{I,0,0}};
S0={{0,0,0},{0,1,0},{0,0,0}};


End[];


(* ::Section::Closed:: *)
(*Visualization Tools*)


(* ::Subsection:: *)
(*Usage Declarations*)


NVForm::usage = "NVForm[H] displayes the NV Hamiltonian matrix in a convenient readable form.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


NVForm[H_,simplifier_:Simplify]:=Module[
	{
		hasNitrogen,
		numCarbon,
		mask,
		labels,
		basis,
		coeffs,
		i
	},
	(* Determine relevant dimensions *)
	If[Mod[Length@H,3]=!=0,Print["Invalid Hamiltonian dimension."];Abort[];];
	hasNitrogen = Mod[Length@H/3,3]===0;
	numCarbon = Log2@If[hasNitrogen, Length@H/9, Length@H/3];

	(* Construct basis and basis labels *)
	basis = {S0,Sz.Sz,Sx,Sy,Sz,Sxp,Syp,Sxx,Syy};
	labels = {"\!\(\*SubscriptBox[\(S\), \(0\)]\)","\!\(\*SuperscriptBox[SubscriptBox[\(S\), \(z\)], \(2\)]\)","\!\(\*SubscriptBox[\(S\), \(x\)]\)","\!\(\*SubscriptBox[\(S\), \(y\)]\)","\!\(\*SubscriptBox[\(S\), \(z\)]\)","\!\(\*SubscriptBox[\(S\), \(x\)]\)'","\!\(\*SubscriptBox[\(S\), \(y\)]\)'","\!\(\*SubscriptBox[\(S\), \(xx\)]\)","\!\(\*SubscriptBox[\(S\), \(yy\)]\)"};
	If[hasNitrogen, 
		basis = Flatten[Outer[KroneckerProduct, basis, basis, 1], 1];
		labels = Flatten[Outer[#1<>"\[CircleTimes]"<>#2&, labels, labels, 1], 1];
	];
	For[i=1, i<=numCarbon, i++,
		basis = Flatten[Outer[KroneckerProduct, basis, {Subscript[\[DoubleStruckOne], 2],X,Y,Z}, 1], 1];
		labels = Flatten[Outer[#1<>"\[CircleTimes]"<>#2&, labels, {"\[DoubleStruckOne]","X","Y","Z"}, 1], 1];
	];

	(* determine the coefficient of each basis element; assumes basis is orthogonal (which it is) *)
	coeffs = simplifier[(Tr[#\[ConjugateTranspose].H]/Tr[#\[ConjugateTranspose].#])&/@basis];

	(* only print non-zero coeffecients. *)
	mask=Not@*PossibleZeroQ/@coeffs;

	(* use a combination of MatrixForm and Row to get the desired form *)
	Inner[
		Row[{MatrixForm[{{#1}}],#2}]&,
		Pick[coeffs,mask],
		Style[#,Bold,Italic,16]&/@Pick[labels,mask],
		Plus
	]
]


End[];


(* ::Section::Closed:: *)
(*Frames and Vectors*)


(* ::Subsection:: *)
(*Usage Declarations*)


Coordinates::usage = "";
Cartesian::usage = "{x,y,z}";
Spherical::usage = "{r,\[Theta],\[Phi]} where \[Theta] is the azimuthal angle.";
Cylindrical::usage = "{\[Rho],\[Theta],z}";


Vector::usage = "";


SphericalToCartesian::usage = "SphericalToCartesian[{r,\[Theta],\[Phi]}] transforms the spherical vector {r,\[Theta],\[Phi]} into the cartesian vector {r Cos[\[Theta]]Sin[\[Phi]], r Sin[\[Theta]]Sin[\[Phi]], r Cos[\[Phi]]}. \[Phi] is the polar angle (down from the z axis) and \[Theta] is the azimual angle (from the x axis towards the y axis). We choose not to use the builtin function CoordinateTransform because it deals with 0 angles in a dumb way.";
CylindricalToCartesian::usage = "CylindricalToCartesian[{\[Rho],\[Theta],z}] transforms the cylindrical vector {\[Rho],\[Theta],z} into the cartesian vector {\[Rho] Cos[\[Theta]], \[Rho] Sin[\[Theta]], z}. \[Theta] is the azimual angle (from the x axis towards the y axis). We choose not to use the builtin function CoordinateTransform because it deals with 0 angles in a dumb way.";
CartesianToSpherical::usage = "CartesianToSpherical[{x,y,z}] transforms the cartesian vector {x,y,z} into the spherical vector {r,\[Theta],\[Phi]}. \[Phi] is the polar angle (down from the z axis) and \[Theta] is the azimual angle (from the x axis towards the y axis). If both x and y are 0, then the ambiguous azimuthal angle \[Theta] is taken to be 0.";
CylindricalToSpherical::usage = "CylindricalToSpherical[{\[Rho],\[Theta],z}] transforms the cylindrical vector {\[Rho],\[Theta],z} into the spherical vector {r,\[Theta],\[Phi]}. \[Phi] is the polar angle (down from the z axis) and \[Theta] is the azimual angle (from the x axis towards the y axis). If both \[Rho] and z are 0, then the ambiguous polar angle \[Phi] is taken to be 0.";
CartesianToCylindrical::usage = "CartesianToCylindrical[{x,y,z}] transforms the cartesian vector {x,y,z} into the cylindrical vector {\[Rho],\[Theta],z}. \[Theta] is the azimual angle (from the x axis towards the y axis). If both x and y are 0, then the ambiguous azimuthal angle \[Theta] is taken to be 0.";
SphericalToCylindrical::usage = "SphericalToCylindrical[{r,\[Theta],\[Phi]}] transforms the spherical vector {r,\[Theta],\[Phi]} into the cylindrical vector {\[Rho],\[Theta],z}.  is the polar angle (down from the z axis) and \[Theta] is the azimual angle (from the x axis towards the y axis).";
ChangeCoordinates::usage = "ChangeCoordinates[v, fromCoords, toCoords] changes the coordinates of the length-3 vector. For example, ChangeCoordinates[{1,1,0},Cartesian,Spherical].";


Frame::usage = "";


FrameMatrix::usage = "FrameMatrix[frame] returns the orthogonal 3x3 matrix, M, corresponding to the given frame. This is the matrix that transforms from the canonical basis to the basis of the frame, all in Cartesian coordinates (Frames entered in non-Cartesian cooridinates will be automatically converted). So, for example, M.{1,0,0}=x, where frame=Frame[x,y,z,Cartesian].";
FrameChangeMatrix::usage = "";


FrameInverse::usage = "FrameInverse[frame]";
FrameChange::usage = "FrameChange[M,fromFrame,toFrame] or FRame[v,fromFrame,ToFrame]";
FrameCompose::usage = "FrameCompose[framen,...,frame2,frame1] returns the resulting frame when all of the input frames are composed. That is, we know frame1 is written in the coordinates of the canonical basis, and if frame2 is written in the coordinates of frame1, and frame3 is written in the coordinates of frame2, etc, then the resulting frame is the composition of all frames written in the canonical coordinates.";


IdentityFrame::usage = "IdentityFrame";
BondFrame::usage = "BondFrame[{a,b,c},coords_]";
NVEulerAngles::usage = "NVEulerAngles[\[Theta]z1,\[Theta]y,\[Theta]z2] returns a Frame corresponding to rotating the IdenityFrame by the extrinsic ZYZ Euler angles \[Theta]z1,\[Theta]y,\[Theta]z2.";


PlotFrame::usage = "PlotFrame[frame1,frame2,...] plots each Frame given as an argement on the same figure.";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection::Closed:: *)
(*Coordinate Conversions*)


Begin["`Private`"];


(* ::Text:: *)
(*We define these conversion functions by hand so we don't have to deal with any opaque problems with the builtin function CoordinateTransform.*)


SphericalToCartesian[{r_,\[Theta]_,\[Phi]_}]:={r*Cos[\[Theta]]*Sin[\[Phi]],r*Sin[\[Theta]]*Sin[\[Phi]],r*Cos[\[Phi]]}
CylindricalToCartesian[{\[Rho]_,\[Theta]_,z_}]:={\[Rho]*Cos[\[Theta]],\[Rho]*Sin[\[Theta]],z}
CartesianToSpherical[{x_,y_,z_}]:=With[{r=Sqrt[x^2+y^2+z^2]},{r,If[PossibleZeroQ[x]&&PossibleZeroQ[y],0,ArcTan[x,y]],If[PossibleZeroQ[r],0,ArcCos[z/r]]}]
CylindricalToSpherical[{\[Rho]_,\[Theta]_,z_}]:={Sqrt[\[Rho]^2+z^2],\[Theta],If[PossibleZeroQ[\[Rho]]&&PossibleZeroQ[z],0,ArcTan[z,\[Rho]]]}
CartesianToCylindrical[{x_,y_,z_}]:={Sqrt[x^2+y^2],If[PossibleZeroQ[x]&&PossibleZeroQ[y],0,ArcTan[x,y]],z}
SphericalToCylindrical[{r_,\[Theta]_,\[Phi]_}]:={r*Sin[\[Phi]],\[Theta],r*Cos[\[Phi]]}


ChangeCoordinates[v_?Vector3Q,coords_,coords_]:=v
ChangeCoordinates[v_?Vector3Q,Spherical,Cartesian]:=SphericalToCartesian[v]
ChangeCoordinates[v_?Vector3Q,Cylindrical,Cartesian]:=CylindricalToCartesian[v]
ChangeCoordinates[v_?Vector3Q,Cartesian,Spherical]:=CartesianToSpherical[v]
ChangeCoordinates[v_?Vector3Q,Cylindrical,Spherical]:=CylindricalToSpherical[v]
ChangeCoordinates[v_?Vector3Q,Cartesian,Cylindrical]:=CartesianToCylindrical[v]
ChangeCoordinates[v_?Vector3Q,Spherical,Cylindrical]:=SphericalToCylindrical[v]


ChangeCoordinates[v_Vector,newCoords_]:=Vector[ChangeCoordinates[Value@v,Coordinates@v,newCoords],newCoords]


End[];


(* ::Subsubsection:: *)
(*Vectors*)


Begin["`Private`"];


(* ::Text:: *)
(*Functions to extract the vector and coordinates out of a Vector.*)


Vector/:Coordinates[Vector[_,coords_]]:=coords
Vector/:Value[Vector[v_,_]]:=v


(* ::Text:: *)
(*Define shorthand notation for converting coordinates.*)


Cartesian[v_Vector]:=ChangeCoordinates[v,Cartesian]
Spherical[v_Vector]:=ChangeCoordinates[v,Spherical]
Cylindrical[v_Vector]:=ChangeCoordinates[v,Cylindrical]


(* ::Text:: *)
(*To add Vectors, simply convert to Cartesian before doing usual addition. Similar for the dot product.*)


Vector/:Plus[v1_Vector,v2_Vector,rest___]:=Plus[Vector[Value@Cartesian@v1+Value@Cartesian@v2,Cartesian],rest]
Vector/:Dot[v1_Vector,v2_Vector]:=(Value@Cartesian@v1).(Value@Cartesian@v2);


Vector/:Times[\[Lambda]1___,v_Vector,\[Lambda]2___]:=Vector[Times[\[Lambda]1,\[Lambda]2,Value@Cartesian@v],Cartesian]


Vector/:Minus[v_Vector]:=Vector[-Value@Cartesian@v,Cartesian]


Vector/:Norm[v_Vector]:=
	Which[
		Coordinates@v===Spherical,
			First@Value@v,
		Coordinates@v===Cartesian,
			Norm[Value@v],
		Coordinates@v===Cylindrical,
			Sqrt[(First@Value@v)^2+(Last@Value@v)^2]
	]


End[];


(* ::Subsubsection:: *)
(*Frames*)


Begin["`Private`"];


(* ::Text:: *)
(*Frame is a head used in many plotting functions. Since it doesn't have any up or down values, it is safe to give it a double meaning. First, we unprotect it.*)


Unprotect@Frame;


(* ::Text:: *)
(*We assume a Cartesian coordinate system when none is specified.*)


Frame[x_,y_,z_]:=Frame[x,y,z,Cartesian]


(* ::Text:: *)
(*The columns of a 3x3 matrix are taken to be coordinates. The serves as a sort of inverse to the function FrameMatrix.*)


Frame[m_?Matrix3Q]:=Frame[m\[Transpose][[1]],m\[Transpose][[2]],m\[Transpose][[3]],Cartesian]


(* ::Text:: *)
(*Use a up values to transform a given frame to cartesian coordinates.*)


Frame/:Cartesian[f:Frame[a_,b_,c_,Cartesian]]:=f
Frame/:Cartesian[f:Frame[a_,b_,c_,Spherical]]:=Frame[SphericalToCartesian[a], SphericalToCartesian[b], SphericalToCartesian[c], Cartesian]
Frame/:Cartesian[f:Frame[a_,b_,c_,Cylindrical]]:=Frame[CylindricalToCartesian[a], CylindricalToCartesian[b], CylindricalToCartesian[c], Cartesian]


Frame/:Coordinates[Frame[_,_,_,coords_]]:=coords


(* ::Text:: *)
(*Returns the orthogonal matrix corresponding to the input frame.*)


FrameMatrix[f_Frame]:=((List@@Cartesian[f])[[1;;3]])\[Transpose]


End[];


(* ::Subsubsection:: *)
(*Changing, Inverting and Composing Frames*)


Begin["`Private`"];


(* ::Text:: *)
(*You can think of converting between frames as first converting to the canonical basis, and then from this converting to the desired frame.*)


FrameChangeMatrix[fromFrame_,toFrame_]:=FrameMatrix[toFrame]\[Transpose].FrameMatrix[fromFrame]


(* ::Text:: *)
(*We have the frame change matrix, so now it's just a matter of actually impementing it.*)


FrameChange[v_?Vector3Q,fromFrame_,toFrame_]:=FrameChangeMatrix[fromFrame,toFrame].v
FrameChange[M_?Matrix3Q,fromFrame_,toFrame_]:=With[{F=FrameChangeMatrix[fromFrame,toFrame]},F.M.F\[Transpose]]
FrameChange[v_Vector,fromFrame_,toFrame_]:=ChangeCoordinates[Vector[FrameChangeMatrix[fromFrame,toFrame].(Value@Cartesian@v),Cartesian],Coordinates@v]


(* ::Text:: *)
(*Since the matrices are orthogonal, the inverse is given by the transpose.*)


FrameInverse[f_Frame]:=Frame[FrameMatrix[f]\[Transpose]]


(* ::Text:: *)
(*Composition is of course just matrix multiplication in Cartesian coordinates.*)


FrameCompose[f_Frame]:=f
FrameCompose[fa_Frame,fb_Frame,rest___]:=FrameCompose[Frame[Simplify[FrameMatrix[fa].FrameMatrix[fb]]],rest]


End[];


(* ::Subsubsection:: *)
(*Special Frames*)


Begin["`Private`"];


(* ::Text:: *)
(*The canonical basis.*)


IdentityFrame=Frame[{1,0,0},{0,1,0},{0,0,1},Cartesian];


(* ::Text:: *)
(*We use the extrinsic ZYZ convention.*)


NVEulerAngles[\[Theta]z1_,\[Theta]y_,\[Theta]z2_]=Frame[RotationMatrix[\[Theta]z2, {0,0,1}].RotationMatrix[\[Theta]y, {0,1,0}].RotationMatrix[\[Theta]z1, {0,0,1}]];


(* ::Text:: *)
(*The bond frame is most easily describable in spherical coordinates, so convert first. Note that all "bond frame means" is some frame in which the z vector is parallel to the input vector of BondFrame; the x-y vectors are chosen sort of arbitrarily, but it doesn't matter because the tensor should be cylindrically symmetric.*)


BondFrame[v_Vector]:=With[{s=Value@ChangeCoordinates[v,Spherical]},NVEulerAngles[0,s[[3]],s[[2]]]]


End[];


(* ::Subsubsection::Closed:: *)
(*Plotting Frames*)


Begin["`Private`"];


(* ::Text:: *)
(*Just need to loop through each frame and draw three arrows.*)


PlotFrame[frames__]:=
	Graphics3D[
		MapIndexed[
			Module[
				{f=Cartesian[#1],x,y,z},
				x=f[[1]];
				y=f[[2]];
				z=f[[3]];
				{
					ColorData[1][First[#2]],
					Arrow[{{0,0,0},x}],Arrow[{{0,0,0},y}],Arrow[{{0,0,0},z}],
					Text["x",x*1.1],Text["y",y*1.1],Text["z",z*1.1]
				}
			]&,
			List[frames]
		],
		BoxRatios->{1,1,1},
		PlotRange->1.1*{{-1,1},{-1,1},{-1,1}}
	]


End[];


(* ::Section::Closed:: *)
(*NV and Lattice Geometry*)


(* ::Subsection:: *)
(*Usage Declarations*)


NVOrientationToFrame::usage = "NVOrientationToFrame[n] returns a Frame (with respect to the crystal frame, of course) corresponding to the n'th NV orientation. n should be one of the values 1,2,3,4,5,6,7,8. Here, 1 is along the positive z direction, 2 is on the x-z plane, rotated an angle ArcCos[-1/3] from the z axis, and orientations 3 and 4 are right-handed Z rotations of orientation 2 by angles 2\[Pi]/3 and 4\[Pi]/3 respectively. The orientations 5,6,7,8 are respectively anti-parallel to 1,2,3,4.";


E0::usage = "The 0th diamond lattice vector. This direction is parallel to the ZFS. This vector has unit length 1, and so represents a length equal to one bond length.";
E1::usage = "The 1st diamond lattice vector. This direction is rotated down from the ZFS axis by angle ArcCos[-1/3]. This vector has unit length 1, and so represents a length equal to one bond length.";
E2::usage = "The 2nd diamond lattice vector. This direction is rotated down from the ZFS axis by angle ArcCos[-1/3], and clockwise about z by 2\[Pi]/3. This vector has unit length 1, and so represents a length equal to one bond length.";
E3::usage = "The 3rd diamond lattice vector. This direction is rotated down from the ZFS axis by angle ArcCos[-1/3], and clockwise about z by 4\[Pi]/3. This vector has unit length 1, and so represents a length equal to one bond length.";


UnitCell::usage = "UnitCell is a list of Vectors specifying the location of the 18 atoms in a diamond lattice. The coordinates are normalized such that the distance between nearest atoms is 1. The width of the cell is therefore 4/Sqrt[3], and the diameter is 4.";
UnitCellWidth::usage = "The width of the cube containing the diamond unit cell, given that the distance between atoms is normalized to 1. Thus, this value is 4/Sqrt[3].";
UnitCellGraphic::usage = "UnitCellGraphic is a graphic of all atoms in the diamond unit cell. See also LaticePositionsGraphic.";
LatticePositions::usage = "LatticePositions[radius] returns a list of Vectors of the positions of all atoms in a diamond lattice no farther than radius from the origin. Radius is in units of nm, or, can contain the symbol \[Lambda], which is the distance between adjacent atoms in the lattice. If radius contains floating point, then the output is likewise numeric, otherwise, it is exact, and takes much longer.";
LatticePositionsGraphic::usage = "LatticePositionsGraphic[radius_,plotNV_:True,carbonOpacity_:1,OptionsPattern[Graphiccs3D]] plots all of the lattice positions from LatticePositions[radius]. See also UnitCellGraphic.";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*NV Orientations*)


Begin["`Private`"];


NVOrientationToFrame[1]=IdentityFrame;
NVOrientationToFrame[2]=NVEulerAngles[0,ArcCos[-1/3],0];
NVOrientationToFrame[3]=NVEulerAngles[0,ArcCos[-1/3],2\[Pi]/3];
NVOrientationToFrame[4]=NVEulerAngles[0,ArcCos[-1/3],4\[Pi]/3];
NVOrientationToFrame[5]=NVEulerAngles[0,\[Pi],0];
NVOrientationToFrame[6]=NVEulerAngles[0,ArcCos[-1/3]-\[Pi],0];
NVOrientationToFrame[7]=NVEulerAngles[0,ArcCos[-1/3]-\[Pi],2\[Pi]/3];
NVOrientationToFrame[8]=NVEulerAngles[0,ArcCos[-1/3]-\[Pi],4\[Pi]/3];


End[]


(* ::Subsubsection::Closed:: *)
(*Lattice Positions*)


Begin["`Private`"];


(* ::Text:: *)
(*The lattice directions.*)


E0=Vector[{0,0,1},Cartesian];
E1=Cartesian@Vector[{1,0,ArcCos[-1/3]},Spherical];
E2=Cartesian@Vector[{1,2\[Pi]/3,ArcCos[-1/3]},Spherical];
E3=Cartesian@Vector[{1,4\[Pi]/3,ArcCos[-1/3]},Spherical];


(* ::Text:: *)
(*The standard diamond lattice. It is rotated, such that the cube enclosing the cell is not square to cartesian coordinates.*)


UnitCell = {
	0*E0,E0,
	E0-E1,E0-E2,E0-E3,
	2*E0-E1, 2*E0-E2, 2*E0-E3,
	2*E0-2*E1, 2*E0-2*E2, 2*E0-2*E3,
	2*E0-E1-E2, 2*E0-E2-E3, 2*E0-E3-E1,
	E0-E2+E3-E1, E0-E3+E1-E2, E0-E1+E2-E3,
	3*E0-E1-E2-E3
};


(* ::Text:: *)
(*Evaluates to 4/Sqrt[3].*)


UnitCellWidth = Norm[E0-E2+E3-E1];


UnitCellGraphic := Module[{lines, positions},
	positions = Value/@UnitCell;
	lines = Flatten[Outer[If[Norm[#1-#2]==1,Tube[{##}],{}]&,positions,positions,1]];
	lines = Join[lines, Flatten[Outer[If[MemberQ[positions[[{1,18,17,16,15,11,10,9}]],#1]&&Norm[#1-#2]==4/Sqrt[3],Tube[{##}],{}]&,positions,positions,1]]];
	Graphics3D[{GrayLevel[.5],Tooltip[Sphere[#,0.2],Simplify@#]&/@positions,White,lines},Boxed->False]
]


(* ::Text:: *)
(*This implementation is rather computationally ineffecient.*)


LatticePositions[radius_] := Module[{positions, ncells, unitrad, numeric, uc, output,ex,ey,ez},
	numeric = Not[(radius/.x_Real:>I)===radius];
	unitrad = InsertConstants@Simplify[radius / \[Lambda]];
	ncells = Round[1.1*InsertConstants[unitrad / (UnitCellWidth)]];
	uc = If[numeric, N@UnitCell, UnitCell];
	(*these are the unit vectors in the coordinate frame of the cell*)
	ex = Vector[{-1/Sqrt[6],1/Sqrt[2],1/Sqrt[3]},Cartesian];
	ey=Vector[{Sqrt[2/3],0,1/Sqrt[3]},Cartesian];
	ez=Vector[{-1/Sqrt[6],-1/Sqrt[2],1/Sqrt[3]},Cartesian];
	output = Select[DeleteDuplicates[
		Flatten@Table[
			((UnitCellWidth*(x*ex+y*ey+z*ez)+#)&/@uc),
			{x,-ncells,ncells},
			{y,-ncells,ncells},
			{z,-ncells,ncells}
		]
	], Norm[Value@#]<=unitrad&];
	output
]


LatticePositionsGraphic[radius_,plotNV_:True,carbonOpacity_:1,opt:OptionsPattern[Graphics3D]]:=Module[{O={0,0,0},r=.3,tr=1/10,tc=Yellow,lp,tubelist},
	lp=Value/@LatticePositions[radius];
	tubelist=DeleteDuplicates[Select[Flatten[Outer[If[0<Norm[#1-#2]<=1.1,{#1,#2},None]&,lp,lp,1],1],#=!=None&],#1==#2||#1==Reverse[#2]&];
	If[plotNV,
		lp=Select[lp,Norm[#-{0,0,1}]>0&&Norm[#]>0&];
	];
	Graphics3D[
		{
			If[plotNV,
				{Yellow,Map[{Opacity[#[[1]]],Sphere[O,#[[2]]*r]}&,{{1,0.6},{0.6,0.8},{0.2,1},{0.1,1.2}},1]},
				{}
			],
			Opacity[1],Blue,Sphere[Value@E0,r],
			GrayLevel[.2],Opacity[carbonOpacity],Map[Sphere[#,r]&,lp],
			Yellow,Tube[#,tr]&/@tubelist
		},
		opt,
		Boxed->False
	]
]


End[];


(* ::Section::Closed:: *)
(*Options*)


(* ::Subsection:: *)
(*Usage Declarations*)


NVHamiltonian::usage = "";
NVOrientation::usage = "";
ZeroFieldSplitting::usage = "";
StaticField::usage = "";
StaticFieldFrame::usage = "";
NVSpin::usage = "";
OutputFrame::usage = "";
CrystalOrientation::usage = "";
Numerical::usage = "";
AngularUnits::usage = "";


LabFrame::usage = "";
ZFSFrame::usage = "";
CrystalFrame::usage =  "";
ZeemanFrame::usage = "";


CheckOptions::usage = "CheckOptions[opts] checks opts for any malformed/unacceptable inputs and Aborts if any are found. opt is an OptionsPattern for NVHamiltonian. This function is called inside of NVHamiltonian.";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*The Options*)


Begin["`Private`"];


Options[NVHamiltonian]=
	{
		ZeroFieldSplitting -> \[CapitalDelta],
		NVOrientation -> 1,
		StaticField -> Vector[{0,0,0},Cartesian],
		StaticFieldFrame -> LabFrame,
		NVSpin -> 1,
		OutputFrame -> ZFSFrame,		
		CrystalOrientation -> IdentityFrame,
		Numerical -> Automatic,
		AngularUnits -> Automatic
	};


(* ::Text:: *)
(*This loops through all of the options in NVHamiltonian and protects them.*)


Protect[Evaluate[Sequence@@(Options[NVHamiltonian][[All,1]])]];


End[];


(* ::Subsubsection:: *)
(*The Options Checker*)


Begin["`Private`"];


CheckOptions[OptionsPattern[NVHamiltonian]]:=
	Module[
		{abort},

		If[Not[MemberQ[{True,False,Automatic},OptionValue[Numerical]]],Message[Numerical::badinput];abort=True];
		If[Not[MemberQ[{True,False,Automatic},OptionValue[AngularUnits]]],Message[AngularUnits::badinput];abort=True];

		If[ListQ[OptionValue[ZeroFieldSplitting]]&&Not[Vector2Q[OptionValue[ZeroFieldSplitting]]||Vector3Q[OptionValue[ZeroFieldSplitting]]],Message[ZeroFieldSplitting::badinput];abort=True];

		If[Not[ValidReferenceFrameQ[OptionValue[OutputFrame]]],Message[OutputFrame::badframe];abort=True];

		If[Not[ValidReferenceFrameQ[OptionValue[StaticFieldFrame]]],Message[StaticField::badframe];abort=True];
		If[OptionValue[StaticFieldFrame]===ZeemanFrame,Message[StaticField::badframe];abort=True];
		If[Head[OptionValue[StaticField]]=!=Vector,Message[StaticField::notvector];abort=True];

		If[abort,Abort[];]
	]



End[];


(* ::Section:: *)
(*Nuclei*)


(* ::Subsection:: *)
(*Usage Declarations*)


Carbon::usage = "Carbon[{Apar,Aperp},vector]";
Nitrogen::usage = "Nucleus[isotope,{Apar,Aperp},quadrapolar]";


DipoleCarbon::usage = "DipoleCarbon[vector] returns a Carbon head, whose hyperfine interaction is just an electron-carbon13 dipole interaction.";


Tensor::usage = "";
Location::usage = "";
QuadrapoleTensor::usage = "";
Isotope::usage = "";
GyromagneticRatio::usage = "";
SpinValue::usage = "";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*Carbon*)


Begin["`Private`"];


Carbon/:Tensor[Carbon[{Apar_,Aperp_},v_Vector]]:=FrameChange[DiagonalMatrix[{Aperp,Aperp,Apar}],BondFrame[v],IdentityFrame]
Carbon/:Tensor[Carbon[{Axx_,Ayy_,Azz_},v_Vector]]:=FrameChange[DiagonalMatrix[{Axx,Ayy,Azz}],BondFrame[v],IdentityFrame]
Carbon/:Tensor[Carbon[A_?Matrix3Q]]:=A


Carbon/:Location[Carbon[_,v_Vector]]:=v
Carbon/:Location[Carbon[_]]:=Vector[{0,0,0},Cartesian]


Carbon/:Isotope[Carbon[___]]:=13


Carbon/:SpinValue[Carbon[___]]:=1/2


Carbon/:SpinDim[Carbon[___]]:=2


Carbon/:GyromagneticRatio[Carbon[___]]:=\[Gamma]c


End[];


(* ::Subsubsection:: *)
(*Dipole Carbon*)


Begin["`Private`"];


(* ::Text:: *)
(*We want the units to be in Hz; this is why the factor of 2\[Pi] appears. (Remember that the gyromagnetic ratios are not in angular units). As a check to make sure these numbers are right, we know that the dipolar coupling between two hydrogen atoms separated by .2nm is 15kHz (Levitt pg 212).*)


DipoleCarbon[vector_Vector]:=
	With[{R=\[Lambda] Norm[vector]},
		Carbon[2 \[Pi] ((-\[Mu]0 \[Gamma]e \[Gamma]c \[HBar])/(4 \[Pi] R^3)){2,-1},vector]
	]


End[];


(* ::Subsubsection:: *)
(*Nitrogen*)


Begin["`Private`"];


Nitrogen/:Tensor[Nitrogen[_,{Apar_,Aperp_},___]]:=DiagonalMatrix[{Aperp,Aperp,Apar}]


Nitrogen/:Location[Nitrogen[___]]:=E0


Nitrogen/:QuadrapoleTensor[Nitrogen[_,_,Q_?Matrix3Q]]:=Q
Nitrogen/:QuadrapoleTensor[Nitrogen[_,_,Q_?Vector3Q]]:=DiagonalMatrix[Q]
Nitrogen/:QuadrapoleTensor[Nitrogen[_,_,Q_]]:=DiagonalMatrix[{0,0,Q}]


Nitrogen/:Isotope[Nitrogen[isotope_,___]]:=isotope


Nitrogen/:SpinValue[Nitrogen[x___]]:=If[Isotope[Nitrogen[x]]===15,1/2,1]


Nitrogen/:SpinDim[n:Nitrogen[___]]:=SpinDim[SpinValue[n]]


Nitrogen/:GyromagneticRatio[Nitrogen[x___]]:=If[Isotope[Nitrogen[x]]===15,\[Gamma]n15,\[Gamma]n14]


End[];


(* ::Section:: *)
(*Nuclear Database*)


(* ::Text:: *)
(*In this section, for convenience, we catalog various measurements and calculations of nuclear hyperfine tensors from literature.*)


(* ::Subsection:: *)
(*Usage Declarations*)


NuclearDatabase::usage = "NVHamiltonian has included a number of measured and simulated values of hyperfine tensors of nuclei surrounding the NV defect in the diamond lattice. Each paper is given its own function which when called, returns either a Carbon or Nitrogen instance. The format of these function names is AuthorYearNucleus[...]. This function serves to remind the user of the names of these functions. Call NuclearDatabase[] to return a table of these names, and some more info.";


Gali08Nucleus::usage = 
"Gali08Nucleus[distance, index] returns the index'th nucleus at the specified distance from the vacancy, 
as calculated using ab initio supercell simulations. No numeric or alphanumeric labels were attached to
the various nuclei, so we reference them by their quoted distance to the vacancy, which happen to be 
unique in their simulation. Data is taken from Gali et al., Table II, \"Ab initio supercell calculations 
on nitrogen-vacancy center in diamond: Electronic structure and hyperfine tensors\", 
DOI: 10.1103/PhysRevB.77.155206.  Since only the eigenvalues of the tensors are provided in the paper, it 
is assumed that the eigenvectors lie on the line joining the nucleus and the vacancy. This will not in 
general be true. In the case where there is ambiguity about which nuclei correspond to which distance 
from the vacancy, a best guess was made (The 2.49, 3.86 and 5.00 carbons were chosen to be the ones 
below the equator instead of above, the carbons at 2.92 and 2.93 might be switched )."


Felton09Nucleus::usage = 
"Felton09Nucleus[label, index] returns the index'th nucleus of the given label, where the label is one of 
the following strings: \"14N\", \"15N\", \"Ca\", or \"Cg\". Data is taken from Felton et al., Tables II 
and III, \"Hyperfine interaction in the ground state of the negatively charged nitrogen vacancy center 
in diamond\", DOI: 10.1103/PhysRevB.79.075203. Since only the eigenvalues of the tensors are provided in 
the paper, it is assumed that the eigenvectors lie on the line joining the nucleus and the vacancy. This 
will not in general be true.";


Childress06Nucleus::usage =
"Childress06Nucleus[label] returns the carbon corresponding to the specified label, which can either be 
the string \"D\" or \"E\". These data is taken from Childress et al., Supplementary Online Material,
\"Coherent Dynamics of Coupled Electron and Nuclear Spin Qubits in Diamond\", DOI: 
10.1126/science.1131871.";


Shim13Nucleus::usage = 
"Shim13Nucleus[index] returns the index'th carbon, where index can be 1, 2, or 3; the three carbon sites 
closest to the vacancy. This data is taken from Shim et al., \"Characterization of hyperfine interaction 
between single electron and single nuclear spins in diamond assisted by quantum beat from the nuclear 
spin\", arXiv:1307.0257.";


Taminiau12Nucleus::usage = 
"Taminiau12Nucleus[label] returns the carbon of the specified label, where label can be any of the 
integers 1,2,3,4,5, or 6. These are distant carbons, and only their 0th order hyperfine parameters were 
measured, hence the carbon that is returned will have a hyperfine tensor with only the bottom row 
populated. These data are taken frome Taminiau et al., Table I, \"Detection and Control of Individual 
Nuclear Spins Using a Weakly Coupled Electron Spin\", DOI: 10.1103/PhysRevLett.109.137602.";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*Nuclear Database*)


Begin["`Private`"];


NuclearDatabase[] :=
	Grid[
		{
			{"Function Name", "Arguments", "Possible Label Values (Shell Count)", "Description", "Paper", "DOI"},
			{Gali08Nucleus, "label, index", "1.68 (1), 1.61 (3), 2.47 (6), 2.49 (3), 2.9 (6), 2.92 (3), 2.93 (3), 3.85 (6), 3.86 (3), 4.99 (6), 5.0 (3)", "Computer simulated hyperfine values", "Ab initio supercell calculations on nitrogen-vacancy center in diamond: Electronic structure and hyperfine tensors", "10.1103/PhysRevB.77.155206"},
			{Felton09Nucleus, "label, index", "\"14N\" (1), \"15N\" (1), \"Ca\" (3), \"Cg\" (6)", "ESR measurements of bulk sample", "Hyperfine interaction in the ground state of the negatively charged nitrogen vacancy center in diamond", "10.1103/PhysRevB.79.075203"},
			{Childress06Nucleus, "label", "\"D\", \"E\"", "Fits from ESEEM experiment. Arbitrary choices from fit degeneracy.", "Coherent Dynamics of Coupled Electron and Nuclear Spin Qubits in Diamond", "10.1126/science.1131871"},
			{Shim13Nucleus, "index", "(3)", "Precision measurement fitting to ground state shift due to off-axis magnet", "Characterization of hyperfine interaction between single electron and single nuclear spins in diamond assisted by quantum beat from the nuclear spin", "arXiv:1307.0257"},
			{Taminiau12Nucleus, "label", "1, 2, 3, 4, 5, 6", "Distant carbon spins, only 0th order pieces measured.", "Detection and Control of Individual Nuclear Spins Using a Weakly Coupled Electron Spin", "10.1103/PhysRevLett.109.137602"}
		},
		Frame->All,
		ItemStyle->Directive["Text", FontSize->10, ParagraphIndent->0, LineIndent->0],
		Alignment->Left,
		Spacings->{1,1},
		Background->{None,{LightGreen, None}}
	]


End[];


(* ::Subsubsection::Closed:: *)
(*Gali et al., 2008*)


Begin["`Private`"];


Gali08Nucleus[label_] := Gali08Nucleus[label, 1]


Gali08Nucleus[1.68, 1] = Nitrogen[14, 10^6*{-1.7,-1.7}, 0];


Gali08Nucleus[1.61, 1] = Carbon[10^6*{19.5, 110.2, 185.4}, E1];
Gali08Nucleus[1.61, 2] = Carbon[10^6*{19.5, 110.2, 185.4}, E2];
Gali08Nucleus[1.61, 3] = Carbon[10^6*{19.5, 110.2, 185.4}, E3];


Gali08Nucleus[2.47, 1] = Carbon[10^6*{-4.8,-3.7,-1.5}, E1-E2];
Gali08Nucleus[2.47, 2] = Carbon[10^6*{-4.8,-3.7,-1.5}, E1-E3];
Gali08Nucleus[2.47, 3] = Carbon[10^6*{-4.8,-3.7,-1.5}, E2-E1];
Gali08Nucleus[2.47, 4] = Carbon[10^6*{-4.8,-3.7,-1.5}, E2-E3];
Gali08Nucleus[2.47, 5] = Carbon[10^6*{-4.8,-3.7,-1.5}, E3-E1];
Gali08Nucleus[2.47, 6] = Carbon[10^6*{-4.8,-3.7,-1.5}, E3-E2];


Gali08Nucleus[2.49, 1] = Carbon[10^6*{-7.4,-7.3,-5.8}, E1-E0];
Gali08Nucleus[2.49, 2] = Carbon[10^6*{-7.4,-7.3,-5.8}, E2-E0];
Gali08Nucleus[2.49, 3] = Carbon[10^6*{-7.4,-7.3,-5.8}, E3-E0];


Gali08Nucleus[2.9, 1] = Carbon[10^6*{2.8,3.3,4.6}, E1-E2+E0];
Gali08Nucleus[2.9, 2] = Carbon[10^6*{2.8,3.3,4.6}, E1-E3+E0];
Gali08Nucleus[2.9, 3] = Carbon[10^6*{2.8,3.3,4.6}, E2-E1+E0];
Gali08Nucleus[2.9, 4] = Carbon[10^6*{2.8,3.3,4.6}, E2-E3+E0];
Gali08Nucleus[2.9, 5] = Carbon[10^6*{2.8,3.3,4.6}, E3-E1+E0];
Gali08Nucleus[2.9, 6] = Carbon[10^6*{2.8,3.3,4.6}, E3-E2+E0];


Gali08Nucleus[2.92, 1] = Carbon[10^6*{1.4,2.4,2.9}, E1-E2+E3];
Gali08Nucleus[2.92, 2] = Carbon[10^6*{1.4,2.4,2.9}, E2-E3+E1];
Gali08Nucleus[2.92, 3] = Carbon[10^6*{1.4,2.4,2.9}, E3-E1+E2];


Gali08Nucleus[2.93, 1] = Carbon[10^6*{3.4,4.7,4.9}, E1-E0+E2];
Gali08Nucleus[2.93, 2] = Carbon[10^6*{3.4,4.7,4.9}, E2-E0+E3];
Gali08Nucleus[2.93, 3] = Carbon[10^6*{3.4,4.7,4.9}, E3-E0+E1];


Gali08Nucleus[3.85, 1] = Carbon[10^6*{13.5,14.2,19.4}, 2*E1-E2];
Gali08Nucleus[3.85, 2] = Carbon[10^6*{13.5,14.2,19.4}, 2*E1-E3];
Gali08Nucleus[3.85, 3] = Carbon[10^6*{13.5,14.2,19.4}, 2*E2-E1];
Gali08Nucleus[3.85, 4] = Carbon[10^6*{13.5,14.2,19.4}, 2*E2-E3];
Gali08Nucleus[3.85, 5] = Carbon[10^6*{13.5,14.2,19.4}, 2*E3-E1];
Gali08Nucleus[3.85, 6] = Carbon[10^6*{13.5,14.2,19.4}, 2*E3-E2];


Gali08Nucleus[3.86, 1] = Carbon[10^6*{12.8,12.8,18}, 2*E1-E0];
Gali08Nucleus[3.86, 2] = Carbon[10^6*{12.8,12.8,18}, 2*E2-E0];
Gali08Nucleus[3.86, 3] = Carbon[10^6*{12.8,12.8,18}, 2*E3-E0];


Gali08Nucleus[4.99, 1] = Carbon[10^6*{2.6,2.7,3.8}, 2*E1-2*E2];
Gali08Nucleus[4.99, 2] = Carbon[10^6*{2.6,2.7,3.8}, 2*E1-2*E3];
Gali08Nucleus[4.99, 3] = Carbon[10^6*{2.6,2.7,3.8}, 2*E2-2*E1];
Gali08Nucleus[4.99, 4] = Carbon[10^6*{2.6,2.7,3.8}, 2*E2-2*E3];
Gali08Nucleus[4.99, 5] = Carbon[10^6*{2.6,2.7,3.8}, 2*E3-2*E1];
Gali08Nucleus[4.99, 6] = Carbon[10^6*{2.6,2.7,3.8}, 2*E3-2*E2];


Gali08Nucleus[5.0, 1] = Carbon[10^6*{1.5,1.5,2.2}, 2*E1-2*E0];
Gali08Nucleus[5.0, 2] = Carbon[10^6*{1.5,1.5,2.2}, 2*E2-2*E0];
Gali08Nucleus[5.0, 3] = Carbon[10^6*{1.5,1.5,2.2}, 2*E3-2*E0];


End[];


(* ::Subsubsection::Closed:: *)
(*Felton et al., 2009*)


Begin["`Private`"];


Felton09Nucleus[label_] := Felton09Nucleus[label, 1]


Felton09Nucleus["14N", 1] = Nitrogen[14, 10^6*{-2.14,-2.7}, 10^6*-5.01];


Felton09Nucleus["15N", 1] = Nitrogen[15, 10^6*{3.03,3.65}, 0];


Felton09Nucleus["Ca", 1] = Carbon[10^6*{198.2,120.8}, E1];
Felton09Nucleus["Ca", 2] = Carbon[10^6*{198.2,120.8}, E2];
Felton09Nucleus["Ca", 3] = Carbon[10^6*{198.2,120.8}, E3];


Felton09Nucleus["Cg", 1] = Carbon[10^6*{18.49,13.26}, 2*E1-E2];
Felton09Nucleus["Cg", 2] = Carbon[10^6*{18.49,13.26}, 2*E1-E3];
Felton09Nucleus["Cg", 3] = Carbon[10^6*{18.49,13.26}, 2*E2-E1];
Felton09Nucleus["Cg", 4] = Carbon[10^6*{18.49,13.26}, 2*E2-E3];
Felton09Nucleus["Cg", 5] = Carbon[10^6*{18.49,13.26}, 2*E3-E1];
Felton09Nucleus["Cg", 6] = Carbon[10^6*{18.49,13.26}, 2*E3-E2];


End[];


(* ::Subsubsection::Closed:: *)
(*Childress et al., 2006*)


Begin["`Private`"];


Childress06Nucleus["D"] = Carbon[10^6*{{0.4,-2.2,-2.1},{-2.2,2.6,-0.4},{-2.1,-0.4,3.5}}];


Childress06Nucleus["E"] = Carbon[10^6*{{5,-6.3,-2.9},{-6.3,4.2,-2.3},{-2.9,-2.3,8.2}}];


End[];


(* ::Subsubsection::Closed:: *)
(*Shim et al., 2013*)


Begin["`Private`"];


Shim13Nucleus[1] = Carbon[10^6*{30.3,122.9,226.6}, Vector[{1,0,123.5*\[Pi]/180}, Spherical]];
Shim13Nucleus[2] = Carbon[10^6*{30.3,122.9,226.6}, Vector[{1,2*\[Pi]/3,123.5*\[Pi]/180}, Spherical]];
Shim13Nucleus[3] = Carbon[10^6*{30.3,122.9,226.6}, Vector[{1,4*\[Pi]/3,123.5*\[Pi]/180}, Spherical]];


End[];


(* ::Subsubsection::Closed:: *)
(*Taminiau et al., 2012*)


(* ::Text:: *)
(*TODO: In this case it might make more sense to search the carbon lattice for the closest dipolar interaction, and then scale it to reproduce the same coupling.*)


Begin["`Private`"];


Taminiau12Nucleus[1] = Carbon[10^6*{{0,0,0},{0,0,0},{0.0838 Sin[21*\[Pi]/180],0,0.0838 Cos[21*\[Pi]/180]}}];


Taminiau12Nucleus[2] = Carbon[10^6*{{0,0,0},{0,0,0},{0.047 Sin[30*\[Pi]/180],0,0.047 Cos[30*\[Pi]/180]}}];


Taminiau12Nucleus[3] = Carbon[10^6*{{0,0,0},{0,0,0},{0.055 Sin[54*\[Pi]/180],0,0.055 Cos[54*\[Pi]/180]}}];


Taminiau12Nucleus[4] = Carbon[10^6*{{0,0,0},{0,0,0},{0.019 Sin[133*\[Pi]/180],0,0.019 Cos[133*\[Pi]/180]}}];


Taminiau12Nucleus[5] = Carbon[10^6*{{0,0,0},{0,0,0},{0.033 Sin[132*\[Pi]/180],0,0.033 Cos[132*\[Pi]/180]}}];


Taminiau12Nucleus[6] = Carbon[10^6*{{0,0,0},{0,0,0},{0.0251 Sin[51*\[Pi]/180],0,0.0251 Cos[51*\[Pi]/180]}}];


End[];


(* ::Section:: *)
(*Hamiltonians*)


(* ::Subsection:: *)
(*Usage Declarations*)


ClearAll[ZeemanHamiltonian,HyperfineHamiltonian];


ZFSHamiltonian::usage = "ZFSHamiltonian[{Dx_,Dy_,Dz_},nvSpin_], ZFSHamiltonian[{Dpar, Dperp},nvSpin_], or ZFSHamiltonian[D,nvSpin_]";


ZeemanHamiltonian::usage = "";


HyperfineHamiltonian::usage = "HyperfineHamiltonian[spin1,spin2,A]";


DipoleDipoleHamiltonian::usage = "";


QuadrapolarHamiltonian::usage = "";


NVHamiltonian::usage = "";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*ZFS Hamlitonian*)


Begin["`Private`"];


ZFSHamiltonian[{Dx_,Dy_,Dz_},nvSpin_]:= Dx Spin[1][nvSpin].Spin[1][nvSpin] + Dy Spin[2][nvSpin].Spin[2][nvSpin] + Dz Spin[3][nvSpin].Spin[3][nvSpin];
ZFSHamiltonian[{Dpar_,Dperp_},nvSpin_]:= Dperp( Spin[1][nvSpin].Spin[1][nvSpin] - Spin[2][nvSpin].Spin[2][nvSpin]) + Dpar Spin[3][nvSpin].Spin[3][nvSpin];
ZFSHamiltonian[D_,nvSpin_]:= D Spin[3][nvSpin].Spin[3][nvSpin];


End[];


(* ::Subsubsection:: *)
(*Zeeman Hamiltonian*)


Begin["`Private`"];


ZeemanHamiltonian[\[Mu]_,{Bx_,By_,Bz_},spin_] := \[Mu](Bx Spin[1][spin] + By Spin[2][spin] + Bz Spin[3][spin])


End[];


(* ::Subsubsection:: *)
(*Hyperfine Hamiltonian*)


Begin["`Private`"];


(* Need with-injection to avoid Spin HoldForm craziness *)
HyperfineHamiltonian[spin1_,spin2_,A_]:= Sum[With[{i=ii,j=jj},A[[i,j]]*(Spin[i][spin1]\[CircleTimes]Spin[j][spin2])],{ii,3},{jj,3}]


End[];


(* ::Subsubsection:: *)
(*Dipole Dipole Hamiltonian*)


Begin["`Private`"];


DipoleDipoleHamiltonian[spin1_,spin2_,\[Gamma]1_,\[Gamma]2_,v1_Vector,v2_Vector]:=
	With[{R=Sqrt[(v1-v2).(v1-v2)],e=Value@(v1-v2)},
		If[PossibleZeroQ[R],
			Message[DipoleDipoleHamiltonian::equallocations];
			ConstantArray[0,{SpinDim[spin1]*SpinDim[spin2],SpinDim[spin1]*SpinDim[spin2]}],
			(2 \[Pi] (-\[Mu]0 \[Gamma]1 \[Gamma]2 \[HBar])/(4\[Pi] R^3))*(3(Total@(Array[Spin[#][spin1]&,3]*e))\[CircleTimes](Total@(Array[Spin[#][spin2]&,3]*e))/R^2-Spin[1][spin1]\[CircleTimes]Spin[1][spin2]-Spin[2][spin1]\[CircleTimes]Spin[2][spin2]-Spin[3][spin1]\[CircleTimes]Spin[3][spin2])
		]
	]


End[];


(* ::Subsubsection:: *)
(*Quadrapolar Hamiltonian*)


Begin["`Private`"];


QuadrapolarHamiltonian[spin_,A_?Matrix3Q]:=Sum[With[{i=ii,j=jj},A[[i,j]]*Spin[i][spin].Spin[j][spin]],{ii,3},{jj,3}];


End[];


(* ::Subsubsection:: *)
(*Combined Hamiltonian*)


Begin["`Private`"];


NVHamiltonian[nuclei___,opt:OptionsPattern[]]:=
	Module[
		{
			cartesianB,sphericalB,
			dimNV,dimN,hasN,dimC,numC,dimList,
			nucleiList,
			nvSpin,
			termList,Term,ExpandTerm,HyperfineTerm,NuclearZeemanTerm,DipoleDipoleTerm,
			labFrame,crystalFrame,zfsFrame,zeemanFrame,frames,
			angularUnits,preFactor,numerical
		},

		(* Check the options for any malformed/incorrect/etc inputs *)
		CheckOptions[opt];

		(* Convert the B field to various coordinate systems *)
		cartesianB = Value@ChangeCoordinates[OptionValue[StaticField],Cartesian];
		sphericalB = Value@ChangeCoordinates[OptionValue[StaticField],Spherical];

		(* We will need to be able to convert the placeholder frame names to actual Frames *)
		frames:={LabFrame->labFrame, CrystalFrame->crystalFrame, ZFSFrame->zfsFrame, ZeemanFrame->zeemanFrame};

		(* Deal with all of the frames first (except nuclei frames) *)
		(* The idea is simple: set the OutputFrame to be the canonical basis, and compose all other frames with respect to this. *)
		(* There is the natural ordering LabFrame > CrystalFrame > ZFSFrame. The ZeemanFrame is sort of independent of this 
		higherarchy, and is determined by the option StaticFieldFrame *)
		Which[
			OptionValue[OutputFrame]===LabFrame,
				labFrame = IdentityFrame;
				crystalFrame = OptionValue[CrystalOrientation];
				zfsFrame = FrameCompose[NVOrientationToFrame[OptionValue[NVOrientation]],crystalFrame];
				zeemanFrame = FrameCompose[NVEulerAngles[0,sphericalB[[3]],sphericalB[[2]]],OptionValue[StaticFieldFrame]/.frames];,
			OptionValue[OutputFrame]===CrystalFrame,
				crystalFrame = IdentityFrame;
				labFrame = FrameInverse[OptionValue[CrystalOrientation]];
				zfsFrame = NVOrientationToFrame[OptionValue[NVOrientation]];
				zeemanFrame = FrameCompose[NVEulerAngles[0,sphericalB[[3]],sphericalB[[2]]],OptionValue[StaticFieldFrame]/.frames];,
			OptionValue[OutputFrame]===ZFSFrame,
				zfsFrame = IdentityFrame;
				crystalFrame = FrameInverse[NVOrientationToFrame[OptionValue[NVOrientation]]];
				labFrame = FrameCompose[crystalFrame,FrameInverse[OptionValue[CrystalOrientation]]];
				zeemanFrame = FrameCompose[OptionValue[StaticFieldFrame]/.frames,NVEulerAngles[0,sphericalB[[3]],sphericalB[[2]]]];,
			OptionValue[OutputFrame]===ZeemanFrame,
				zeemanFrame = IdentityFrame;
				(* The ZeemanFrame case is special because the StaticField has the option of being written in any frame. *)
				Which[
					OptionValue[StaticFieldFrame]===LabFrame,
						labFrame = FrameInverse[NVEulerAngles[0,sphericalB[[3]],sphericalB[[2]]]];
						crystalFrame = FrameCompose[OptionValue[CrystalOrientation],labFrame];
						zfsFrame = FrameCompose[NVOrientationToFrame[OptionValue[NVOrientation]],crystalFrame];,
					OptionValue[StaticFieldFrame]===CrystalFrame,
						crystalFrame = FrameInverse[NVEulerAngles[0,sphericalB[[3]],sphericalB[[2]]]];
						labFrame = FrameCompose[FrameInverse[OptionValue[CrystalOrientation]],crystalFrame];
						zfsFrame = FrameCompose[NVOrientationToFrame[OptionValue[NVOrientation]],crystalFrame];,
					OptionValue[StaticFieldFrame]===ZFSFrame,
						zfsFrame = FrameInverse[NVEulerAngles[0,sphericalB[[3]],sphericalB[[2]]]];
						crystalFrame = FrameCompose[FrameInverse[NVOrientationToFrame[OptionValue[NVOrientation]]],zfsFrame];
						labFrame = FrameCompose[FrameInverse[OptionValue[CrystalOrientation]],crystalFrame];
				];
		];

		(* Now write the magnetic field in the appropriate frame *)
		cartesianB = FrameChange[cartesianB, OptionValue[StaticFieldFrame]/.frames, IdentityFrame];

		(* Determine the spin of the NV center. *)
		nvSpin = OptionValue[NVSpin];

		(* Find the nitrogen (if it exists) and move it to the front of the list. *)
		(* All other spins stay in the same order. *)
		nucleiList = List[nuclei];
		nucleiList = Join[Select[nucleiList,NitrogenQ,1],Select[nucleiList,CarbonQ]];
		hasN = Length[nucleiList]>0 && NitrogenQ@First@nucleiList;
		numC = Length[nucleiList] - If[hasN,1,0];

		(* Determine the dimensions of each subsystem. *)
		dimNV = SpinDim[nvSpin];
		dimList = Prepend[SpinDim/@nucleiList,dimNV];
		dimN = If[hasN, SpinDim@First@nucleiList, 1];
		dimC = 2^numC;

		(* We will store all terms in termList, at first, with the head Term *)
		termList = {};

		(* Start by adding the ZFS term*)
		AppendTo[termList,Term[ZFSHamiltonian[OptionValue[ZeroFieldSplitting],nvSpin],1]];

		(* Now add the NV Zeeman term *)
		AppendTo[termList,Term[ZeemanHamiltonian[\[Gamma]e,cartesianB,nvSpin],1]];

		(* Add all of the nuclei Zeeman terms *)
		NuclearZeemanTerm[nucleus_,index_]:=Term[ZeemanHamiltonian[GyromagneticRatio[nucleus],cartesianB,SpinValue[nucleus]],index+1];
		termList = Join[termList,MapIndexed[NuclearZeemanTerm[#1,First@#2]&, nucleiList]];

		(* Calculate and store all of the hyperfine interaction terms *)
		HyperfineTerm[nucleus_,index_]:=Term[HyperfineHamiltonian[nvSpin,SpinValue[nucleus],FrameChange[Tensor[nucleus],IdentityFrame,zfsFrame]], 1, index+1];
		termList = Join[termList,MapIndexed[HyperfineTerm[#1,First@#2]&, nucleiList]];

		(* Add the quadrapolar term for the Nitrogen *)
		If[hasN,
			AppendTo[termList,Term[QuadrapolarHamiltonian[SpinValue@First@nucleiList,QuadrapoleTensor@First@nucleiList],2]];
		];

		(* Add all of the dipole-dipole interactions *)
		DipoleDipoleTerm[nucleus1_,nucleus2_,index1_,index2_]:=
			Term[DipoleDipoleHamiltonian[SpinValue[nucleus1],SpinValue[nucleus2],GyromagneticRatio[nucleus1],GyromagneticRatio[nucleus2],\[Lambda]*Location[nucleus1],\[Lambda]*Location[nucleus2]],index1+1,index2+1];
		termList = Join[termList,Flatten@Table[DipoleDipoleTerm[nucleiList[[i]],nucleiList[[j]],i,j],{i,Length[nucleiList]},{j,i+1,Length[nucleiList]}]];

		(* Define rules for turning a Term into a matrix on the full Hilbert space. *)
		(* We only have 1-local and 2-local terms, so make a rule for each one manually. *)
		Term/:ExpandTerm[Term[localHam_,ind_]] :=
			IdentityMatrix[Times@@dimList[[1;;ind-1]]]\[CircleTimes]localHam\[CircleTimes]IdentityMatrix[Times@@dimList[[ind+1;;-1]]];
		Term/:ExpandTerm[Term[localHam_,ind1_,ind2_]] :=
			IdentityInsert[localHam, 
							dimList[[ind1]], dimList[[ind2]], 
							Times@@dimList[[1;;ind1-1]], Times@@dimList[[ind1+1;;ind2-1]], Times@@dimList[[ind2+1;;-1]]
			];

		(* Decide whether or not a numerical output is desired *)
		If[
			OptionValue[Numerical]===Automatic,
				(* If any term contains something with the Head Real or Complex, set numerical to True. *)
				numerical = Not[FreeQ[termList,_Real]&&FreeQ[termList,_Complex]];,
				(* Otherwise, the user will have said what to do. *)
				numerical = OptionValue[Numerical];
		];
		(* Decide whether or not to multiply everything by 2\[Pi]; the Automatic logic is the same as numerical *)
		If[OptionValue[AngularUnits]===Automatic, angularUnits = numerical;, angularUnits = OptionValue[AngularUnits];];
		preFactor = If[angularUnits, preFactor=2*\[Pi], 1];

		(* Not that it _really _ matters, but do the search and replace for constants before expanding to the
		   full Hilbert space because it should be a bit faster, at least, asymtotically.  *)
		If[numerical,
			termList = InsertConstants[termList];
		];

		(* Finally, expand all the terms and sum them up. *)
		(* We are carfeful to do this in a way that doesn't use more memory than needed. *)
		Fold[#1+preFactor*ExpandTerm[#2]&, 0, termList]
	]


End[];


(* ::Section:: *)
(*Approximation Tools*)


(* ::Subsection:: *)
(*Usage Declarations*)


NVAverageHamiltonian::usage="NVAverageHamiltonian[order,\[Omega]rot,...] returns the order'th order average NV Hamiltonian in the \[Omega]rot*Sz^2 frame. Always specify \[Omega]rot in terms of frequency; it is multipled by 2\[Pi] to get angular frequencey if the Hamiltonian uses angular frequency units. Here, '...' represents any arguments you would pass into the function NVHamiltonian. If order is negative, the effective Hamiltonian is returned.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Knowing what happens when we conjugate by the rotating frame unitary, we can precompute the average Hamiltonians a bit in terms of the floquet coefficients. In partucular, we can precompute all of the intergrals, and are just left with sums of nested comutators. These are the formulas used below.*)


NVAverageHamiltonian[order_,\[Omega]rot_,nuclei___,opt:OptionsPattern[NVHamiltonian]]:=
	Module[
		{H,H0,U,ndim,sgn,t,\[Omega],numerical,angularUnits,Heff,Hm,Hp,Hout,z,rotation},

		(* First calculate the lab frame hamiltonian *)
		H = NVHamiltonian[nuclei, opt];

		(* Decide whether or not a numerical output is desired. Same logic as in NVHamiltonian *)
		If[OptionValue[Numerical]===Automatic,
			numerical = Not[FreeQ[H,_Real]&&FreeQ[H,_Complex]];,
			numerical = OptionValue[Numerical];
		];
		If[OptionValue[AngularUnits]===Automatic, angularUnits = numerical;, angularUnits = OptionValue[AngularUnits];];

		(* Get some numbers we will need. *)
		ndim = Times@@(SpinDim/@{nuclei});
		\[Omega]=If[numerical,InsertConstants[\[Omega]rot],\[Omega]rot];
		If[angularUnits,\[Omega]=2*\[Pi]*\[Omega]];

    	(* Finally, do the actual calculation *)
		If[order < 0,
			H0 = ZFSHamiltonian[\[Omega], OptionValue[NVSpin]]\[CircleTimes]IdentityMatrix[ndim];
			U[sgn_] = MatrixExp[sgn I t H0];
			Heff[s_]=U[1].(H-H0).U[-1]/.t->s;Heff,

			z=ConstantArray[0,{ndim,ndim}];
			rotation=\[Omega]*IdentityMatrix[ndim];
			(* Compute the three Floquet coefficients *)
			H0 = ArrayFlatten[{{H[[1;;ndim,1;;ndim]]-rotation,z,H[[1;;ndim,2*ndim+1;;-1]]},{z,H[[ndim+1;;2*ndim,ndim+1;;2*ndim]],z},{H[[2*ndim+1;;-1,1;;ndim]],z,H[[2*ndim+1;;-1,2*ndim+1;;-1]]-rotation}}];
			Hm = ArrayFlatten[{{z,z,z},{H[[ndim+1;;2*ndim,1;;ndim]],z,H[[ndim+1;;2*ndim,2*ndim+1;;-1]]},{z,z,z}}];
			Hp = ArrayFlatten[{{z,H[[1;;ndim,ndim+1;;2*ndim]],z},{z,z,z},{z,H[[2*ndim+1;;-1,ndim+1;;2*ndim]],z}}];
			Hout=H0;

			If[order>=1,
				Hout=Hout+(Com[H0,Hp]-Com[H0,Hm]-Com[Hm,Hp])/\[Omega];
			];
			If[order>=2,
				Hout=Hout+(
					 (Com[Hm,Com[H0,Hp]]+Com[Hp,Com[H0,Hm]])
					-2(Com[H0,Com[H0,Hp]]+Com[H0,Com[H0,Hm]])
					+(Com[Hp,Com[Hp,H0]]+Com[Hm,Com[Hm,H0]]-Com[Hm,Com[Hp,H0]]-Com[Hp,Com[Hm,H0]])
					+2(Com[Hp,Com[Hm,Hp]]+Com[Hm,Com[Hp,Hm]]))/(2*\[Omega]^2);
			];
			If[order>=3,Print["Warning: Third order AH not implemented. Truncating to second order instead."];];
			Hout
		]
	]


End[];


(* ::Section::Closed:: *)
(*Epilogue*)


(* ::Text:: *)
(*This section exists so that EndPackage[] isn't just placed in the last section. Not doing this can result in very confusing behaviour when you add another section but forget to move EndPackage[].*)


EndPackage[];
