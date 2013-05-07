(* ::Package:: *)

BeginPackage["NVHamiltonian`"];


(* ::Section:: *)
(*Error Messages*)


(* ::Subsection:: *)
(*Matrices, Bases, and Linear Algebra*)


Spin::badindex = "An invalid spin index was entered.";


IdentityInsert::baddimensions = "The size of the input matrix does not match the product of the specified dimensions.";


(* ::Subsection:: *)
(*Frames*)


Frame::notorthonormal = "The first three arguments of Frame must form an orthonormal basis for R^3.";


(* ::Section:: *)
(*Predicates*)


(* ::Subsection:: *)
(*Usage Declarations*)


Carbon::usage = ""; Nitrogen::usage = "";
CarbonQ::usage = "CarbonQ[c] returns True iff the Head of c is Carbon";
NitrogenQ::usage = "NitrogenQ[c] returns True iff the Head of c is Carbon";
NucleusQ::usage = "CarbonQ[c] returns True iff one of CarbonQ or NitrogenQ holds.";


Vector3Q::usage = "Vector3Q[v] returns True iff v is a vector of length 3.";
Matrix3Q::usage = "Matrix3Q[m] returns True iff m is a matrix of height 3.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


CarbonQ[c_]:=Head[c]===Carbon
NitrogenQ[c_]:=Head[c]===Nitrogen
NucleusQ[c_]:=CarbonQ[c]||NitrogenQ[c]


Vector3Q[v_]:=VectorQ[v]&&Length[v]==3
Matrix3Q[m_]:=MatrixQ[m]&&Length[m]==3


End[];


(* ::Section::Closed:: *)
(*Matrices, Bases, and Linear Algebra*)


(* ::Subsection:: *)
(*Usage Declarations*)


X::usage = "The 2x2 Pauli X operator.";
Y::usage = "The 2x2 Pauli Y operator.";
Z::usage = "The 2x2 Pauli Z operator.";


Sx::usage = "The 3x3 Spin-1 X operator.";
Sy::usage = "The 3x3 Spin-1 Y operator.";
Sz::usage = "The 3x3 Spin-1 Z operator.";


Spin::usage = "Spin[s] returns the spin operators for spin s particles. Spin[Carbon[...]] and Spin[Nitrogen[...]] return the spin operators of the respected nucleus.";
SpinDim::usage = "SpinDim[s] returns the dimension of Hilbert space corresponding to spin number s. For example, s=1/2 will return 2, and s=1 will return 3.";


IdentityInsert::usage = "IdentityInsert[C,dimA,dimB,n1,n2,n3] assumes C is a square matrix acting on a bipartite system with dimensions dimA and dimB and proceeds to add identity operations before, in-between, and after the bipartite system, with dimensions n1, n2, and n3, respectively.";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*Spin operators*)


Begin["`Private`"];


X={{0,1},{1,0}};
Y={{0,-I},{I,0}};
Z={{1,0},{0,-1}};


Sx={{0,1,0},{1,0,1},{0,1,0}}/Sqrt[2];
Sy={{0,-I,0},{I,0,-I},{0,I,0}}/Sqrt[2];
Sz={{1,0,0},{0,0,0},{0,0,-1}};


Spin[1/2]={X,Y,Z}/2;
Spin[1]={Sx,Sy,Sz};


Spin[s_,i_]:=
	If[i>3||i<1,
		Message[Spin::badindex];Abort;,		
		Spin[s][[i]]
	]


SpinDim[s_?NumericQ]:=2*s+1


End[];


(* ::Subsubsection:: *)
(*Kronecker product*)


Begin["`Private`"];


CircleTimes[A_]:=A
CircleTimes[A__]:=KroneckerProduct[A]


End[];


(* ::Subsubsection:: *)
(*Tensor Manipulation*)


Begin["`Private`"];


(* ::Text:: *)
(*First deal with the most difficult case where we want to insert an identity operation into the middle.*)


IdentityInsert[C_,dimA_,dimB_,n_]:=
	If[Length[C]!=dimA*dimB,
		Message[IdentityInsert::baddimensions];Abort;,
		ArrayFlatten[Map[IdentityMatrix[n]\[CircleTimes]#&,Partition[C,{dimA,dimB}],{2}]]
	]


(* ::Text:: *)
(*Now expand the definition to the left and right sides too.*)


IdentityInsert[C_,dimA_,dimB_,n1_,n2_,n3_]:=IdentityMatrix[n1]\[CircleTimes]IdentityInsert[C,dimA,dimB,n2]\[CircleTimes]IdentityMatrix[n3]


End[];


(* ::Section:: *)
(*Frames*)


(* ::Subsection:: *)
(*Usage Declarations*)


SphericalToCartesian::usage = "SphericalToCartesian[{r,\[Theta],\[Phi]}] transforms the spherical vector {r,\[Theta],\[Phi]} into the cartesian vector {r Cos[\[Theta]]Sin[\[Phi]], r Sin[\[Theta]]Sin[\[Phi]], r Cos[\[Phi]]}. \[Phi] is the polar angle (down from the z axis) and \[Theta] is the azimual angle (from the x axis towards the y axis). We choose not to use the builtin function CoordinateTransform because it deals with 0 angles in a dumb way.";
CylindricalToCartesian::usage = "CylindricalToCartesian[{\[Rho],\[Theta],z}] transforms the cylindrical vector {\[Rho],\[Theta],z} into the cartesian vector {\[Rho] Cos[\[Theta]], \[Rho] Sin[\[Theta]], z}. \[Theta] is the azimual angle (from the x axis towards the y axis). We choose not to use the builtin function CoordinateTransform because it deals with 0 angles in a dumb way.";
CartesianToSpherical::usage = "CartesianToSpherical[{x,y,z}] transforms the cartesian vector {x,y,z} into the spherical vector {r,\[Theta],\[Phi]}. \[Phi] is the polar angle (down from the z axis) and \[Theta] is the azimual angle (from the x axis towards the y axis). If both x and y are 0, then the ambiguous azimuthal angle \[Theta] is taken to be 0.";
CylindricalToSpherical::usage = "CylindricalToSpherical[{\[Rho],\[Theta],z}] transforms the cylindrical vector {\[Rho],\[Theta],z} into the spherical vector {r,\[Theta],\[Phi]}. \[Phi] is the polar angle (down from the z axis) and \[Theta] is the azimual angle (from the x axis towards the y axis). If both \[Rho] and z are 0, then the ambiguous polar angle \[Phi] is taken to be 0.";
CartesianToCylindrical::usage = "CartesianToCylindrical[{x,y,z}] transforms the cartesian vector {x,y,z} into the cylindrical vector {\[Rho],\[Theta],z}. \[Theta] is the azimual angle (from the x axis towards the y axis). If both x and y are 0, then the ambiguous azimuthal angle \[Theta] is taken to be 0.";
SphericalToCylindrical::usage = "SphericalToCylindrical[{r,\[Theta],\[Phi]}] transforms the spherical vector {r,\[Theta],\[Phi]} into the cylindrical vector {\[Rho],\[Theta],z}.  is the polar angle (down from the z axis) and \[Theta] is the azimual angle (from the x axis towards the y axis).";
ChangeCoordinates::usage = "ChangeCoordinates[v, fromCoords, toCoords] changes the coordinates of the length-3 vector. For example, ChangeCoordinates[{1,1,0},Cartesian,Spherical].";


Frame::usage = "";


Coordinates::usage = "";
Cartesian::usage = "";
Spherical::usage = "";
Cylindrical::usage = "";


FrameMatrix::usage = "FrameMatrix[frame] returns the orthogonal 3x3 matrix, M, corresponding to the given frame. This is the matrix that transforms from the canonical basis to the basis of the frame, all in Cartesian coordinates (Frames entered in non-Cartesian cooridinates will be automatically converted). So, for example, M.{1,0,0}=x, where frame=Frame[x,y,z,Cartesian].";
FrameChangeMatrix::usage = "";


FrameInverse::usage = "FrameInverse[frame]";
FrameChange::usage = "FrameChange[M,fromFrame,toFrame] or FRame[v,fromFrame,ToFrame]";
FrameCompose::usage = "FrameCompose[framen,...,frame2,frame1] returns the resulting frame when all of the input frames are composed. That is, we know frame1 is written in the coordinates of the canonical basis, and if frame2 is written in the coordinates of frame1, and frame3 is written in the coordinates of frame2, etc, then the resulting frame is the composition of all frames written in the canonical coordinates.";


IdentityFrame::usage = "";
BondFrame::usage = "";
EulerAngles::usage = "EulerAngles[\[Theta]z1,\[Theta]y,\[Theta]z2] returns a Frame corresponding to rotating the IdenityFrame by the extrinsic ZYZ Euler angles \[Theta]z1,\[Theta]y,\[Theta]z2.";


PlotFrame::usage = "PlotFrame[frame1,frame2,...] plots each Frame given as an argement on the same figure.";


(* ::Subsection:: *)
(*Implementations*)


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


(* ::Text:: *)
(*Returns the orthogonal matrix corresponding to the input frame.*)


FrameMatrix[f_Frame]:=((List@@Cartesian[f])[[1;;3]])\[Transpose]


(* ::Text:: *)
(*You can think of converting between frames as first converting to the canonical basis, and then from this converting to the desired frame.*)


FrameChangeMatrix[fromFrame_,toFrame_]:=FrameMatrix[toFrame].FrameMatrix[fromFrame]\[Transpose]


(* ::Text:: *)
(*We have the frame change matrix, so now it's just a matter of actually impementing it.*)


FrameChange[M_?Vector3Q,fromFrame_,toFrame_]:=FrameChangeMatrix[fromFrame,toFrame].M
FrameChange[M_?Matrix3Q,fromFrame_,toFrame_]:=With[{F=FrameChangeMatrix[fromFrame,toFrame]},F.M.F\[Transpose]]


(* ::Text:: *)
(*Since the matrices are orthogonal, the inverse is given by the transpose.*)


FrameInverse[f_Frame]:=Frame[FrameMatrix[f]\[Transpose]]


(* ::Text:: *)
(*Composition is of course just matrix multiplication in Cartesian coordinates.*)


FrameCompose[f_Frame]:=f
FrameCompose[fa_Frame,fb_Frame,rest___]:=FrameCompose[Frame[Simplify[FrameMatrix[fa].FrameMatrix[fb]]],rest]


IdentityFrame=Frame[{1,0,0},{0,1,0},{0,0,1},Cartesian];


EulerAngles[\[Theta]z1_,\[Theta]y_,\[Theta]z2_]=Frame[RotationMatrix[\[Theta]z2, {0,0,1}].RotationMatrix[\[Theta]y, {0,1,0}].RotationMatrix[\[Theta]z1, {0,0,1}]];


(* ::Text:: *)
(*The bond frame is most easily describable in spherical coordinates, so convert first. Note that all "bond frame means" is some frame in which the z vector is parallel to the input vector of BondFrame; the x-y vectors are chosen sort of arbitrarily, but it doesn't matter because the tensor should be cylindrically symmetric.*)


BondFrame[{x_,y_,z_},Cartesian]:=BondFrame[CartesianToSpherical[{x,y,z}],Spherical];
BondFrame[{\[Rho]_,\[Theta]_,z_},Cylindrical]:=BondFrame[CylindricalToSpherical[{\[Rho],\[Theta],z}],Spherical];
BondFrame[{r_,\[Theta]_,\[Phi]_},Spherical]:=EulerAngles[0,\[Phi],\[Theta]]


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


(* ::Section:: *)
(*NV and Lattice Geometry*)


(* ::Subsection:: *)
(*Usage Declarations*)


NVOrientationToFrame::usage = "NVOrientationToFrame[n] returns a Frame (with respect to the crystal frame, of course) corresponding to the n'th NV orientation. n should be one of the values 1,2,3,4,5,6,7,8. Here, 1 is along the positive z direction, 2 is on the x-z plane, rotated an angle ArcCos[-1/3] from the z axis, and orientations 3 and 4 are right-handed Z rotations of orientation 2 by angles 2\[Pi]/3 and 4\[Pi]/3 respectively. The orientations 5,6,7,8 are respectively anti-parallel to 1,2,3,4.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


NVOrientationToFrame[1]=IdentityFrame;
NVOrientationToFrame[2]=EulerAngles[0,ArcCos[-1/3],0];
NVOrientationToFrame[3]=EulerAngles[0,ArcCos[-1/3],2\[Pi]/3];
NVOrientationToFrame[4]=EulerAngles[0,ArcCos[-1/3],4\[Pi]/3];
NVOrientationToFrame[5]=EulerAngles[0,\[Pi],0];
NVOrientationToFrame[6]=EulerAngles[0,ArcCos[-1/3]-\[Pi],0];
NVOrientationToFrame[7]=EulerAngles[0,ArcCos[-1/3]-\[Pi],2\[Pi]/3];
NVOrientationToFrame[8]=EulerAngles[0,ArcCos[-1/3]-\[Pi],4\[Pi]/3];


End[]


(* ::Section:: *)
(*Options*)


(* ::Subsection:: *)
(*Usage Declarations*)


NVHamiltonian::usage = "";
NVOrientation::usage = "";
ZeroFieldSplitting::usage = "";
StaticField::usage = "";
StaticFieldFrame::usage = "";
StaticFieldCoordinates::usage = "";
NVSpin::usage = "";
OutputFrame::usage = "";
CrystalOrientation::usage = "";


LabFrame::usage = "";
ZFSFrame::usage = "";
CrystalFrame::usage =  "";
ZeemanFrame::usage = "";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


Options[NVHamiltonian]=
	{
		ZeroFieldSplitting -> \[CapitalDelta],
		NVOrientation -> 1,
		StaticField -> {0,0,0},
		StaticFieldFrame -> LabFrame,
		StaticFieldCoordinates -> Cartesian,
		NVSpin -> 1,
		OutputFrame -> ZFSFrame,		
		CrystalOrientation -> IdentityFrame
	};


(* ::Text:: *)
(*This loops through all of the options in NVHamiltonian and protects them.*)


Protect[Evaluate[Sequence@@(Options[NVHamiltonian][[All,1]])]];


End[];


(* ::Section:: *)
(*Hyperfine Operations*)


(* ::Subsection:: *)
(*Usage Declarations*)


RotateTensor::usage = "Rotation[A, {\[Theta]azimuth, \[Theta]polar}] performs the rotation RotZ[-\[Theta]azimuth].RotY[-\[Theta]polar].A.RotY[\[Theta]polar].RotZ[\[Theta]azimuth]. If these angles are interpreted as the spherical coordinates of the nucleus in the ZFS frame, then this operation will take you _from_ the bond frame _to_ the ZFS frame.";


HyperfineHamiltonian::usage = "HyperfineHamiltonian[spin1,spin2,A]";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


RotateTensor[A_,{\[Theta]azimuth_,\[Theta]polar_}]:=
	With[{O=RotationMatrix[-\[Theta]azimuth, {0,0,1}].RotationMatrix[-\[Theta]polar, {0,1,0}]},
		O.A.O\[Transpose]
	];


(* ::Text:: *)
(*Start by defining the hyperfine interaction on a bipartite system.*)


HyperfineHamiltonian[spin1_,spin2_,A_]:=
	Total[Table[A[[i,j]]*Spin[spin1,i]\[CircleTimes]Spin[spin2,j],{i,3},{j,3}],2]


(* ::Text:: *)
(*Now expand the definition to include identity operations on systems before, in-between, and after.*)


HyperfineHamiltonian[spin1_,spin2_,A_,n1_,n2_,n3_]:=IdentityInsert[HyperfineHamiltonian[spin1,spin2,A],SpinDim[spin1],SpinDim[spin2],n1,n2,n3]


End[];


(* ::Section:: *)
(*Nuclei*)


(* ::Subsection:: *)
(*Usage Declarations*)


Carbon::usage = "";
Nitrogen::usage = "";


Tensor::usage = "";
QuadrapoleTensor::usage = "";
Isotope::usage = "";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*Carbon*)


Begin["`Private`"];


Carbon/:Tensor[Carbon[{Apar_,Aperp_}]]:=DiagonalMatrix[{Aperp,Aperp,Apar}]
Carbon/:Tensor[Carbon[{Apar_,Aperp_},{\[Theta]azimuth_,\[Theta]polar_}]]:=RotateTensor[Tensor@Carbon[{Apar,Aperp}],{\[Theta]azimuth,\[Theta]polar}]
Carbon/:Tensor[Carbon[A_?MatrixQ]]:=A
Carbon/:Tensor[Carbon[{Azx_,Azy_,Azz_}]]:={{0,0,0},{0,0,0},{Azx,Azy,Azz}}


Carbon/:Isotope[Carbon[___]]:=13


Carbon/:Spin[Carbon[___]]:=1/2


End[];


(* ::Subsubsection:: *)
(*Nitrogen*)


Begin["`Private`"];


Nitrogen/:Tensor[Nitrogen[_,{Apar_,Aperp_},___]]:=DiagonalMatrix[{Aperp,Aperp,Apar}]
Nitrogen/:Tensor[Nitrogen[_,Azz_,___]]:=DiagonalMatrix[{0,0,Azz}]
Nitrogen/:Tensor[Nitrogen[_,A_?MatrixQ,___]]:=A
Nitrogen/:Tensor[Nitrogen[isotope_]]:=Tensor@Nitrogen[isotope,{ANpar,ANperp}]


Nitrogen/:QuadrapoleTensor[Nitrogen[_,_,Q_]]:=DiagonalMatrix[{0,0,Q}]
Nitrogen/:QuadrapoleTensor[Nitrogen[_,_,Q_?MatrixQ]]:=Q


Nitrogen/:Isotope[Nitrogen[isotope_,___]]:=isotope


Nitrogen/:Spin[Nitrogen[x___]]:=Spin[If[Isotope[Nitrogen[x]]===13,1/2,1]]


End[];


(* ::Section:: *)
(*Hamiltonians*)


(* ::Subsection:: *)
(*Usage Declarations*)


ZFSHamiltonian::usage = "";


NVHamiltonian::usage = "";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*ZFS Hamlitonian*)


Begin["`Private`"];


ZFSHamiltonian[OptionsPattern[NVHamiltonian]]:=OptionValue[ZeroFieldSplitting] Sz


End[];


(* ::Subsubsection:: *)
(*Zeeman Hamiltonian*)


Begin["`Private`"];


ZeemanHamiltonian[spin_,{\[Omega]x_,\[Omega]y_,\[Omega]z_}] := \[Omega]x Spin[spin,1] + \[Omega]x Spin[spin,2] + \[Omega]x Spin[spin,3]


End[];


(* ::Subsubsection:: *)
(*Combined Hamiltonian*)


Begin["`Private`"];


NVHamiltonian[nuclei___,opt:OptionsPattern[NVHamiltonian]]:=
	Module[
		{
			cartesianB,sphericalB,
			dimNV,dimN,hasN,dimC,numC,dimList,
			nucleiList,
			nvSpin,
			termList,Term,ExpandTerm,HyperfineTerm,
			labFrame,crystalFrame,zfsFrame,zeemanFrame
		},

		(* Convert the B field to various coordinate systems *)
		cartesianB = ChangeCoordinates[OptionValue[StaticField],OptionValue[StaticFieldCoordinates],Cartesian];
		sphericalB = ChangeCoordinates[OptionValue[StaticField],OptionValue[StaticFieldCoordinates],Spherical];

		(* Deal with all of the frames first (except nuclei frames) *)
		(* The idea is simple: we have a natural ordering of frames: LabFrame>CrystalFrame>ZFSFrame>Zeeman *)
		Which[
			OptionValue[OutputFrame]===LabFrame,
				labFrame = IdentityFrame;
				crystalFrame = OptionValue[CrystalOrientation];
				zfsFrame = FrameCompose[NVOrientationToFrame[OptionValue[NVOrientation]],crystalFrame];
				zeemanFrame = FrameCompose[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]],zfsFrame];
			OptionValue[OutputFrame]===CrystalFrame,
				crystalFrame = IdentityFrame;
				labFrame = FrameInverse[OptionValue[CrystalOrientation]];
				zfsFrame = NVOrientationToFrame[OptionValue[NVOrientation]];
				zeemanFrame = FrameCompose[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]],zfsFrame];
			OptionValue[OutputFrame]===ZFSFrame,
				zfsFrame = IdentityFrame;
				crystalFrame = FrameInverse[NVOrientationToFrame[OptionValue[NVOrientation]]];
				labFrame = FrameCompose[OptionValue[CrystalOrientation],crystalFrame];
				zeemanFrame = EulerAngles[0,sphericalB[[3]],sphericalB[[2]]];
			OptionValue[OutputFrame]===ZeemanFrame,
				zeemanFrame = IdentityFrame;
				zfsFrame = FrameInverse[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]]];
				crystalFrame = FrameCompose[FrameInverse[NVOrientationToFrame[OptionValue[NVOrientation]]],zfsFrame];
				labFrame = FrameCompose[OptionValue[CrystalOrientation],crystalFrame];
		];

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

		(* Calculate and store all of the hyperfine interaction terms *)
		HyperfineTerm[nucleus_,index_]:=Term[HyperfineHamiltonian[nvSpin,Spin[nucleus],Tensor[nucleus]], 1, index];
		termList = Join[termList,MapIndexed[HyperfineTerm[#1,First@#2]&, nucleiList]];

		(* Define rules for turning a Term into a matrix on the full Hilbert space. *)
		(* We only have 1-local and 2-local terms, so make a rule for each one manually. *)
		Term/:ExpandTerm[Term[localHam_,ind_]] :=
			IdentityMatrix[Times@@dimList[[1;;ind-1]]]\[CircleTimes]localHam\[CircleTimes]IdentityMatrix[Times@@dimList[[ind+1;;-1]]];
		Term/:ExpandTerm[Term[localHam_,ind1_,ind2_]] :=
			IdentityInsert[localHam, 
							dimList[[ind1]], dimList[[ind2]], 
							Times@@dimList[[1;;ind1-1]], Times@@dimList[[ind1+1;;ind2-1]], Times@@dimList[[ind2+1;;-1]]
			];

		(* Finally, expand all the terms and sum them up. *)
		(* We are carfeful to do this in a way that doesn't use more memory than needed. *)
		Fold[#1+ExpandTerm[#2]&, 0, termList]
	]


End[];


(* ::Section::Closed:: *)
(*Epilogue*)


EndPackage[];
