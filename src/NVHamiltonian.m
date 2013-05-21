(* ::Package:: *)

BeginPackage["NVHamiltonian`"];


(* ::Text:: *)
(*Todo:*)
(*-Verify BondFrame is correct*)
(*-Change the static field frame if necessary*)


(* ::Section::Closed:: *)
(*Error Messages*)


(* ::Subsection::Closed:: *)
(*Matrices, Bases, and Linear Algebra*)


Spin::badindex = "An invalid spin index was entered.";


IdentityInsert::baddimensions = "The size of the input matrix does not match the product of the specified dimensions.";


(* ::Subsection::Closed:: *)
(*Frames and Vectors*)


Frame::notorthonormal = "The first three arguments of Frame must form an orthonormal basis for R^3.";


(* ::Subsection::Closed:: *)
(*Options*)


Numerical::badinput = "The option Numerical must be set to one of True, False, or Automatic.";


AngularUnits::badinput = "The option AngularUnits must be set to one of True, False, or Automatic.";


ZeroFieldSplitting::badinput = "The ZeroFieldSplitting must be a single number/symbol, or a list of length three (Each one will be multiplied by Sx.Sx, Sy.Sy, and Sz.Sz respectively.).";


OutputFrame::badframe = "The output frame must be one of LabFrame, CrystalFrame, ZFSFrame, or ZeemanFrame.";


StaticField::badframe = "The static field's frame must be one of LabFrame, CrystalFrame, or ZFSFrame.";
StaticField::notvector = "The static field must be input as a Vector, e.g., Vector[{0,0,10}, Cartesian].";


(* ::Subsection::Closed:: *)
(*Hamiltonians*)


DipoleDipoleHamiltonian::equallocations = "The dipole-dipole Hamilotian was requested for two spins at the same physical locatian; this would result in division by 0. Instead, the 0 Hamiltonian has been return as the dipole-dipole Hamiltonian.";


(* ::Section::Closed:: *)
(*Predicates*)


(* ::Subsection::Closed:: *)
(*Usage Declarations*)


CarbonQ::usage = "CarbonQ[c] returns True iff the Head of c is Carbon";
NitrogenQ::usage = "NitrogenQ[c] returns True iff the Head of c is Carbon";
NucleusQ::usage = "CarbonQ[c] returns True iff one of CarbonQ or NitrogenQ holds.";


Vector3Q::usage = "Vector3Q[v] returns True iff v is a vector of length 3.";
Matrix3Q::usage = "Matrix3Q[m] returns True iff m is a matrix of height 3.";


ValidReferenceFrameQ::usage = "ValidReferenceFrameQ[input] returns True iff input is one of LabFrame, CrystalFrame, ZFSFrame, or ZeemanFrame."


(* ::Text:: *)
(*We need to initialize the following things now so that they fall in the right context.*)


LabFrame::usage = "";ZFSFrame::usage = "";CrystalFrame::usage =  "";ZeemanFrame::usage = "";
Carbon::usage = ""; Nitrogen::usage = "";


(* ::Subsection::Closed:: *)
(*Implementations*)


Begin["`Private`"];


CarbonQ[c_]:=Head[c]===Carbon
NitrogenQ[c_]:=Head[c]===Nitrogen
NucleusQ[c_]:=CarbonQ[c]||NitrogenQ[c]


Vector3Q[v_]:=VectorQ[v]&&Length[v]==3
Matrix3Q[m_]:=MatrixQ[m]&&Length[m]==3


ValidReferenceFrameQ[input_]:=MemberQ[{LabFrame,ZFSFrame,CrystalFrame,ZeemanFrame},input]


End[];


(* ::Section:: *)
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


Com::usage = "Com[A,B]=A.B-B.A";
NestCom::usage="NestCom[B,A,n] returns the nested commutator if B and A:
NestCom[B,A,0]=A,
NestCom[B,A,1]=Com[B,A],
...
NestCom[B,A,n]=Com[B,NestCom[B,A,n-1]].";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection::Closed:: *)
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
(*Algebra*)


Begin["`Private`"];


CircleTimes[A_]:=A
CircleTimes[A__]:=KroneckerProduct[A]


Com[A_,B_]:=A.B-B.A


NestCom[B_,A_,0]:=A;
NestCom[B_,A_,n_Integer]:=NestCom[B,A,n]=Com[B,NestCom[B,A,n-1]]


End[];


(* ::Subsubsection::Closed:: *)
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


(* ::Section::Closed:: *)
(*Physical Quantities*)


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
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


InsertConstants[expr_]:=expr/.$constants


End[];


(* ::Section::Closed:: *)
(*Frames and Vectors*)


(* ::Subsection::Closed:: *)
(*Usage Declarations*)


Coordinates::usage = "";
Cartesian::usage = "";
Spherical::usage = "";
Cylindrical::usage = "";


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
EulerAngles::usage = "EulerAngles[\[Theta]z1,\[Theta]y,\[Theta]z2] returns a Frame corresponding to rotating the IdenityFrame by the extrinsic ZYZ Euler angles \[Theta]z1,\[Theta]y,\[Theta]z2.";


PlotFrame::usage = "PlotFrame[frame1,frame2,...] plots each Frame given as an argement on the same figure.";


(* ::Subsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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


End[];


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*Changing, Inverting and Composing Frames*)


Begin["`Private`"];


(* ::Text:: *)
(*You can think of converting between frames as first converting to the canonical basis, and then from this converting to the desired frame.*)


FrameChangeMatrix[fromFrame_,toFrame_]:=FrameMatrix[toFrame].FrameMatrix[fromFrame]\[Transpose]


(* ::Text:: *)
(*We have the frame change matrix, so now it's just a matter of actually impementing it.*)


FrameChange[v_?Vector3Q,fromFrame_,toFrame_]:=FrameChangeMatrix[fromFrame,toFrame].v
FrameChange[M_?Matrix3Q,fromFrame_,toFrame_]:=With[{F=FrameChangeMatrix[fromFrame,toFrame]},F.M.F\[Transpose]]
FrameChange[v_Vector,fromFrame_,toFrame_]:=ChangeCoordinates[FrameChangeMatrix[fromFrame,toFrame].ChangeCoordinates[v,Cartesian],Cartesian,Coordinates@v]


(* ::Text:: *)
(*Since the matrices are orthogonal, the inverse is given by the transpose.*)


FrameInverse[f_Frame]:=Frame[FrameMatrix[f]\[Transpose]]


(* ::Text:: *)
(*Composition is of course just matrix multiplication in Cartesian coordinates.*)


FrameCompose[f_Frame]:=f
FrameCompose[fa_Frame,fb_Frame,rest___]:=FrameCompose[Frame[Simplify[FrameMatrix[fa].FrameMatrix[fb]]],rest]


End[];


(* ::Subsubsection::Closed:: *)
(*Special Frames*)


Begin["`Private`"];


(* ::Text:: *)
(*The canonical basis.*)


IdentityFrame=Frame[{1,0,0},{0,1,0},{0,0,1},Cartesian];


(* ::Text:: *)
(*We use the extrinsic ZYZ convention.*)


EulerAngles[\[Theta]z1_,\[Theta]y_,\[Theta]z2_]=Frame[RotationMatrix[\[Theta]z2, {0,0,1}].RotationMatrix[\[Theta]y, {0,1,0}].RotationMatrix[\[Theta]z1, {0,0,1}]];


(* ::Text:: *)
(*The bond frame is most easily describable in spherical coordinates, so convert first. Note that all "bond frame means" is some frame in which the z vector is parallel to the input vector of BondFrame; the x-y vectors are chosen sort of arbitrarily, but it doesn't matter because the tensor should be cylindrically symmetric.*)


BondFrame[v_Vector]:=With[{s=Value@ChangeCoordinates[v,Spherical]},EulerAngles[0,s[[3]],s[[2]]]]


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


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection::Closed:: *)
(*NV Orientations*)


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


(* ::Subsubsection::Closed:: *)
(*Lattice Positions*)


Begin["`Private`"];


(* ::Text:: *)
(*The lattice directions.*)


E0=Vector[{0,0,1},Cartesian];
E1=Cartesian@Vector[{1,0,ArcCos[-1/3]},Spherical];
E2=Cartesian@Vector[{1,2\[Pi]/3,ArcCos[-1/3]},Spherical];
E3=Cartesian@Vector[{1,4\[Pi]/3,ArcCos[-1/3]},Spherical];


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

		If[ListQ[OptionValue[ZeroFieldSplitting]]&&Not[Vector3Q[OptionValue[ZeroFieldSplitting]]],Message[ZeroFieldSplitting::badinput];abort=True];

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


Tensor::usage = "";
Location::usage = "";
QuadrapoleTensor::usage = "";
Isotope::usage = "";
GyromagneticRatio::usage = "";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection:: *)
(*Carbon*)


Begin["`Private`"];


Carbon/:Tensor[Carbon[{Apar_,Aperp_},v_Vector]]:=FrameChange[DiagonalMatrix[{Aperp,Aperp,Apar}],IdentityFrame,BondFrame[v]]
Carbon/:Tensor[Carbon[A_?Matrix3Q]]:=A


Carbon/:Location[Carbon[_,v_Vector]]:=v


Carbon/:Isotope[Carbon[___]]:=13


Carbon/:Spin[Carbon[___]]:=1/2


Carbon/:SpinDim[Carbon[___]]:=2


Carbon/:GyromagneticRatio[Carbon[___]]:=\[Gamma]c


End[];


(* ::Subsubsection::Closed:: *)
(*Nitrogen*)


Begin["`Private`"];


Nitrogen/:Tensor[Nitrogen[_,{Apar_,Aperp_},___]]:=DiagonalMatrix[{Aperp,Aperp,Apar}]


Nitrogen/:Location[Nitrogen[___]]:=E0


Nitrogen/:QuadrapoleTensor[Nitrogen[_,_,Q_?Matrix3Q]]:=Q
Nitrogen/:QuadrapoleTensor[Nitrogen[_,_,Q_?Vector3Q]]:=DiagonalMatrix[Q]
Nitrogen/:QuadrapoleTensor[Nitrogen[_,_,Q_]]:=DiagonalMatrix[{0,0,Q}]


Nitrogen/:Isotope[Nitrogen[isotope_,___]]:=isotope


Nitrogen/:Spin[Nitrogen[x___]]:=If[Isotope[Nitrogen[x]]===13,1/2,1]


Nitrogen/:SpinDim[n:Nitrogen[___]]:=SpinDim[Spin[n]]


Nitrogen/:GyromagneticRatio[Nitrogen[x___]]:=If[Isotope[Nitrogen[x]]===13,\[Gamma]n13,\[Gamma]n14]


End[];


(* ::Section:: *)
(*Hamiltonians*)


(* ::Subsection::Closed:: *)
(*Usage Declarations*)


ZFSHamiltonian::usage = "";


ZeemanHamiltonian::usage = "";


HyperfineHamiltonian::usage = "HyperfineHamiltonian[spin1,spin2,A]";


DipoleDipoleHamiltonian::usage = "";


QuadrapolarHamiltonian::usage = "";


NVHamiltonian::usage = "";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection::Closed:: *)
(*ZFS Hamlitonian*)


Begin["`Private`"];


ZFSHamiltonian[{Dx_,Dy_,Dz_},nvSpin_]:= Dx Spin[nvSpin,1].Spin[nvSpin,1] + Dy Spin[nvSpin,2].Spin[nvSpin,2] + Dz Spin[nvSpin,3].Spin[nvSpin,3];
ZFSHamiltonian[D_,nvSpin_]:= D Spin[nvSpin,3].Spin[nvSpin,3];


End[];


(* ::Subsubsection::Closed:: *)
(*Zeeman Hamiltonian*)


Begin["`Private`"];


ZeemanHamiltonian[\[Mu]_,{Bx_,By_,Bz_},spin_] := \[Mu](Bx Spin[spin,1] + By Spin[spin,2] + Bz Spin[spin,3])


End[];


(* ::Subsubsection::Closed:: *)
(*Hyperfine Hamiltonian*)


Begin["`Private`"];


HyperfineHamiltonian[spin1_,spin2_,A_]:=
	Sum[A[[i,j]]*Spin[spin1,i]\[CircleTimes]Spin[spin2,j],{i,3},{j,3}]


End[];


(* ::Subsubsection:: *)
(*Dipole Dipole Hamiltonian*)


Begin["`Private`"];


DipoleDipoleHamiltonian[spin1_,spin2_,\[Gamma]1_,\[Gamma]2_,v1_Vector,v2_Vector]:=
	With[{R=Sqrt[(v1-v2).(v1-v2)],e=Value@(v1-v2)},
		If[PossibleZeroQ[R],
			Message[DipoleDipoleHamiltonian::equallocations];
			ConstantArray[0,{SpinDim[spin1]*SpinDim[spin2],SpinDim[spin1]*SpinDim[spin2]}],
			((-\[Mu]0 \[Gamma]1 \[Gamma]2 \[HBar])/(4\[Pi] R^3))*(3(Total@(Spin[spin1]*e))\[CircleTimes](Total@(Spin[spin2]*e))/R^2-Spin[spin1,1]\[CircleTimes]Spin[spin2,1]-Spin[spin1,2]\[CircleTimes]Spin[spin2,2]-Spin[spin1,3]\[CircleTimes]Spin[spin2,3])
		]
	]


End[];


(* ::Subsubsection::Closed:: *)
(*Quadrapolar Hamiltonian*)


Begin["`Private`"];


QuadrapolarHamiltonian[spin_,A_?Matrix3Q]:=Total[Table[A[[i,j]]*Spin[spin,i].Spin[spin,j],{i,3},{j,3}],2];


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
				zeemanFrame = FrameCompose[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]],OptionValue[StaticFieldFrame]/.frames];,
			OptionValue[OutputFrame]===CrystalFrame,
				crystalFrame = IdentityFrame;
				labFrame = FrameInverse[OptionValue[CrystalOrientation]];
				zfsFrame = NVOrientationToFrame[OptionValue[NVOrientation]];
				zeemanFrame = FrameCompose[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]],OptionValue[StaticFieldFrame]/.frames];,
			OptionValue[OutputFrame]===ZFSFrame,
				zfsFrame = IdentityFrame;
				crystalFrame = FrameInverse[NVOrientationToFrame[OptionValue[NVOrientation]]];
				labFrame = FrameCompose[OptionValue[CrystalOrientation],crystalFrame];
				zeemanFrame = FrameCompose[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]],OptionValue[StaticFieldFrame]/.frames];,
			OptionValue[OutputFrame]===ZeemanFrame,
				zeemanFrame = IdentityFrame;
				(* The ZeemanFrame case is special because the StaticField has the option of being written in any frame. *)
				Which[
					OptionValue[StaticFieldFrame]===LabFrame,
						labFrame = FrameInverse[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]]];
						crystalFrame = FrameCompose[OptionValue[CrystalOrientation],labFrame];
						zfsFrame = FrameCompose[NVOrientationToFrame[OptionValue[NVOrientation]],crystalFrame];,
					OptionValue[StaticFieldFrame]===CrystalFrame,
						crystalFrame = FrameInverse[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]]];
						labFrame = FrameCompose[FrameInverse[OptionValue[CrystalOrientation]],crystalFrame];
						zfsFrame = FrameCompose[NVOrientationToFrame[OptionValue[NVOrientation]],crystalFrame];,
					OptionValue[StaticFieldFrame]===ZFSFrame,
						zfsFrame = FrameInverse[EulerAngles[0,sphericalB[[3]],sphericalB[[2]]]];
						crystalFrame = FrameCompose[FrameInverse[NVOrientationToFrame[OptionValue[NVOrientation]]],zfsFrame];
						labFrame = FrameCompose[FrameInverse[OptionValue[CrystalOrientation]],crystalFrame];
				];
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

		(* Start by adding the ZFS term*)
		AppendTo[termList,Term[ZFSHamiltonian[OptionValue[ZeroFieldSplitting],nvSpin],1]];

		(* Now add the NV Zeeman term *)
		AppendTo[termList,Term[ZeemanHamiltonian[\[Gamma]e,cartesianB,nvSpin],1]];

		(* Add all of the nuclei Zeeman terms *)
		NuclearZeemanTerm[nucleus_,index_]:=Term[ZeemanHamiltonian[GyromagneticRatio[nucleus],cartesianB,Spin[nucleus]],index+1];
		termList = Join[termList,MapIndexed[NuclearZeemanTerm[#1,First@#2]&, nucleiList]];

		(* Calculate and store all of the hyperfine interaction terms *)
		HyperfineTerm[nucleus_,index_]:=Term[HyperfineHamiltonian[nvSpin,Spin[nucleus],FrameChange[Tensor[nucleus],IdentityFrame,zfsFrame]], 1, index+1];
		termList = Join[termList,MapIndexed[HyperfineTerm[#1,First@#2]&, nucleiList]];

		(* Add the quadrapolar term for the Nitrogen *)
		If[hasN,
			AppendTo[termList,Term[QuadrapolarHamiltonian[Spin@First@nucleiList,QuadrapoleTensor@First@nucleiList],2]];
		];

		(* Add all of the dipole-dipole interactions *)
		DipoleDipoleTerm[nucleus1_,nucleus2_,index1_,index2_]:=
			Term[DipoleDipoleHamiltonian[Spin[nucleus1],Spin[nucleus2],GyromagneticRatio[nucleus1],GyromagneticRatio[nucleus2],\[Lambda]*Location[nucleus1],\[Lambda]*Location[nucleus2]],index1+1,index2+1];
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

		(* Not that it _really_ matters, but do the search and replace for constants before expanding to the
		   full Hilbert space because it should be a bit faster, at least, asymtotically.  *)
		If[numerical,
			termList = InsertConstants[termList];
		];

		(* Finally, expand all the terms and sum them up. *)
		(* We are carfeful to do this in a way that doesn't use more memory than needed. *)
		Fold[#1+preFactor*ExpandTerm[#2]&, 0, termList]
	]


End[];


(* ::Section::Closed:: *)
(*Plotting Tools*)


(* ::Subsection:: *)
(*Usage Declarations*)


Resolution::usage = "Resolution is an option for SpectrumData which decides when two frequencies are close enough to be considered the same frequency. The default value is $MachineEpsilon. Units of Resolution are MHz.";
SpectrumData::usage = "SpectrumData[nvHamiltonian] accepts a numerical NV Hamiltonian (where the first product space is that of the electron) and returns a list of the form {{freq1, rate1},{freq2,rate2},...}, containing the frequencies and respective rates of transitions under microwave radiation. The rates are normalized to sum to 1. The option Resolution is accepted, which decides when two freqencies are equal.";


SpectrumPlot::usage = "SpectrumPlot[nvhamiltonian,\[Sigma]] displays a simple plot of the spectrum of the input numerical Hamiltonian.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


Options[SpectrumData]={Resolution->0.001, AngularUnits->True};
SpectrumData[ham_,OptionsPattern[]]:=
	Module[{U,D,C,probs,freqs,output},
		If[PossibleZeroQ[Mod[Length[ham],3]],
			C = Spin[1,1]\[CircleTimes]IdentityMatrix[Length[ham]/3];,
			C = Spin[1/2,1]\[CircleTimes]IdentityMatrix[Length[ham]/2];
		];

		{D, U} = Eigensystem[ham/10^6];
		U = Normalize/@U;

		probs = Abs[Flatten[UpperTriangularize[U.C.U\[ConjugateTranspose]]]]^2;
		probs = probs/Total[probs];
		freqs = Abs@Flatten@UpperTriangularize@Outer[Subtract,D,D,1];_

		If[OptionValue[AngularUnits],freqs=freqs/(2\[Pi])];

		output = Select[{freqs,probs}\[Transpose],(Last@#>OptionValue[Resolution])&];
		output = Sort[output,First@#1>First@#2&];
		output = Split[output,Abs[First@#1-First@#2]<=OptionValue[Resolution]&];
		output = Map[{Mean[#[[All,1]]],Total[#[[All,2]]]}&,output,1];
		output
	]


Options[SpectrumPlot]=Join[Options[Plot],Options[SpectrumData],{}];
SpectrumPlot[ham_,\[Sigma]_,opt:OptionsPattern[]]:=
	Module[{spectData,spectrum,minFreq,maxFreq,freqDiff,minPlot,maxPlot,plotOptions},
		spectData = SpectrumData[ham,Resolution->OptionValue[Resolution],AngularUnits->OptionValue[AngularUnits]];
		spectrum[f_]=Total[Last[#]*Exp[-(f-First[#])^2/(2*\[Sigma]^2)]&/@spectData];
		minFreq = First@Last@spectData;
		maxFreq = First@First@spectData;
		freqDiff = maxFreq - minFreq;
		minPlot = minFreq - Max[0.1*freqDiff,5*\[Sigma]];
		maxPlot = maxFreq + Max[0.1*freqDiff,5*\[Sigma]];
		plotOptions = Sequence@@Select[{opt},MemberQ[(Options[Plot])[[All,1]],First[#]]&];
		Plot[spectrum[f],{f,minPlot,maxPlot},AxesOrigin->{minPlot,0},Evaluate@plotOptions]
	]


End[];


(* ::Section:: *)
(*Approximation Tools*)


(* ::Subsection:: *)
(*Usage Declarations*)


MagnusExpansionTerm::usage="MagnusExpansionTerm[{At,t,T},k] Returns the k-th term in the magnus expansions (k=1,2,..),
of the time dependent operator At with time parameter t, integrated over the time period 0 to T.";
MagnusExpansion::usage="MagnusExpansion[{At,T,T},k] returns a list of the first k (k=1,2,3,..) terms of the Magnus Expansion of time dependent operator At with time paramter t, intergrated over the period 0 to T.";


AverageHamiltonianTerm::usage="AverageHamiltonianTerm[{Ht,t,T},k] Returns the k-th term in the average hamiltonian expansion (k=0,1,2,..),
of the time dependent Hamiltonian Ht with time parameter t, integrated over the time period 0 to T.";
AverageHamiltonian::usage="AverageHamiltonian[{Ht,T,T},k] returns the k-th order (k=0,1,2,3,..) Average Hamiltonian expansions of the time dependent Hamiltonian Ht with time paramter t, intergrated over the period 0 to T.";


MagnusConvergenceTest::usage="MagnusConvergenceTest[{H,t,T}] performs the integral \!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(0\)], \(T\)]\)\[RightDoubleBracketingBar]H[t]\[RightDoubleBracketingBar]dt where \[RightDoubleBracketingBar]\[CenterDot]\[RightDoubleBracketingBar] is the operator norm. If this quantity is less than \[Pi], then the Magnus series converges.";


NVAverageHamiltonian::usage="NVAverageHamiltonian[order,...] returns the order'th order average NV Hamiltonian. Here, '...' represents any arguments you would pass into the function NVHamiltonian. If order is negative, the effective Hamiltonian is returned.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*The following is a recursive formula for computing the magnus expanion to arbitrary order. It works with symbolic or numberic matrices, though it can get quite slow for higher orders with symbolic matrices. The recursive formula is based off the one given on the Wikipedia page for the Magnus Expansion. *)
(*We use the convention that the Magnus Expansion states with first term k=1, while the Average Hamiltonian starts with first term k=0, and divided by the appropriate factor of i and the intergration time to make it a Hamiltonian. *)


MagnusGenerator[{At_,t_,T_},{n_,1},fn_:Identity]:=
	MagnusGenerator[{At,t,T},{n,1},fn]=
		fn@Com[MagnusExpansionTerm[{At,t,T},n-1,fn],At/.{t->T}]

MagnusGenerator[{At_,t_,T_},{n_,j_},fn_:Identity]:=
	MagnusGenerator[{At,t,T},{n,j},fn]=
		fn@If[j===n-1,
			NestCom[MagnusExpansionTerm[{At,t,T},1,fn],(At/.{t->T}),n-1],
			Total[
				Com[MagnusExpansionTerm[{At,t,T},#,fn],MagnusGenerator[{At,t,T},{n-#,j-1},fn]]&/@Range[n-j]
			]
		]

MagnusExpansionTerm[{At_,t_,T_},1,fn_:Identity]:=
	MagnusExpansionTerm[{At,t,T},1,fn]=
		fn@Integrate[At,{t,0,T}]

MagnusExpansionTerm[{At_,t_,T_},k_,fn_:Identity]:=
	MagnusExpansionTerm[{At,t,T},k,fn]=
		fn@Total[BernoulliB[#]*Integrate[MagnusGenerator[{At,t,t1},{k,#},fn],{t1,0,T}]/#!&/@Range[k-1]]

MagnusExpansion[{At_,t_,T_},order_]:=MagnusExpansionTerm[{At,t,T},#]&/@Range[order];


(* ::Text:: *)
(*The following allows for optional arguments which can speed up symbolic calculations, in general using the Simply->True option will be the fastest with symbolic matrices.*)


MagnusExpansion[{At_,t_,T_},order_,Simplify->True]:=MagnusExpansionTerm[{At,t,T},#,Simplify]&/@Range[order];
MagnusExpansion[{At_,t_,T_},order_,FullSimplify->True]:=MagnusExpansionTerm[{At,t,T},#,FullSimplify]&/@Range[order];
MagnusExpansion[{At_,t_,T_},order_,Chop->True]:=MagnusExpansionTerm[{At,t,T},#,Chop]&/@Range[order];


(* ::Text:: *)
(*Finally we include an AverageHamiltonian function which sums up the magnus terms, correctly renormalized to be the average Hamiltonian over the period T. Further we use the Physics convention of starting the AverageHamiltonian with 0th order term, instead of the lowest order being 1st as in the Magnus Expansion.*)


AverageHamiltonianTerm[{Ht_,t_,T_},k_,fn_:Identity]:=(-I)^k*MagnusExpansionTerm[{Ht,t,T},k+1,fn]/T

AverageHamiltonian[{Ht_,t_,T_},order_]:=Total[Array[(-I)^(#-1)/T&,order+1]*MagnusExpansion[{Ht,t,T},order+1]]
AverageHamiltonian[{Ht_,t_,T_},order_,fn_->True]:=Total[Array[(-I)^(#-1)/T&,order+1]*MagnusExpansion[{Ht,t,T},order+1,fn->True]]


MagnusConvergenceTest[{At_,t_,T_}]:=Integrate[Abs[Eigenvalues[At,1]],{t,0,T}][[1]]


NVAverageHamiltonian[order_,nuclei___,opt:OptionsPattern[NVHamiltonian]]:=
	Module[
		{H,H0,U,ndim,sgn,t,zfs,numerical,angularUnits,fn,Heff},

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
		zfs = OptionValue[ZeroFieldSplitting];
		If[numerical,zfs=InsertConstants[zfs]];
		If[angularUnits,zfs=2*\[Pi]*zfs];

		(* Decide whether to Simplify or Chop between each magnus order calculation*)
		If[numerical,fn=Simplify,fn=Chop];

		(* Find the change of frame required *)
		H0 = ZFSHamiltonian[zfs, OptionValue[NVSpin]]\[CircleTimes]IdentityMatrix[ndim];
		U[sgn_] = MatrixExp[sgn I t H0];

		(* Finally, do the actual calculation *)
		If[order < 0,
			Heff[s_]=U[1].(H-H0).U[-1]/.t->s;Heff,
			AverageHamiltonian[{U[1].(H-H0).U[-1],t,2\[Pi]/zfs},order,fn->True]
		]
	]


End[];


(* ::Section::Closed:: *)
(*Epilogue*)


(* ::Text:: *)
(*This section exists so that EndPackage[] isn't just placed in the last section. Not doing this can result in very confusing behaviour when you add another section but forget to move EndPackage[].*)


EndPackage[];
