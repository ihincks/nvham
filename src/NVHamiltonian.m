(* ::Package:: *)

BeginPackage["NVHamiltonian`"];


(* ::Section:: *)
(*Error Messages*)


(* ::Subsection:: *)
(*Matrices, Bases, and Linear Algebra*)


Spin::badindex = "An invalid spin index was entered.";


IdentityInsert::baddimensions = "The size of the input matrix does not match the product of the specified dimensions.";


(* ::Section:: *)
(*Predicates*)


(* ::Subsection:: *)
(*Usage Declarations*)


Carbon::usage = ""; Nitrogen::usage = "";
CarbonQ::usage = "CarbonQ[c] returns True iff the Head of c is Carbon";
NitrogenQ::usage = "NitrogenQ[c] returns True iff the Head of c is Carbon";
NucleusQ::usage = "CarbonQ[c] returns True iff one of CarbonQ or NitrogenQ holds.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


CarbonQ[c_]:=Head[c]===Carbon
NitrogenQ[c_]:=Head[c]===Nitrogen
NucleusQ[c_]:=CarbonQ[c]||NitrogenQ[c]


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


Frame::usage = "";
LabFrame::usage = "";
ZFSFrame::usage = "";
CrystalFrame::usage =  "";
ZeemanFrame::usage = "";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


End[];


(* ::Section:: *)
(*Options*)


(* ::Subsection:: *)
(*Usage Declarations*)


NVHamiltonian::options = "";
ZeroFieldSplitting::usage = "";
StaticField::usage = "";
CrystalOrientation::usage = "";
NVSpin::usage = "";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


Options[NVHamiltonian]=
	{
		ZeroFieldSplitting->\[CapitalDelta],
		StaticField->{0,0,0},
		CrystalOrientation->{0,0,0},
		NVSpin->1
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
		{dimNV,dimN,hasN,dimC,numC,dimList,nucleiList,nvSpin,termList,Term,ExpandTerm,HyperfineTerm},

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
