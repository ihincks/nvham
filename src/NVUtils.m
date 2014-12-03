(* ::Package:: *)

BeginPackage["NVUtils`"];


(* ::Text:: *)
(*This package contains mostly duplicate functions from QuantumUtils` which are needed in the nvsim software package. Do not include this package if you are using QuantumUtils`, but you will need to include this package otherwise. Needs should be called on this package before Needs is called on NVHamiltonian`.*)


(* ::Section::Closed:: *)
(*Error Messages*)


(* ::Subsection:: *)
(*Matrices, Bases, and Linear Algebra*)


Spin::badindex = "An invalid spin index was entered.";


IdentityInsert::baddimensions = "The size of the input matrix does not match the product of the specified dimensions.";


(* ::Section:: *)
(*Matrices, Bases, and Linear Algebra*)


(* ::Subsection:: *)
(*Usage Declarations*)


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


(* ::Subsubsection:: *)
(*Spin operators*)


Begin["`Private`"];


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


(* ::Subsubsection:: *)
(*Tensor Manipulate*)


Begin["`Private`"];


(* ::Text:: *)
(*First deal with the most difficult case where we want to insert an identity operation into the middle.*)


IdentityInsert[C_,dimA_,dimB_,n_]:=
	If[Length[C]!=dimA*dimB,
		Message[IdentityInsert::baddimensions];Abort;,
		ArrayFlatten[Map[IdentityMatrix[n]\[CircleTimes]#&,Partition[C,{dimB,dimB}],{2}]]
	]


(* ::Text:: *)
(*Now expand the definition to the left and right sides too.*)


IdentityInsert[C_,dimA_,dimB_,n1_,n2_,n3_]:=IdentityMatrix[n1]\[CircleTimes]IdentityInsert[C,dimA,dimB,n2]\[CircleTimes]IdentityMatrix[n3]


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


(* ::Section::Closed:: *)
(*Epilog*)


EndPackage[];
