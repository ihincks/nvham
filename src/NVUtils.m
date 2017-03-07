(* ::Package:: *)

BeginPackage["NVUtils`"];


(* ::Text:: *)
(*This package contains mostly duplicate functions from QuantumUtils` which are needed in the nvsim software package. Do not include this package if you are using QuantumUtils`, but you will need to include this package otherwise. Needs should be called on this package before Needs is called on NVHamiltonian`.*)


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
