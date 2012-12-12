(* ::Package:: *)

BeginPackage["NVSim`", {"QuantumUtils`","Predicates`","NVHamiltonian`"}]


(* ::Section::Closed:: *)
(*Predicates*)


(* ::Subsection:: *)
(*Usage Declarations*)


PulseShapeFileQ::usage = "PulseShapeFileQ[str] returns True iff str is a string pointing to a text file containg a pulse.";
PulseShapeMatrixQ::usage = "PulseShapeMatrixQ[M] returns True iff M is a 2D matrix";
PulseShapeQ::usage = "PulseShapeQ[in] returns True iff one of PulseShapeFileQ or PulseShapeMatrixQ is True";


ShapedPulseQ::usage = "ShapedPulse[p] returns True iff p is of the form {pulse,{Hcontrol1,Hcontrol2,..}} where the Hcontrols are the control Hamiltonians and pulse satisfies PulseShapeQ.";
DriftPulseQ::usage = "DriftPulseQ[p] returns True iff p is a real number indicating the amount of time to evolve under the drift Hamiltonian.";
InstantaneousPulseQ::usage = "InstantaneousPulseQ[p] returns True iff p is a matrix.";
PulseQ::usage = "PulseQ[p] returns True iff p satisfies at least one of the NVSim predicates ending in PulseQ (ShapedPulseQ,InstantaneousPulseQ,etc)";


PulseSequenceQ::usage = "PulseSequenceQ[seq] returns True iff seq is a list where each element satisfies PulseQ.";


DriftHamConstQ::usage = "DriftHamConstQ[H] returns True iff H is a square matrix.";
DriftHamNonConstQ::usage = "DriftHamNonConst[H] returns True iff H is a function which accepts a number and returns a square matrix.";
DriftHamQ::usage = "DriftHamQ[H] returns True iff one of DriftHamConstQ[H] or DriftHamNonConstQ[H] is True.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


PulseShapeFileQ[str_]:=StringQ[str]


PulseShapeMatrixQ[M_]:=MatrixQ[M]


PulseShapeQ[in_]:=PulseShapeFileQ[in]||PulseShapeMatrixQ[in]


ShapedPulseQ[p_]:=ListQ[p]&&(Length[p]==2)&&PulseShapeQ[p[[1]]]


DriftPulseQ[p_]:=NumericQ[p]&&p\[Element]Reals


InstantaneousPulseQ[p_]:=SquareMatrixQ[p]


PulseQ[p_]:=Or@@(Through[{ShapedPulseQ,DriftPulseQ,InstantaneousPulseQ}[p]])


PulseSequenceQ[seq_]:=ListQ[seq]&&(And@@(PulseQ/@seq))


DriftHamConstQ[H_]:=SquareMatrixQ[H]


DriftHamNonConstQ[H_]:=Module[{t},SquareMatrixQ[H[t]]]


DriftHamQ[H_]:=DriftHamConstQ[H]||DriftHamNonConstQ[H]


End[];


(* ::Section::Closed:: *)
(*Options and Helper Functions*)


(* ::Subsection:: *)
(*Usage Declarations*)


PulseOptions::usage = "PulseOptions is a dummy function which stores the options for various NVSim functions. Use the command Options[PulseOptions] to view these options and their default values.";


SimTrace::usage = "SimTrace is a simulation option that chooses which in-sequence information of the pulse sequence to store.";
SimMethod::usage = "SimMethod is a simulation option that chooses which method of simulation to use (Automatic, Naive, Compiled, CompiledC, etc.)";
StepSize::usage = "StepSize is a simulation option that chooses the time discretization when the internal Hamiltonian is time dependent. Can be set to Automatic.";
IntegrationMethod::usage = "IntegrationMethod is a simulation option that chooses which integration method to use in the case of a time dependent Hamiltonian.";


CompiledC::usage = "CompiledC is a SimMethod option value.";
Naive::usage = "Naive is a SimMethod option value.";


GetPulseShapeMatrix::usage = "GetPulseShapeMatrix[in] returns in if in is a matrix, but if in is a file name, returns the contents of that file as a matrix.";
GetStepSize::usage = "GetStepSize[H,stepsize:Automatic] returns a fifth of the biggest element of H[0] if stepsize is Automatic, and stepsize otherwise.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


Protect[SimTrace,SimMethod,StepSize,IntegrationMethod,CompiledC,Naive];


Options[PulseOptions]={
	SimTrace->False,
	SimMethod->Automatic,
	StepSize->Automatic,
	IntegrationMethod->Automatic
};


GetPulseShapeMatrix[in_?PulseShapeFileQ]:=With[{out=Import[in]//N},Pick[out,Length[#]>1&/@out]]
GetPulseShapeMatrix[in_?PulseShapeMatrixQ]:=in//N


(* ::Text:: *)
(*If the step size is set to Automatic, let it be a tenth of the biggest element of H maximized over time, otherwise, the user has given the stepsize.*)


GetStepSize[H_,stepsize_:Automatic]:=
	If[stepsize===Automatic,
		Module[{t,max},
			max=NMaximize[Max[H[t]],t][[1]];
			1/(10*max)
		],
		stepsize
	]


Else=True;


End[];


(* ::Section::Closed:: *)
(*Barebone Simulators*)


(* ::Subsection:: *)
(*Usage Declarations*)


DotExpNaive::usage = "DotExpNaive=Function[{Hint,Hctls,dt,amps},expr]";
DotExpCompiled::usage = "DotExpCompiled=Compile[{{Hint,_Complex,2},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},expr]";
DotExpCompiledC::usage = "DotExpCompiledC=Compile[{{Hint,_Complex,2},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},expr]";


TimeDepDriftNaive::usage = "TimeDepDriftNaive=Function[{T,dt},expr]";
TimeDepDriftCompiled::usage = "TimeDepDriftCompiled=Compile[{{T,_Real},{dt,_Real}},expr]";
TimeDepDriftCompiledC::usage = "TimeDepDriftCompiledC=Compile[{{T,_Real},{dt,_Real}},expr]";


DotExpTimeDepNaive::usage = "DotExpTimeDepNaive=Function[{subdt,Hctls,dt,amps},expr]";
DotExpTimeDepCompiled::usage "DotExpTimeDepCompiled=Compile[{{subdt,_Real},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},expr]";
DotExpTimeDepCompiledC::usage "DotExpTimeDepCompiledC=Compile[{{subdt,_Real},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},expr]";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Declare a private symbol. We will set this symbol to be time dependent Hamiltonians so that we can use them inside of compiled functions.*)


Hfun::usage = "A private variable storing the current time dependent Hamiltonian";


(* ::Subsubsection::Closed:: *)
(*DotExp functions*)


(* ::Text:: *)
(*Takes a time independent drift Hamiltonian, the control Hamiltonians, the time steps, and the amplitudes and multiplies all of the matrix exponentials together, just using the usual Fold method.*)


DotExpNaive=Function[{Hint,Hctls,dt,amps},
	With[{ampind=If[Length[amps[[1]]]>1,{2,-1},{2}]},
		Fold[
			(MatrixExp[-I*(#2[[1]])*(Hint+Total[Hctls*(#2[[ampind]])])].#1)&,
			IdentityMatrix[Length[Hint]],
			Join[{dt},amps\[Transpose]]\[Transpose]
		]
	]
];


(* ::Text:: *)
(*Takes a time independent drift Hamiltonian, the control Hamiltonians, the time steps, and the amplitudes and multiplies all of the matrix exponentials together, using a compiled function.*)


DotExpCompiled=Compile[{{Hint,_Complex,2},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},
	With[{num=Length[dt],dim=Length[Hint]},
		Module[{out=(1.0+0.0 I)*IdentityMatrix[dim],H=Table[0.0 I,{i,dim},{j,dim}],k},
			For[k=1,k<=num,k++,
				H=-I*dt[[k]]*(Hint+Total[Hctls*amps[[k,All]]]);
				out=MatrixExp[H].out;
			];
			out
		]
	]
];


(* ::Text:: *)
(*Same as the above function, except compiles to C. This should be faster in the limit of larger Hilbert space dimension, but has not been tested yet.*)


DotExpCompiledC=Compile[{{Hint,_Complex,2},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},
	With[{num=Length[dt],dim=Length[Hint]},
		Module[{out=(1.0+0.0 I)*IdentityMatrix[dim],H=Table[0.0 I,{i,dim},{j,dim}],k},
			For[k=1,k<=num,k++,
				H=-I*dt[[k]]*(Hint+Total[Hctls*amps[[k,All]]]);
				out=MatrixExp[H].out;
			];
			out
		]
	],
	CompilationTarget->"C"
];


(* ::Subsubsection::Closed:: *)
(*TimeDepDrift functions*)


(* ::Text:: *)
(*Compute the propagator of the time dependent Hamiltonian Hfun for a total time t and stepsize dt. Use the midpoint method.*)


TimeDepDriftNaive=Function[{T,dt,dim},
	Module[{num=Floor[T/dt],out},
		out=Fold[
			(MatrixExp[-I*dt*NVSim`Private`Hfun[(#2-0.5)*dt]].#1)&,
			IdentityMatrix[dim],
			Table[n,{n,num}]
		];
		out=MatrixExp[-I*(T-num*dt)*NVSim`Private`Hfun[0.5*(T+num*dt)]].out;
		out
	]
];


TimeDepDriftCompiled=Compile[{{T,_Real},{dt,_Real},{dim,_Integer}},
	Module[{num,out,H,k},
		num = Floor[T/dt];
		out = (1.0+0.0 I)*IdentityMatrix[dim];
		H=ConstantArray[0.0 I,{dim,dim}];
		For[k=1,k<=num,k++,
			H=-I*dt*NVSim`Private`Hfun[(k-0.5)*dt];
			out=MatrixExp[H].out
		];
		H=-I*(T-num*dt)*NVSim`Private`Hfun[0.5*(T+num*dt)];
		out = MatrixExp[H].out;
		out
	],
	{{NVSim`Private`Hfun[_],_Complex,2}}
];


TimeDepDriftCompiledC=Compile[{{T,_Real},{dt,_Real},{dim,_Integer}},
	Module[{num,out,H,k},
		num = Floor[T/dt];
		out = (1.0+0.0 I)*IdentityMatrix[dim];
		H=ConstantArray[0.0 I,{dim,dim}];
		For[k=1,k<=num,k++,
			H=-I*dt*NVSim`Private`Hfun[(k-0.5)*dt];
			out=MatrixExp[H].out
		];
		H=-I*(T-num*dt)*NVSim`Private`Hfun[0.5*(T+num*dt)];
		out = MatrixExp[H].out;
		out
	],
	{{NVSim`Private`Hfun[_],_Complex,2}},
	CompilationTarget->"C"
];


(* ::Subsubsection:: *)
(*DotExpTimeDep functions*)


(* ::Text:: *)
(*DotExpTimeDep functions is essentially a TimeDepDrift function inserted into the for loop of a DotExp function: on each control step, the evolution is broken down into steps of length subdt*)


(* TODO: turn this into a Fold or something *)
DotExpTimeDepNaive=Function[{subdt,Hctls,dt,amps},
	Module[{dim,num,subnum,out,H,Hctl,T,k,l},
		dim=Length[Hctls[[1]]];
		num=Length[dt];
		out=(1.0+0.0 I)*IdentityMatrix[dim];
		H=ConstantArray[0.0 I,{dim,dim}];
		(* In each control step the controls do not change, and so we store it in the variable Hctl *)
		Hctl=ConstantArray[0.0 I,{dim,dim}];
		(* T records the time at the end of the previous control step *)
		T=0.0;
		(* k loops through the control sequuence *)
		For[k=1,k<=num,k++,
			subnum = Floor[dt[[k]]/subdt];
			Hctl=Total[Hctls*amps[[k,All]]];
			(* l loops through substeps to account for drift Hamiltonian time dependence *)
			For[l=1,l<=subnum,l++,
				H=-I*subdt*(NVSim`Private`Hfun[T+(l-0.5)*subdt]+Hctl);
				out=MatrixExp[H].out
			];
			(* Tack on the remaining time excluded from the above loop *)
			H=-I*(dt[[k]]-subnum*subdt)*(NVSim`Private`Hfun[T+0.5*(dt[[k]]+subnum*subdt)]+Hctl);
			out=MatrixExp[H].out;
			T=T+dt[[k]]
		];
		out
	]
];


DotExpTimeDepCompiled=Compile[{{subdt,_Real},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},
	Module[{dim,num,subnum,out,H,Hctl,T,k,l},
		dim=Length[Hctls[[1]]];
		num=Length[dt];
		out=(1.0+0.0 I)*IdentityMatrix[dim];
		H=ConstantArray[0.0 I,{dim,dim}];
		(* In each control step the controls do not change, and so we store it in the variable Hctl *)
		Hctl=ConstantArray[0.0 I,{dim,dim}];
		(* T records the time at the end of the previous control step *)
		T=0.0;
		(* k loops through the control sequuence *)
		For[k=1,k<=num,k++,
			subnum = Floor[dt[[k]]/subdt];
			Hctl=Total[Hctls*amps[[k,All]]];
			(* l loops through substeps to account for drift Hamiltonian time dependence *)
			For[l=1,l<=subnum,l++,
				H=-I*subdt*(NVSim`Private`Hfun[T+(l-0.5)*subdt]+Hctl);
				out=MatrixExp[H].out
			];
			(* Tack on the remaining time excluded from the above loop *)
			H=-I*(dt[[k]]-subnum*subdt)*(NVSim`Private`Hfun[T+0.5*(dt[[k]]+subnum*subdt)]+Hctl);
			out=MatrixExp[H].out;
			T=T+dt[[k]]
		];
		out
	],
	{{NVSim`Private`Hfun[_],_Complex,2}}
];


DotExpTimeDepCompiledC=Compile[{{subdt,_Real},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},
	Module[{dim,num,subnum,out,H,Hctl,T,k,l},
		dim=Length[Hctls[[1]]];
		num=Length[dt];
		out=(1.0+0.0 I)*IdentityMatrix[dim];
		H=ConstantArray[0.0 I,{dim,dim}];
		(* In each control step the controls do not change, and so we store it in the variable Hctl *)
		Hctl=ConstantArray[0.0 I,{dim,dim}];
		(* T records the time at the end of the previous control step *)
		T=0.0;
		(* k loops through the control sequuence *)
		For[k=1,k<=num,k++,
			subnum = Floor[dt[[k]]/subdt];
			Hctl=Total[Hctls*amps[[k,All]]];
			(* l loops through substeps to account for drift Hamiltonian time dependence *)
			For[l=1,l<=subnum,l++,
				H=-I*subdt*(NVSim`Private`Hfun[T+(l-0.5)*subdt]+Hctl);
				out=MatrixExp[H].out
			];
			(* Tack on the remaining time excluded from the above loop *)
			H=-I*(dt[[k]]-subnum*subdt)*(NVSim`Private`Hfun[T+0.5*(dt[[k]]+subnum*subdt)]+Hctl);
			out=MatrixExp[H].out;
			T=T+dt[[k]]
		];
		out
	],
	{{NVSim`Private`Hfun[_],_Complex,2}},
	CompilationTarget->"C"
];


(* ::Subsubsection::Closed:: *)
(*Epilog*)


End[];


(* ::Section::Closed:: *)
(*Single Pulse Evaluator*)


(* ::Subsection:: *)
(*TODO*)


(* ::Text:: *)
(*Complete SimTrace=True options*)


(* ::Text:: *)
(*Add new SimTrace option: hermitian observable list*)


(* ::Subsection:: *)
(*Usage Declarations*)


EvalPulse::usage = "EvalPulse[]";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Subsubsection::Closed:: *)
(*Shaped Pulse Evaluators*)


EvalPulse[H_?DriftHamConstQ,p_?ShapedPulseQ,OptionsPattern[PulseOptions]]:=
	Module[{simFunction,pulse,dt,Hctls,amps,Uj},
		pulse=GetPulseShapeMatrix[p[[1]]];
		dt = pulse[[All,1]];
		amps = If[Length[pulse[[1]]]>2,pulse[[All,{2,-1}]],pulse[[All,{2}]]];
		Hctls = p[[2]];
		Which[
			OptionValue[SimTrace]==False,
				simFunction=OptionValue[SimMethod]/.{Automatic->DotExpCompiled,Compiled->DotExpCompiled,CompiledC->DotExpCompiledC,Naive->DotExpNaive};									
				Apply[simFunction,{H,Hctls,dt,amps}],
			OptionValue[SimTrace]==True,
				Uj=Table[MatrixExp[-I*dt[[k]]*(H+Total[Hctls*amps[[k]]])],{k,Length[dt]}];
				FoldList[ #2.#1&,IdentityMatrix[Dimensions[H]],Uj]
		]
	]


EvalPulse[H_?DriftHamNonConstQ,p_?ShapedPulseQ,OptionsPattern[PulseOptions]]:=
	Module[{simFunction,pulse,dt,d\[Tau],Hctls,amps},
		pulse=GetPulseShapeMatrix[p[[1]]];
		dt = pulse[[All,1]];
		amps = If[Length[pulse[[1]]]>2,pulse[[All,{2,-1}]],pulse[[All,{2}]]];
		Hctls = p[[2]];
		NVSim`Private`Hfun = H;
		d\[Tau]=GetStepSize[H,OptionValue[StepSize]];
		Which[
			OptionValue[SimTrace]==False,
				simFunction=OptionValue[SimMethod]/.{Automatic->DotExpTimeDepCompiled,Compiled->DotExpTimeDepCompiled,CompiledC->DotExpTimeDepCompiledC,Naive->DotExpTimeDepNaive};
				Apply[simFunction,{d\[Tau],Hctls,dt,amps}],
			OptionValue[SimTrace]==True,
				{}
		]
	]


(* ::Subsubsection::Closed:: *)
(*Instantaneous Pulse Evaluators*)


EvalPulse[H_?DriftHamQ,p_?InstantaneousPulseQ,OptionsPattern[PulseOptions]]:=
	With[{st=OptionValue[SimTrace]}, 
		Which[
			st==False,
				p,
			st==True,
				{p}
		]
	]


(* ::Subsubsection::Closed:: *)
(*Drift Pulse Evaluators*)


EvalPulse[H_?DriftHamConstQ,p_?DriftPulseQ,OptionsPattern[PulseOptions]]:=
	With[{U=MatrixExp[-I*p*H]},
		Which[
			OptionValue[SimTrace]===False,
				U,
			OptionValue[SimTrace]===True,
				{U}
		]
	]


EvalPulse[H_?DriftHamNonConstQ,p_?DriftPulseQ,OptionsPattern[PulseOptions]]:=
	Module[{simFunction,d\[Tau],T},
		d\[Tau]=GetStepSize[H,OptionValue[StepSize]];
		T = p;
		NVSim`Private`Hfun = H;
		Which[
			OptionValue[SimTrace]===False,
				simFunction=OptionValue[SimMethod]/.{Automatic->TimeDepDriftCompiled,Compiled->TimeDepDriftCompiled,CompiledC->TimeDepDriftCompiledC,Naive->TimeDepDriftNaive};
				Apply[simFunction, {T,d\[Tau],Length[H[0]]}],
			OptionValue[SimTrace]===True,
				{}
		]
	]


(* ::Subsubsection::Closed:: *)
(*Epilog*)


End[];


(* ::Section::Closed:: *)
(*Pulse Sequence Evaluator*)


(* ::Subsection:: *)
(*UsageDeclarations*)


EvalPulseSequence::usage = "EvalPulseSequence[]";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


EvalPulseSequence[H_?DriftHamQ,seq_?PulseSequenceQ,options:OptionsPattern[PulseOptions]]:=
	Module[{allPulses},
		allPulses=EvalPulse[H,#,options]&/@seq;
		Which[
			OptionValue[SimTrace]===False,
				Dot@@Reverse[allPulses],
			OptionValue[SimTrace]===True,
				Flatten[allPulses,1]
		]
	]


End[];


(* ::Section::Closed:: *)
(*Sequence Drawing*)


(* ::Subsection:: *)
(*Usage Declarations*)


DrawSequence::usage = "DrawSequence[seq] outputs a grahpical representation of the pulse sequence seq.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*A helper function to print out pulse times in appropriate units.*)


ShowTime[t_?NumberQ]:=
	With[{O=Function[{x},ToString[Round[x]]]},
		Which[
			10^-3<=t<1,O[t*10^3]<>"ms",
			10^-6<=t<10^-3,O[t*10^6]<>"\[Mu]s",
			10^-9<=t<10^-6,O[t*10^9]<>"ns",
			10^-12<=t<10^-9,O[t*10^12]<>"ps",
			True,ToString[PaddedForm[t,{3,1}]]<>"s"
		]
	]


(* ::Text:: *)
(*For shaped pulses we just pick a nice looking shape instead of using the actual data.*)


DrawPulse[p_?ShapedPulseQ,width_,height_,offset_]:=
	Module[{minx,maxx,n,d,data,totaltime},
		minx = -2.1;
		maxx = 3;
		n = 30;
		d=(maxx-minx)/(n-1);
		data = (Exp[-(#-1)^2]+Exp[-3*(#+1)^2]/2)&/@Range[minx,maxx,d];
		data = height*data/Max[data];
		d = width/n;
		totaltime = Total[GetPulseShapeMatrix[p[[1]]][[All,1]]];
		{
			Table[Line[{{(k-1)*d+offset,data[[k]]},{k*d+offset,data[[k]]}}],{k,n}],
			Table[Line[{{k*d+offset,0},{k*d+offset,Max[data[[k]],data[[k+1]]]}}],{k,1,n-1}],
			Text[p[[1]],{d*n/2+offset,-height/8}],
			Text[ShowTime[totaltime],{d*n/2+offset,-height/4}]
		}
	]


DrawPulse[p_?InstantaneousPulseQ,width_,height_,offset_]:=
	{
		Line[{{offset,0},{offset,height},{width+offset,height},{width+offset,0}}],
		Text[If[Length[p]<3,p//N,"Instant Pulse"],{width/2+offset,-height/6}]
	}


DrawPulse[p_?DriftPulseQ,width_,height_,offset_]:=
	{
		Text["\!\(\*SubscriptBox[\(\[ScriptCapitalH]\), \(drift\)]\)",{width/2+offset,-height/8}],
		Text[ShowTime[p],{width/2+offset,-height/4}]
	}


(* ::Text:: *)
(*We give each kind of pulse its own width weight for aesthetics.*)


DrawSequence[seq_?PulseSequenceQ]:=
	Module[{shapedFrac=0.3,instFrac=0.07,driftFrac=0.63,width=500,height=100,widths},
		widths=shapedFrac*(ShapedPulseQ/@seq)/.{True->1,False->0};
		widths=widths+instFrac*(InstantaneousPulseQ/@seq)/.{True->1,False->0};
		widths=widths+driftFrac*(DriftPulseQ/@seq)/.{True->1,False->0};
		widths=width*widths/Total[widths];
		Graphics[{
			Arrowheads[0.02],
			Arrow[{{0,0},{width,0}},-width/20],
			Text["t",{width+width/15,0}],
			Table[DrawPulse[seq[[k]],widths[[k]],height*0.8,Total[Take[widths,k]]-widths[[k]]],{k,Length[seq]}]
		}]
	]


End[];


(* ::Section::Closed:: *)
(*Epilog*)


EndPackage[];
