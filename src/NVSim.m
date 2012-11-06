(* ::Package:: *)

BeginPackage["NVSim`", {"QuantumUtils`","Predicates`","NVHamiltonian`"}]


(* ::Section:: *)
(*Predicates*)


(* ::Subsection:: *)
(*Usage Declarations*)


PulseShapeFileQ::usage = "PulseShapeFileQ[str] returns True iff str is a string pointing to a text file containg a pulse.";
PulseShapeMatrixQ::usage = "PulseShapeMatrixQ[M] returns True iff M is a 2D matrix";
PulseShapeQ::usage = "PulseShapeQ[in] returns True iff one of PulseShapeFileQ or PulseShapeMatrixQ is True";


ShapedPulseQ::usage = "ShapedPulse[p] returns True iff p is of the form {filename,{Hcontrol1,Hcontrol2,..}} where the Hcontrols are the control Hamiltonians and filename is a string specifying the location of a pulse file.";
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


DriftPulseQ[p_]:=NumberQ[p]&&p\[Element]Reals


InstantaneousPulseQ[p_]:=SquareMatrixQ[p]


PulseQ[p_]:=Or@@(Through[{ShapedPulseQ,DriftPulseQ,InstantaneousPulseQ}[p]])


PulseSequenceQ[seq_]:=ListQ[seq]&&(And@@(PulseQ/@seq))


DriftHamConstQ[H_]:=SquareMatrixQ[H]


DriftHamNonConstQ[H_]:=Block[{t},SquareMatrixQ[H[t]]]


DriftHamQ[H_]:=DriftHamConstQ[H]||DriftHamNonConstQ[H]


End[];


(* ::Section:: *)
(*Options and Helper Functions*)


(* ::Subsection:: *)
(*Usage Declarations*)


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
(*If the step size is set to Automatic, let it be a 5th of the biggest element of H at time 0, otherwise, the user has given the stepsize.*)


GetStepSize[H_,stepsize_:Automatic]:=
	If[stepsize==Automatic,
		1/(5*Max[Flatten[H[0]]]),
		stepsize
	]


Else=True;


End[];


(* ::Section:: *)
(*Simulators*)


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Decclare a private symbol. We will set this symbol to be time dependent Hamiltonians so that we can use them inside of compiled functions.*)


Hfun;
Hdim;


(* ::Text:: *)
(*Takes a time independent drift Hamiltonian, the control Hamiltonians, the time steps, and the amplitudes and multiplies all of the matrix exponentials together, just using the usual Fold method.*)


DotExpNaive=Function[{Hint,Hctls,dt,amps},
	Fold[
		(MatrixExp[-I*(#2[[1]])*(Hint+Total[Hctls*(#2[[{2,-1}]])])].#1)&,
		IdentityMatrix[Length[Hint]],
		Join[{dt},amps\[Transpose]]\[Transpose]
	]
]


(* ::Text:: *)
(*Takes a time independent drift Hamiltonian, the control Hamiltonians, the time steps, and the amplitudes and multiplies all of the matrix exponentials together, using a compiled function.*)


DotExpCompiled=Compile[{{Hint,_Complex,2},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},
	With[{num=Length[dt],dim=Length[Hint]},
		Block[{out=(1.0+0.0 I)*IdentityMatrix[dim],H=Table[0.0 I,{i,dim},{j,dim}]},
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


DotExpCompiledC:=Compile[{{Hint,_Complex,2},{Hctls,_Complex,3},{dt,_Real,1},{amps,_Real,2}},
	With[{num=Length[dt],dim=Length[Hint]},
		Block[{out=(1.0+0.0 I)*IdentityMatrix[dim],H=Table[0.0 I,{i,dim},{j,dim}]},
			For[k=1,k<=num,k++,
				H=-I*dt[[k]]*(Hint+Total[Hctls*amps[[k,All]]]);
				out=MatrixExp[H].out;
			];
			out
		]
	],
	CompilationTarget->"C"
];


(* ::Text:: *)
(*Compute the propagator of the time dependent Hamiltonian Hfun for a total time t and stepsize dt. Use the midpoint method.*)


TimeDepDriftCompiled=Compile[{{T,_Real},{dt,_Real}},
	Block[{N,out,H},
		N = Floor[T/dt];
		H=Table[0.0 I,{i,Hdim},{j,Hdim}];
		out = (1.0+0.0 I)*IdentityMatrix[Hdim];
		(*For[k=1,k<=N,k++,
			H = -I*dt*Hfun[(k-0.5)*dt];
			(*out = MatrixExp[H].out;      evolve for time dt *)
		];
		(*out = MatrixExp[-I*(T-N*dt)*Hfun[0.5*(T+N*dt)]].out;  add the extra bit on *)
		out*)
	]
];


End[];


(* ::Section:: *)
(*Single Pulse Evaluator*)


(* ::Subsection:: *)
(*Usage Declarations*)


EvalPulse::usage = "";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


EvalPulse[H_?DriftHamConstQ,p_?ShapedPulseQ,OptionsPattern[PulseOptions]]:=
	Block[{st,sm,simFunction,pulse,dt,Hctls,amps,Uj},
		st=OptionValue[SimTrace];
		sm=OptionValue[SimMethod];
		pulse=GetPulseShapeMatrix[p[[1]]];
		dt = pulse[[All,1]];
		amps = pulse[[All,{2,-1}]];
		Hctls = p[[2]];
		Which[
			st==False,
				simFunction=sm/.{Automatic->DotExpCompiled,Compiled->DotExpCompiled,CompiledC->DotExpCompiledC,Naive->DotExpNaive};
				Apply[simFunction,{H,Hctls,dt,amps}],
			st==True,
				Uj=Table[MatrixExp[-I*dt[[k]]*(H+Total[Hctls*amps[[k]]])],{k,Length[dt]}];
				FoldList[ #2.#1&,IdentityMatrix[Dimensions[H]],Uj]
		]
	]


EvalPulse[H_?DriftHamNonConstQ,p_?ShapedPulseQ,OptionsPattern[PulseOptions]]:=
	Block[{st,sm,simFunction,pulse,dt,Hctls,amps},
		st=OptionValue[SimTrace];
		sm=OptionValue[SimMethod];
		pulse=GetPulseShapeMatrix[p[[1]]];
		dt = pulse[[All,1]];
		amps = pulse[[All,{2,-1}]];
		Hctls = p[[2]];
		{}
	]


EvalPulse[H_?DriftHamQ,p_?InstantaneousPulseQ,OptionsPattern[PulseOptions]]:=
	With[{st=OptionValue[SimTrace]}, 
		Which[
			st==False,
				p,
			st==True,
				{p}
		]
	]


EvalPulse[H_?DriftHamConstQ,p_?DriftPulseQ,OptionsPattern[PulseOptions]]:=
	With[{st=OptionValue[SimTrace]}, 
		Which[
			st==False,
				MatrixExp[-I*p*H],
			st==True,
				{MatrixExp[-I*p*H]}
		]
	]


EvalPulse[H_?DriftHamNonConstQ,p_?DriftPulseQ,OptionsPattern[PulseOptions]]:=
	Block[{st,sm,simFunction,d\[Tau],T},
		st=OptionValue[SimTrace];
		sm=OptionValue[SimMethod];
		d\[Tau]=GetStepSize[H,OptionValue[StepSize]];
		T = p;
		Hfun = H;
		Hdim = Length[Hfun[0]];
		Which[
			st==False,
				simFunction=sm/.{Automatic->TimeDepDriftCompiled,Compiled->TimeDepDriftCompiled,CompiledC->NotImplemented,Naive->NotImplemented};
				Apply[simFunction, {T,d\[Tau]}],
			st==True,
				{}
		]
	]


End[];


(* ::Section::Closed:: *)
(*Pulse Sequence Evaluator*)


(* ::Subsection:: *)
(*UsageDeclarations*)


EvalPulseSequence::usage = "";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


EvalPulseSequence[H_?DriftHamQ,seq_?PulseSequenceQ,OptionsPattern[PulseOptions]]:=
	With[{st=OptionValue[SimTrace],sm=OptionsValue[SimMethod]},
		Which[
			st==False,
				Dot@@Reverse[EvalPulse[H,#,SimTrace->st,SimMethod->sm]&/@seq],
			st==True,
				Flatten[EvalPulse[H,#,SimTrace->st,SimMethod->sm]&/@seq,1]
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
	Block[{minx,maxx,n,d,data,totaltime},
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
	Block[{shapedFrac=0.3,instFrac=0.07,driftFrac=0.63,width=500,height=100,widths},
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


EndPackage[];
