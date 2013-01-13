(* ::Package:: *)

BeginPackage["NVSim`", {"QuantumUtils`","Predicates`","NVHamiltonian`"}]


(* ::Section:: *)
(*Predicates*)


(* ::Subsection::Closed:: *)
(*Usage Declarations*)


PulseShapeFileQ::usage = "PulseShapeFileQ[str] returns True iff str is a string pointing to a text file containg a pulse.";
PulseShapeMatrixQ::usage = "PulseShapeMatrixQ[M] returns True iff M is a 2D matrix";
PulseShapeQ::usage = "PulseShapeQ[in] returns True iff one of PulseShapeFileQ or PulseShapeMatrixQ is True";


ShapedPulseQ::usage = "ShapedPulse[p] returns True iff p is of the form {pulse,{Hcontrol1,Hcontrol2,..}} where the Hcontrols are the control Hamiltonians and pulse satisfies PulseShapeQ.";
DriftPulseQ::usage = "DriftPulseQ[p] returns True iff p is a real number indicating the amount of time to evolve under the drift Hamiltonian.";
InstantaneousPulseQ::usage = "InstantaneousPulseQ[p] returns True iff p is of the form {U,t} where U is a matrix and t is how \"long\" you want the instantaneous pulse to take.";
PulseQ::usage = "PulseQ[p] returns True iff p satisfies at least one of the NVSim predicates ending in PulseQ (ShapedPulseQ,InstantaneousPulseQ,etc)";


PulseSequenceQ::usage = "PulseSequenceQ[seq] returns True iff seq is a list where each element satisfies PulseQ.";


DriftHamConstQ::usage = "DriftHamConstQ[H] returns True iff H is a square matrix.";
DriftHamNonConstQ::usage = "DriftHamNonConst[H] returns True iff H is a function which accepts a number and returns a square matrix.";
DriftHamQ::usage = "DriftHamQ[H] returns True iff one of DriftHamConstQ[H] or DriftHamNonConstQ[H] is True.";


DensityMatrixQ::usage = "DensityMatrixQ[\[Rho]] returns True iff \[Rho] is a square matrix.";
ObservableListQ::usage = "ObservableListQ[obs] returns True iff obs is a list of square matrices.";
FunctionListQ::usage = "FunctionListQ[lst] retruns True iff lst is a List.";


Else=True;


(* ::Subsection::Closed:: *)
(*Implementations*)


Begin["`Private`"];


PulseShapeFileQ[str_]:=StringQ[str]


PulseShapeMatrixQ[M_]:=MatrixQ[M]


PulseShapeQ[in_]:=PulseShapeFileQ[in]||PulseShapeMatrixQ[in]


ShapedPulseQ[p_]:=ListQ[p]&&(Length[p]==2)&&ListQ[p[[2]]]&&PulseShapeQ[p[[1]]]


DriftPulseQ[p_]:=NumericQ[p]&&p\[Element]Reals


InstantaneousPulseQ[p_]:=ListQ[p]&&SquareMatrixQ[p[[1]]]&&(p[[2]]\[Element]Reals)


PulseQ[p_]:=Or@@(Through[{ShapedPulseQ,DriftPulseQ,InstantaneousPulseQ}[p]])


PulseSequenceQ[seq_]:=ListQ[seq]&&(And@@(PulseQ/@seq))


DriftHamConstQ[H_]:=SquareMatrixQ[H]


DriftHamNonConstQ[H_]:=SquareMatrixQ[H[0.0]]


DriftHamQ[H_]:=DriftHamConstQ[H]||DriftHamNonConstQ[H]


ObservableListQ[obs_]:=ListQ[obs]&&(And@@(SquareMatrixQ/@obs))


DensityMatrixQ[\[Rho]_]:=SquareMatrixQ[\[Rho]]


FunctionListQ[lst_]:=ListQ[lst]


End[];


(* ::Section:: *)
(*Options and Helper Functions*)


(* ::Subsection:: *)
(*Usage Declarations*)


SimulationOptions::usage = "SimulationOptions is a dummy function which stores the options for various NVSim functions. Use the command Options[SimulationOptions] to view these options and their default values.";


StepSize::usage = "StepSize is a simulation option that chooses the time discretization when the internal Hamiltonian is time dependent. Can be set to Automatic.";
PollingInterval::usage = "PollingInterval is a simulation option that specifies the time interval at which results of the simulation should be returned. The default value is Off.";
InitialState::usage = "InitialState is a simulation option. Set this option to the initial density matrix of your system. The default value is None.";
SimulationOutput::usage = "SimulationOutput is a simulation option. Set this option to be one of, or a nonempty subset of, the list {Unitaries,States,Observables,Functions}. These will be the values (and order) of the simulation output. The default value is Automatic.";
SequenceMode::usage = "SequenceMode is a simulation option. If set to True, EvalPulse returns the final state (or None) in addition to the usual output.";


TimeVector::usage = "TimeVector is a function Head used in a simulation's output. TimeVector[data] can also be used to extract the TimeVector from data.";
Unitaries::usage = "Unitaries is two things: (1) an element of the SimulationOutput list, and (2) a function Unitaries[data] (or Unitaries[data,t]) which exctracts the unitaries out of data (or extracts the unitary which happens at the closest time to t calculated), where data is in the form as outputed by EvalPulse.";
States::usage = "States is two things: (1) an element of the SimulationOutput list, and (2) a function States[data] (or States[data,t]) which exctracts the states out of data (or extracts the state which happens at the closest time to t calculated), where data is in the form as outputed by EvalPulse.";
Functions::usage = "Functions is three things: (1) simulation option which can be set to a list of functions which take a square matrix as input, (2) an element of the SimulationOutput list, and (3) a function Functions[data] (Functions[data,n]) which extracts all function values (n'th function values) from data, where data is in the format that EvalPulse outputs.";
Observables::usage = "Observables is three things: (1) simulation option which can be set to a list of observables (list of hermitian matrices), (2) an element of the SimulationOutput list, and (3) a function Observables[data] (Observables[data,n]) which extracts all observable values (n'th observable values) from data, where data is in the format that EvalPulse outputs.";


GetPulseShapeMatrix::usage = "GetPulseShapeMatrix[in] returns in if in is a matrix, but if in is a file name, returns the contents of that file as a matrix.";
GetStepSize::usage = "GetStepSize[H,stepsize:Automatic] returns a fifth of the biggest element of H[0] if stepsize is Automatic, and stepsize otherwise.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*Options and Input Handling*)


Options[SimulationOptions]={
	StepSize->Automatic,
	PollingInterval->Off,
	InitialState->None,
	Observables->None,
	Functions->None,
	SimulationOutput->Automatic,
	SequenceMode->False
};


GetPulseShapeMatrix[in_?PulseShapeFileQ]:=With[{out=Import[in]//N},Pick[out,Length[#]>1&/@out]]
GetPulseShapeMatrix[in_?PulseShapeMatrixQ]:=in//N


(* ::Text:: *)
(*If the step size is set to Automatic, let it be a tenth of the biggest element of H maximized over time, otherwise, the user has given the stepsize.*)


GetStepSize[H_,stepsize_:Automatic]:=
	If[stepsize===Automatic,
		Module[{t,max},
			max=Sqrt[NMaximize[Max[Abs[H[t]]^2],t][[1]]];
			1/(10*max)
		],
		stepsize
	]


FormatOutputAndReturn::usage="FormatOutputAndReturn[] returns all privately stored simulation data.";


(* ::Subsubsection:: *)
(*Private Variables*)


staVar::usage = "A private variable to store the initial state.";
obsVar::usage = "A private variable to store the observables.";
funVar::usage = "A private variable to store the functions.";


uniVals::usage = "A private variable to store the unitaries that will be returned.";
staVals::usage = "A private variable to store the states that will be returned.";
obsVals::usage = "A private variable to store the observable values that will be returned.";
funVals::usage = "A private variable to store the function values that will be returned.";
timVals::usage = "A private variable to store the polling times.";


returnKey::usage = "A private variable to encode which things our simulation is returning.";
outputList::usage = "A private variable to store which things our simulation is returning.";


AppendReturnables::usage = "A private function for appending new data to the already collected data.";


InitializePrivateVariables[OptionsPattern[SimulationOptions]]:=
	(
		(* initialize output storage containers*)
		uniVals={};
		staVals={};
		obsVals={};
		funVals={};
		timVals={};

		(* initialize the initial state, obserables, and monitoring functions *)
		staVar=OptionValue[InitialState];
		obsVar=OptionValue[Observables];
		funVar=OptionValue[Functions];

		(* initialize the outputList *)
		outputList=FormatSimulationOutput[OptionsPattern[SimulationOptions]];
		If[OptionValue[SimulationOutput]===Automatic,
			If[DensityMatrixQ[staVar],
				outputList={};
				If[ObservableListQ[obsVar],AppendTo[outputList,Observables]];
				If[FunctionListQ[funVar],AppendTo[outputList,Functions]];
				If[Length[outputList]===0,outputList={States}];,
				outputList={Unitaries};
			];,
			outputList=OptionValue[SimulationOutput];
			If[Not[ListQ[outputList]],outputList={outputList}];
			If[Not[DensityMatrixQ[staVar]]&&(MemberQ[outputList,States]||MemberQ[outputList,Observables]||MemberQ[outputList,Functions]),
				Print["Warning: You asked for an output which requires an InitialState, but no InitialState was specified."]
			];
		];
		AppendTo[outputList,TimeVector];

		(* initialize the returnKey *)
		returnKey=0;
		If[MemberQ[outputList,Unitaries],returnKey+=2^0];
		If[MemberQ[outputList,States],returnKey+=2^1];
		If[MemberQ[outputList,Observables],returnKey+=2^2];
		If[MemberQ[outputList,Functions],returnKey+=2^3];

		(* for each returnKey we need a different AppendReturnables *)
		Which[
			returnKey==1,AppendReturnables[U_,t_]:=(AppendTo[timVals,t];AppendTo[uniVals,U];),
			returnKey==2,AppendReturnables[U_,t_]:=(AppendTo[timVals,t];AppendTo[staVals,U.staVar.U\[ConjugateTranspose]];),
			returnKey==3,AppendReturnables[U_,t_]:=(AppendTo[timVals,t];AppendTo[uniVals,U];AppendTo[staVals,U.staVar.U\[ConjugateTranspose]];),
			returnKey==4,AppendReturnables[U_,t_]:=(AppendTo[timVals,t];AppendTo[obsVals,Re[Tr[#.U.staVar.U\[ConjugateTranspose]]]&/@obsVar];),
			returnKey==5,AppendReturnables[U_,t_]:=(AppendTo[timVals,t];AppendTo[uniVals,U];AppendTo[obsVals,Re[Tr[#.U.staVar.U\[ConjugateTranspose]]]&/@obsVar];),
			returnKey==6,AppendReturnables[U_,t_]:=With[{\[Rho]=U.staVar.U\[ConjugateTranspose]},AppendTo[timVals,t];AppendTo[staVals,\[Rho]];AppendTo[obsVals,Re[Tr[#.\[Rho]]]&/@obsVar];],
			returnKey==7,AppendReturnables[U_,t_]:=With[{\[Rho]=U.staVar.U\[ConjugateTranspose]},AppendTo[timVals,t];AppendTo[uniVals,U];AppendTo[staVals,\[Rho]];AppendTo[obsVals,Re[Tr[#.\[Rho]]]&/@obsVar];],
			returnKey==8,AppendReturnables[U_,t_]:=(AppendTo[timVals,t];AppendTo[funVals,(#[U.staVar.U\[ConjugateTranspose]])&/@funVar];),
			returnKey==9,AppendReturnables[U_,t_]:=(AppendTo[timVals,t];AppendTo[uniVals,U];AppendTo[funVals,(#[U.staVar.U\[ConjugateTranspose]])&/@funVar];),
			returnKey==10,AppendReturnables[U_,t_]:=With[{\[Rho]=U.staVar.U\[ConjugateTranspose]},AppendTo[timVals,t];AppendTo[staVals,\[Rho]];AppendTo[funVals,(#[\[Rho]])&/@funVar];],
			returnKey==11,AppendReturnables[U_,t_]:=With[{\[Rho]=U.staVar.U\[ConjugateTranspose]},AppendTo[timVals,t];AppendTo[uniVals,U];AppendTo[staVals,\[Rho]];AppendTo[funVals,(#[\[Rho]])&/@funVar];],
			returnKey==12,AppendReturnables[U_,t_]:=With[{\[Rho]=U.staVar.U\[ConjugateTranspose]},AppendTo[timVals,t];AppendTo[obsVals,Re[Tr[#.\[Rho]]]&/@obsVar];AppendTo[funVals,(#[\[Rho]])&/@funVar];],
			returnKey==13,AppendReturnables[U_,t_]:=With[{\[Rho]=U.staVar.U\[ConjugateTranspose]},AppendTo[timVals,t];AppendTo[uniVals,U];AppendTo[obsVals,Re[Tr[#.\[Rho]]]&/@obsVar];AppendTo[funVals,(#[\[Rho]])&/@funVar];],
			returnKey==14,AppendReturnables[U_,t_]:=With[{\[Rho]=U.staVar.U\[ConjugateTranspose]},AppendTo[timVals,t];AppendTo[staVals,\[Rho]];AppendTo[obsVals,Re[Tr[#.\[Rho]]]&/@obsVar];AppendTo[funVals,(#[\[Rho]])&/@funVar];],
			returnKey==15,AppendReturnables[U_,t_]:=With[{\[Rho]=U.staVar.U\[ConjugateTranspose]},AppendTo[timVals,t];AppendTo[uniVals,U];AppendTo[staVals,\[Rho]];AppendTo[obsVals,Re[Tr[#.\[Rho]]]&/@obsVar];AppendTo[funVals,(#[\[Rho]])&/@funVar];],
			Else,AppendReturnables[U_,t_]:=Null
		];
	)


(* ::Subsubsection:: *)
(*Output Formatting*)


FormatOutputAndReturn[]:=
	(
		outputList/.
			{Unitaries->{Unitaries,uniVals},
			States->{States,staVals},
			Observables->{Observables,obsVals},
			Functions->{Functions,funVals},
			TimeVector->{TimeVector,timVals}}
	)


TimeVector[data_]:=Select[data,(#[[1]]===TimeVector)&,1][[1,2]]


Unitaries[data_]:=Select[data,(#[[1]]===Unitaries)&,1][[1,2]]
Unitaries[data_,t_]:=With[{minpos=Ordering[Abs[t-#]&/@TimeVector[data],1][[1]]},Unitaries[data][[minpos]]]


States[data_]:=Select[data,(#[[1]]===States)&,1][[1,2]]
States[data_,t_]:=With[{minpos=Ordering[Abs[t-#]&/@TimeVector[data],1][[1]]},States[data][[minpos]]]


Observables[data_,OptionsPattern[{TimeVector->False}]]:=
	With[{obs=Transpose[Select[data,(#[[1]]===Observables)&,1][[1,2]]]},
		If[OptionValue[TimeVector]&&Length[obs]>0,
			With[{tv=TimeVector[data]},{tv,#}\[Transpose]&/@obs],
			obs
		]
	]
Observables[data_,n_,opt:OptionsPattern[{TimeVector->False}]]:=Observables[data,opt][[n]]


Functions[data_,OptionsPattern[{TimeVector->False}]]:=
	With[{obs=Transpose[Select[data,(#[[1]]===Functions)&,1][[1,2]]]},
		If[OptionValue[TimeVector]&&Length[obs]>0,
			With[{tv=TimeVector[data]},{tv,#}\[Transpose]&/@obs],
			obs
		]
	]
Functions[data_,n_,opt:OptionsPattern[{TimeVector->False}]]:=Functions[data,opt][[n]]


Minimize


(* ::Subsubsection::Closed:: *)
(*Epilog*)


Protect[PollingInterval,StepSize,IntitialState,Observables,Functions,SimulationOutput,Unitaries,States,TimeVector];


End[];


(* ::Section:: *)
(*Single Pulse Evaluator*)


(* ::Subsection:: *)
(*Usage Declarations*)


EvalPulse::usage = "EvalPulse[H,p] is the work house of the simulator. H is the Hamiltonion, either a matrix or a function accepting one real argument and returning a matrix, and p is a pulse. p must satisfy one of ShapedPulseQ, InstantaneousPulseQ, or DriftPulseQ. Type Options[SimulationOptions] to see all of the possible options.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*Shaped Pulse Evaluators*)


EvalPulse[H_?DriftHamConstQ,p_?ShapedPulseQ,opts:OptionsPattern[SimulationOptions]]:=
	Module[{dt,ds,\[Epsilon],upper,lower,pt,t,U,W,V,dim,n,pulse,amps,Hctls},
		InitializePrivateVariables[opts];
		pulse=GetPulseShapeMatrix[p[[1]]];
		dt = pulse[[All,1]];
		n=Length[dt];
		dim=Length[H];
		amps = If[Length[pulse[[1]]]>2,pulse[[All,{2,-1}]],pulse[[All,{2}]]];
		Hctls = p[[2]];
		pt=OptionValue[PollingInterval]//N;
		
		If[pt===Off,
			AppendReturnables[U=IdentityMatrix[dim],0];
			Table[
				U=MatrixExp[-I*dt[[k]]*(H+Total[Hctls*amps[[k]]])].U;,
				{k,n}
			];
			AppendReturnables[U,Total[dt]];,

			\[Epsilon]=10*$MachineEpsilon;
			AppendReturnables[U=IdentityMatrix[dim],0];
			t=0;
			Table[
				If[(lower=Ceiling[t/pt+\[Epsilon]])<=(upper=Floor[(t+dt[[k]])/pt+\[Epsilon]]),
					V=U;
					Table[
						If[(m==lower)||(m==upper)||(m==lower+1),
							ds=Min[t+dt[[k]],m pt]-Max[t,(m-1)pt];
							W=MatrixExp[-I*ds*(H+Total[Hctls*amps[[k]]])];
						];
						If[Abs[t+dt[[k]]-m*pt]>\[Epsilon],
							AppendReturnables[V=W.V,m*pt];
						];,
						{m,lower,upper}
					];
				];
				U=MatrixExp[-I*dt[[k]]*(H+Total[Hctls*amps[[k]]])].U;
				t=t+dt[[k]];,
				{k,n}
			];

			AppendReturnables[U,t]
		];
		If[OptionValue[SequenceMode],
			With[{\[Rho]=OptionValue[InitialState]},{If[\[Rho]===None,None,U.\[Rho].U\[ConjugateTranspose]],FormatOutputAndReturn[]}],
			FormatOutputAndReturn[]
		]
	]


EvalPulse[H_?DriftHamNonConstQ,p_?ShapedPulseQ,OptionsPattern[SimulationOptions]]:=
	Module[{simFunction,pulse,dt,d\[Tau],Hctls,amps},
		pulse=GetPulseShapeMatrix[p[[1]]];
		dt = pulse[[All,1]];
		amps = If[Length[pulse[[1]]]>2,pulse[[All,{2,-1}]],pulse[[All,{2}]]];
		Hctls = p[[2]];
		NVSim`Private`Hfun = H;
		d\[Tau]=GetStepSize[H,OptionValue[StepSize]];
		Which[
			OptionValue[PollingInterval]==False,
				simFunction=OptionValue[SimMethod]/.{Automatic->DotExpTimeDepCompiled,Compiled->DotExpTimeDepCompiled,CompiledC->DotExpTimeDepCompiledC,Naive->DotExpTimeDepNaive};
				Apply[simFunction,{d\[Tau],Hctls,dt,amps}],
			OptionValue[PollingInterval]==True,
				{},
			DensityMatrixQ[OptionValue[PollingInterval]],
				Print["This SimTrace option has not yet been implemented."];Abort[]
		]
	]


(* ::Subsubsection::Closed:: *)
(*Instantaneous Pulse Evaluators*)


EvalPulse[H_?DriftHamQ,p_?InstantaneousPulseQ,opts:OptionsPattern[SimulationOptions]]:=
	(
		InitializePrivateVariables[opts];
		AppendReturnables[IdentityMatrix[Length[p[[1]]]],0];
		AppendReturnables[p[[1]],p[[2]]];
		If[OptionValue[SequenceMode],
			With[{\[Rho]=OptionValue[InitialState]},{If[\[Rho]===None,None,p[[1]].\[Rho].p[[1]]\[ConjugateTranspose]],FormatOutputAndReturn[]}],
			FormatOutputAndReturn[]
		]
	)


(* ::Subsubsection::Closed:: *)
(*Drift Pulse Evaluators*)


EvalPulse[H_?DriftHamConstQ,T_?DriftPulseQ,opts:OptionsPattern[SimulationOptions]]:=
	Module[{dt,U,W},
		InitializePrivateVariables[opts];
		dt=OptionValue[PollingInterval]//N;
		Which[
			dt===Off||dt>T,
				AppendReturnables[IdentityMatrix[Length[H]],0];
				AppendReturnables[W=MatrixExp[-I*T*H],T];,
			Else,
				W=IdentityMatrix[Length[H]];
				U=MatrixExp[-I*dt*H];
				AppendReturnables[W,0];
				Table[
					W=W.U;
					AppendReturnables[W,k*dt];,
					{k,Floor[T/dt]}
				];
				dt=T-Floor[T/dt]*dt;
				AppendReturnables[W=W.MatrixExp[-I*dt*H],T];
		];
		If[OptionValue[SequenceMode],
			With[{\[Rho]=OptionValue[InitialState]},{If[\[Rho]===None,None,W.\[Rho].W\[ConjugateTranspose]],FormatOutputAndReturn[]}],
			FormatOutputAndReturn[]
		]
	]


EvalPulse[H_?DriftHamNonConstQ,T_?DriftPulseQ,opts:OptionsPattern[SimulationOptions]]:=
	Module[{dt,pt,U,dim,n},
		InitializePrivateVariables[opts];
		dt=GetStepSize[H,OptionValue[StepSize]];
		pt=OptionValue[PollingInterval]//N;
		dim=Length[H[0]];
		If[pt===Off,
			U=IdentityMatrix[dim];
			AppendReturnables[U,0];
			Table[U=MatrixExp[-I*dt*H[t-dt/2]].U;,{t,dt,T,dt}];
			U=MatrixExp[-I*(T-Floor[T/dt]*dt)*H[T-Floor[T/dt]*dt/2]].U;
			AppendReturnables[U,T];,

			If[pt<=dt,
				dt=pt;
				U=IdentityMatrix[dim];
				AppendReturnables[U,0];
				Table[AppendReturnables[U=MatrixExp[-I*dt*H[t-dt/2]].U,t],{t,dt,T,dt}];
				U=MatrixExp[-I*(T-Floor[T/dt]*dt)*H[T-Floor[T/dt]*dt/2]].U;
				AppendReturnables[U,T];,
				
				(* In this last case, we change dt to divide pt evenly without increasing its size *)
				dt=pt/Ceiling[pt/dt];
				U=IdentityMatrix[dim];
				AppendReturnables[U,0];
				Table[
					Table[U=MatrixExp[-I*dt*H[s-dt/2]].U;,{s,dt,pt,dt}];
					AppendReturnables[U,t];,
					{t,pt,T,pt}
				];
				(* We need to tack on two things now *)
				Table[U=MatrixExp[-I*dt*H[s-dt/2]].U;,{s,dt,T-Floor[T/pt]*pt,dt}];
				U=MatrixExp[-I*(T-Floor[T/dt]*dt)*H[T-Floor[T/dt]*dt/2]].U;
				AppendReturnables[U,T];
			];
		];
		If[OptionValue[SequenceMode],
			With[{\[Rho]=OptionValue[InitialState]},{If[\[Rho]===None,None,U.\[Rho].U\[ConjugateTranspose]],FormatOutputAndReturn[]}],
			FormatOutputAndReturn[]
		]
	]


(* ::Subsubsection::Closed:: *)
(*Epilog*)


End[];


(* ::Section::Closed:: *)
(*Pulse Sequence Evaluator*)


(* ::Subsection:: *)
(*UsageDeclarations*)


EvalPulseSequence::usage = "EvalPulseSequence[H,{p1,p2,p3,...}] evaluates the pulse sequence {p1,p2,p3,...} by evaluting each of EvalPulse[H,pi] where everything is properly tied together, ie., the initial state for one pulse is taken to be the final state of the previous pulse, etc. etc.";


(* ::Subsection:: *)
(*Implementations*)


Begin["`Private`"];


EvalPulseSequence[H_?DriftHamQ,seq_?PulseSequenceQ,options:OptionsPattern[SimulationOptions]]:=
	Module[{JoinTwoFields,JoinTwoEvalPulses},
		(* Define how to join each kind of output *)
		JoinTwoFields[{TimeVector,f1_},{TimeVector,f2_}]:={TimeVector,Join[f1,Last[f1]+Rest[f2]]};
		JoinTwoFields[{Unitaries,f1_},{Unitaries,f2_}]:={Unitaries,Join[f1,(#.Last[f1])&/@Rest[f2]]};
		JoinTwoFields[{type_,f1_},{type_,f2_}]:={type,Join[f1,Rest[f2]]};
		JoinTwoEvalPulses[p1_,p2_]:=MapThread[JoinTwoFields,{p1,p2},1];
		
		(* Now iteravely join each pulse. We need to deal with the slight complication of passing 
		   the previous final state to the new pulse as InitialState. This is only possible with 
			the SequenceMode option. *)
		Fold[
			With[{newEval=EvalPulse[H,#2,InitialState->#1[[1]],SequenceMode->True,options]},
				{newEval[[1]],JoinTwoEvalPulses[#1[[2]],newEval[[2]]]}
			]&,
			EvalPulse[H,First[seq],SequenceMode->True,options],
			Rest[seq]
		][[2]]
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
		Text[If[Length[p[[1]]]<3,p[[1]]//N,"Instant Pulse"],{width/2+offset,-height/6}]
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
