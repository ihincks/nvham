(* ::Package:: *)

BeginPackage["NVCalibrationTools`",{"NVHamiltonian`","NVSim`"}]


(* ::Section:: *)
(*Predicates*)


(* ::Section::Closed:: *)
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
(*Magnetic Fields*)


(* ::Subsection:: *)
(*Usage Declarations*)


Save2DField::usage = "Save2DField[filename,{rmin,rmax},{zmin,zmax},magmat] saves the magnetic field data stored in magmat to disk, where rmin, rmax, zmin, and zmax indicate the physical extent of the data in millimetres.";
Load2DField::usage = "Load2DField[filename,M:1000] is the inverse of Save2DField, with the output format {{rmin,rmax},{zmin,zmax},M*magmat}.";
Save3DField::usage = "Save3DField[filename,{xmin,xmax},{ymin,ymax},{zmin,zmax},magmat] saves the magnetic field data stored in magmat to disk, where xmin,xmax,ymin,ymax,zmin, and zmax indicate the physical extent of the data in millimetres.";
Load3DField::usage = "Load3DField[filename,M:1000] is the inverse of Save3DField, with the output format {{xmin,xmax},{ymin,ymax},{zmin,zmax},magmat}.";


DiskComputeFieldAtPoint::usage = "ComputeFieldPoint[{r,\[Phi],z},m:1000,d:7.5,h:2.5] computes the magnetic field at cylindrical coordinates {r,\[Phi],z} due to a solid cylindrical magnet at the origin of diameter d and height h, with a magnetization factor m. The return value is also in cylindrical coordinates. {0,0,0} is returned for locations within the magnet.";
ComputeDiskFieldAndSave::usage = "ComputeDiskFieldAndSave[filename,{rmin,rmax,dr},{zmin,zmax,dz},m:1000,d:7.5,h:2.5,parallel:False] computes the field using DiskComputeFieldAtPoint on each point of the rectangular array defined by {rmin,rmax,dr},{zmin,zmax,dz}. The results are stored to file using Save2DField. Calculation will be done in parallel if you set the parallel flag to True.";
Plot2DField::usage = "Plot2DField[fun2d] (OR Plot2DField[{fun2d1,fun2d2,...}] (doesn't work)) takes a InterpolatingFunction which is a 2D vector field and plots it. You probably want fun2d to be the output of Interpolate2DField.";


CubeComputeFieldAtPoint::usage = "CubeComputeFieldAtPoint[{x,y,z}, m:1000, l:5, w:5, h:5] computes the magnetic field due to a rectangular magnet located at the origin at the coordinate {x,y,z}. m is an overall scaling factor which you can think of as the density of ferromagnetism, or something physicsy like that.";
ComputeCubeFieldAndSave::usage = "ComputeCubeFieldAndSave[filename,{xmin,xmax,dx},{ymin,ymax,dy},{zmin,zmax,dz},m:1000,l:5,w:5,h:5,parallel:False] computes the field using CubeComputeFieldAtPoint on each point of the rectangular array defined by {xmin,xmax,dx},{ymin,ymax,dy},{zmin,zmax,dz}. The results are stored to file using Save3DField. Calculation will be done in parallel if you set the parallel flag to True.";
CubeMagnet::usage = "CubeMagnet[{x,y,z},fun3dup,fun3dright,orientation] returns the field at position {x,y,z} due to a rectangular magnet centered at the origin whose magnetic field is given by fun3dup or fun3dright (output of Interpolate3DField). This function has two purposes: (1) extend the domain of the field from the positive cone to all octants of space (this lets the data we save take up 1/8th of the space), and (2), flip the field if we want orientation 3 or 4, so we don't need to store redundant data.The magnetic field only needs to be computable in the positive cone; the symetry is taken into account. fun3dup is used for orientations 1 and 3, and fun3dright is used for orientations 2 and 4. Orientation should be one of the integers 1,2,3 or 4.";
HalbachTopField::usage = "HalbachTopField[{x,y,z},{topPos,bottomPos,azimuth},fun3dup,fun3dright,gap,m,l,w,h,topOrientations]";
HalbachBottomField::usage = "HalbachBottomField[{x,y,z},{topPos,bottomPos,azimuth},fun3dup,fun3dright,gap,m,l,w,h,bottomOrientations]";
HalbachField::usage = "HalbachField[{x,y,z},{topPos,bottomPos,azimuth},fun3dup,fun3dright,gap,m,l,w,h,bottomOrientations,topOrientations]";


Interpolate2DField::usage = "Interpolate2DField[data] accepts the output of Load2DField, and interpolates the matrix to third order, returning an InterpolationFunciton.";
CylindricalFieldExtrapolation::usage = "FieldInterpolation[{x,y,z},fun2d] returns the field at position {x,y,z} of a cylindrically symmetric, radial field defined by the 2D vector field fun2d, which is a slice of the field through the z-x plane. You almost certainly want to set fun2d to the output of Interpolate2DField.";
Interpolate3DField::usage = "Interpolate2DField[data] accepts the output of Load3DField, and interpolates the matrix to third order, returning an InterpolationFunciton.";


CompareSolutionFields::usage = "CompareSolutionFields[{z1,z2,z3},fun2d,{zz1,zz2,zz3},funn2d] takes the normalized inner product between solution translations {z1,z2,z3} and {zz1,zz2,zz3} over the possible volume swept by the magnet motors.";


(* ::Subsection:: *)
(*Implementations*)


(* ::Subsubsection::Closed:: *)
(*Loading and Saving*)


Begin["`Private`"];


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


Save3DField[filename_,xbounds_,ybounds_,zbounds_,magmat_]:=
	Block[{xmin,xmax,ymin,ymax,zmin,zmax,mat},
		xmin=xbounds[[1]];
		xmax=xbounds[[2]];
		ymin=ybounds[[1]];
		ymax=ybounds[[2]];
		zmin=zbounds[[1]];
		zmax=zbounds[[2]];
		mat=magmat;
		Save[filename,{xmin,xmax,ymin,ymax,zmin,zmax,mat}];
		Print["Save Complete."];
	]


Load3DField[filename_String,M_]:=
	Block[{xmin,xmax,ymin,ymax,zmin,zmax,mat},
		Get[filename];
		{{xmin,xmax},{ymin,ymax},{zmin,zmax},(M/1000)*mat}
	]


End[];


(* ::Subsubsection::Closed:: *)
(*Field of a Disk*)


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


End[];


(* ::Subsubsection::Closed:: *)
(*Field of a Halbach Array*)


Begin["`Private`"];


CubeComputeFieldAtPoint[{x_?NumberQ,y_?NumberQ,z_?NumberQ},m_:1000,l_:5,w_:5,h_:5]:=
	If[
		Abs[x]<=l/2&&Abs[y]<=w/2&&Abs[z]<=h/2, {0,0,0},
		Module[{fun,x0,y0,z0},
			fun[x0_,y0_,z0_]:=With[{r=Sqrt[(x0-x)^2+(y0-y)^2+(z0-z)^2]},
				(m/r^3)*(3*(z-z0)*{x-x0,y-y0,z-z0}/r^2 - {0,0,1})
			];
			NIntegrate[fun[x0,y0,z0],{x0,-l/2,l/2},{y0,-w/2,w/2},{z0,-h/2,h/2},AccuracyGoal->5]
		]
	]


ComputeCubeFieldAndSave[filename_,{xmin_,xmax_,dx_},{ymin_,ymax_,dy_},{zmin_,zmax_,dz_},m_:1000,l_:5,w_:5,h_:5,parallel_:False]:=
	Module[{data,tablefcn,x,y,z},
		tablefcn=If[parallel,ParallelTable[##,Method->"CoarsestGrained"]&,Table];
		data=tablefcn[
				CubeComputeFieldAtPoint[{x,y,z},m,l,w,h],
				{x,xmin,xmax,dx},
				{y,ymin,ymax,dy},
				{z,zmin,zmax,dz}
			];
		Save3DField[
			filename,
			{xmin,xmax},
			{ymin,ymax},
			{zmin,zmax},
			data
		]
	]


CubeMagnet[{x_,y_,z_},fun3dup_,fun3dright_,orientation_]:=
	Which[
		orientation===1,
			{Sign[x]*Sign[z],Sign[y]*Sign[z],1}*fun3dup[Abs[x],Abs[y],Abs[z]],
		orientation===2,
			{1,Sign[y]*Sign[x],Sign[z]*Sign[x]}*fun3dright[Abs[x],Abs[y],Abs[z]],
		orientation===3,
			-{Sign[x]*Sign[z],Sign[y]*Sign[z],1}*fun3dup[Abs[x],Abs[y],Abs[z]],
		orientation===4,
			-{1,Sign[y]*Sign[x],Sign[z]*Sign[x]}*fun3dright[Abs[x],Abs[y],Abs[z]]
	]


HalbachTopField[{x_,y_,z_},{topPos_,bottomPos_,azimuth_},fun3dup_,fun3dright_,gap_:1.4125,m_:1000,l_:6.35,w_:6.35,h_:6.35/4,orientations_:{2,1,4,3,2,1}]:=
	Module[{numMagnets,B},
		numMagnets=Length[orientations];
		(* put the top magnets 4*h+gap above the origin *)
		B=Sum[CubeMagnet[{x,y,z}-{l/2+(k-1)*l+topPos,0,4*h+gap+h/2},fun3dup,fun3dright,orientations[[k]]],{k,numMagnets}];
		m*(B)/1000
	]


HalbachBottomField[{x_,y_,z_},{topPos_,bottomPos_,azimuth_},fun3dup_,fun3dright_,gap_:1.4125,m_:1000,l_:6.35,w_:6.35,h_:6.35/4,orientations_:{3,2,1,4,3,2}]:=
	Module[{numMagnets,B},
		numMagnets=Length[orientations];
		(* stack four magnets vertically for each orientation *)
		B=Sum[CubeMagnet[{x,y,z}-{l/2+(k-1)*l+bottomPos,0,h/2+(j-1)*h},fun3dup,fun3dright,orientations[[k]]],{k,numMagnets},{j,4}];
		m*(B)/1000
	]


HalbachField[{x_,y_,z_},{topPos_,bottomPos_,azimuth_},fun3dup_,fun3dright_,gap_:1.4125,m_:1000,l_:6.35,w_:6.35,h_:6.35/4,bottomOrientations_:{3,2,1,4,3,2},topOrientations_:{2,1,4,3,2,1}]:=
	(
		HalbachTopField[{x,y,z},{topPos,bottomPos,azimuth},fun3dup,fun3dright,gap,m,l,w,h,topOrientations]+
		HalbachBottomField[{x,y,z},{topPos,bottomPos,azimuth},fun3dup,fun3dright,gap,m,l,w,h,bottomOrientations]
	)


End[];


(* ::Subsubsection:: *)
(*Interpolating Fields*)


Begin["`Private`"];


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
		out={Cos[\[Phi]]*val[[1]],Sin[\[Phi]]*val[[1]],val[[2]]};
		If[z>=0,out,out*{-1,-1,1}]
	]


Interpolate3DField[data_]:=
	Block[{xmin,xmax,ymin,ymax,zmin,zmax,mat,L,W,H},
		xmin = data[[1,1]];
		xmax = data[[1,2]];
		ymin = data[[2,1]];
		ymax = data[[2,2]];
		zmin = data[[3,1]];
		zmax = data[[3,2]];
		L = Dimensions[data[[4]]][[1]];
		W = Dimensions[data[[4]]][[2]];
		H = Dimensions[data[[4]]][[3]];
		Interpolation[
			Flatten[
				Table[
					{
						xmin+(xmax-xmin)(i-1)/(L-1), 
						ymin+(ymax-ymin)(j-1)/(W-1), 
						zmin+(zmax-zmin)(k-1)/(H-1), 
						data[[4,i,j,k]]
					},
					{i,H},{j,W},{k,H}
				],
			2]
		]
	]


End[];


(* ::Subsubsection::Closed:: *)
(*Plotting and Comparing (might be partially broken)*)


Begin["`Private`"];


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
NVZeemanSplitting::usage = "NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},nvOrientation] computes the Zeeman splitting of the NV electron in MHz given the magnetic field, the crystal orientation, and the nv orientation. Uses eigenvales; not the secular approximation.";
NVCenterFrequency::usage = "NVCenterFrequency[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},nvOrientation] computes the center frequency of the NV electron in MHz given the magnetic field, the crystal orientation, and the nv orientation. Uses eigenvales; not the secular approximation.";


FakeData::usage = "FakeData[\[Sigma],points,nvOrientation_Integer,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m]
generates fake data of the form {{x1,y1,z1,splitting1,nvorientation1},...} by sampling from a Gaussian of variance \[Sigma] at each point of the list points={{x1,y1,z1},{x2,y2,z3},...}, with the given frame transformation z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3, offset field {B1,B2,B3} and magnitization m.
FakeData[\[Sigma],points,nvOrientations_List,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m]
generates fake data of the form {{x1,y1,z1,splitting1,nvorientation1},...} by sampling from a Gaussian of variance \[Sigma] at each point of the list points={{x1,y1,z1},{x2,y2,z3},...} with respective nv orientations nvOrientations={nv1,nv2,...}, with the given frame transformation z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3, offset field {B1,B2,B3} and magnetization m.
FakeData[\[Sigma],points,nvOrientations_List,solution] calls FakeData[\[Sigma],points,nvOrientations_List,Sequence@@solution]
FakeData, in all cases, returns an output which is compatible with the data input of both PlotSplittingData and CalibrateField.";
PlotSplittingData::usage = "PlotSplittingData[data,minsplit:0] 
shows you the splitting at each point in space in your data set. 
All data points which have a splitting less than minsplit will be Black, otherwise, each NV orientation gets a different colour.";


MostInformativeSpot::usage = "MostInformativeSpot[{z11,z12,z13},{\[Alpha]11,\[Alpha]12,\[Alpha]13},{B11,B12,B13},m1,{z21,z22,z23},{\[Alpha]21,\[Alpha]22,\[Alpha]23},{B21,B22,B23},m2,nvOrientation,minsplitting:15,maxsplitting:280,minz:6]
computes where you should do your next measurement to maximally distinguish (which means to maximize the difference of the splittings) the two solutions given by the two z's, \[Alpha]'s, B's and m's. You can specifiy the maximum and minimum splitting you can tolerate, as well as the lowest z3 value
MostInformativeSpot[solution1,solution2,nvOrientation,minsplitting:15,maxsplitting:280,minz:6] calls MostInformativeSpot[Sequence@@solution1,Sequence@@solution2,nvOrientation,minsplitting:15,maxsplitting:280,minz:6]";
MagnetFittingCostFunction::usage = "MagnetFittingCostFunction[{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m,data]
runs the cost function used by CalibrateField on the potential solution {z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3}, where m is the magnetization factor and {B1,B2,B3} is the offset field, and data is of the same form as accepted by CalibrateField.
MagnetFittingCostFunction[solution,data] returns MagnetFittingCostFunction[Sequence@@solution,data]";
CalibrateField::usage = "CalibrateField[data,m,constraints:{{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]}},method:\"NelderMead\",minsplit:13]
 does a numerical minimization to find the best frame transformation between the magnet origin and the NV PAS given the data.
The data should be a list of the form {{x1,y1,z1,S1,nv1},{x2,y2,z2,S2,nv2},...} where the x,y,z are the positions of the magnet motor, the S's are the splittings in MHz, and the nv's are the NV orientations. m is the magenetization factor, and minsplit specifies which data points to throw out.
The output of the function is of the form {mincostfunctionvalue,{{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m}} where mincostfunctionvalue is in units \!\(\*SuperscriptBox[\(MHz\), \(2\)]\), {z1,z2,z3} is the translation of the stage, {\[Alpha]1,\[Alpha]2,\[Alpha]3} are the Euler angles, {B1,B2,B3} is a constant offset field in Gauss, and m is the magnetization (the same number you input).";


$typicalSolution::usage = "$typicalSolution is the sort of output you expect to get from CalibrateField. This is useful for feeding into functions like FakeData and PredictSplitting.";


PositionGivenField::usage = "PositionGivenField[{z1,z2,z3},{B1,B2,B3},m,Bvector,minz:6]
finds the position you should set the motors to in order to get the magenetic field Bvector (Bvector=Vector[{Bx,By,Bz},Coords]) in the lab frame. You can demand that the z-motor be greater than minz.
{z1,z2,z3}, {B1,B2,B3}, and m, are part of the calibration solution.
The output is in the form {costfunctionvalue,{x,y,z}}, where costfunctionvalue has units of MHz, {x,y,z} is where you should put the motors.
PositionGivenField[solution,B,BCoords:Cartesian,minz:6] returns PositionGivenField[solution[[1]],solution[[3]],solution[[4]],Bvector,minz:6].";
AlignField::usage = "AlignField[{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m,nvOrientation,absB,minz:6] invokes PositionGivenField to find the motor position at which the magnetic field is aligned with the given nvOrienation, and with field strength Babs in Gauss. WARNING: The sign of the field matters. For example, if you ask for an alignment with orientation 4, it will try to find a field which is parallel (and not anti-parallel) to that axis. You can input negative values of Babs if you want to search for anti-aligned solutions.
AlignField[solution,nvOrientation,absB,minz:6] returns AlignField[Sequence@@solution,nvOrientation,absB,minz:6]";


PredictSplitting::usage = "PredictSplitting[{r1,r2,r3},nvOrientation,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m] 
predicts what the splitting will be at the magnet position {r1,r2,r3} given the the solution specified by z, \[Alpha], B and m.
PredictSplitting[{r1,r2,r3},nvOrientation,solution] returns PredictSplitting[{r1,r2,r3},nvOrientation,Sequence@@solution].";
PredictCenterFreq::usage = "PredictCenterFreq[{r1,r2,r3},nvOrientation,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m] 
predicts what the center frequency will be at the magnet position {r1,r2,r3} given the the solution specified by z, \[Alpha], B and m.
PredictCenterFreq[{r1,r2,r3},nvOrientation,solution] returns PredictCenterFreq[{r1,r2,r3},nvOrientation,Sequence@@solution].";
PredictNVPASVector::usage = "PredictNVPASVector[{r1,r2,r3},nvOrientation,{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m] 
predicts what the magnetic field vector will be at the magnet position {r1,r2,r3} in the NV PAS given the the solution specified by z, \[Alpha], B and m.
PredictNVPASVector[{r1,r2,r3},nvOrientation,solution] returns PredictNVPASVector[{r1,r2,r3},nvOrientation,Sequence@@solution].";


(* ::Subsubsection:: *)
(*Implementations*)


Begin["`Private`"];


$typicalSolution = {{1.34,23.5,-5.3},{0,-0.977,-0.04},{0,0,0},.71};


ImportSplittingData[filename_String,scale_List:{1,1,1,1,1}]:=
	Catch[
		If[FindFile[filename]==$Failed,Throw["File "<>filename<>" not found."];];
		Block[{x=Import[filename]},
			If[Length[Dimensions[x]]>2,x=Flatten[x,1];];
			If[Dimensions[x][[2]]!=5,Throw["Your data doesn't appear to include NV orientations. Call ImportSplittingData[filename,nvOrientation] to specify all your orientations as equal."]];
			x=(#*scale)&/@x;
			ReplacePart[x\[Transpose],5->Round[x[[All,5]]]]\[Transpose]
		]
	]
ImportSplittingData[filename_String,nvOrientation_Integer,scale_List:{1,1,1,1,1}]:=
	Catch[
		If[FindFile[filename]==$Failed,Throw["File "<>filename<>" not found."];];
		Block[{x=Import[filename]},
			If[Length[Dimensions[x]]>2,x=Flatten[x,1];];
			x=Join[x[[All,{1,2,3,4}]]\[Transpose],{Table[nvOrientation,{k,Length[x]}]}]\[Transpose];
			(#*scale)&/@x
		]
	]


ExtractPositionArray[data_]:={#[[1]],#[[2]],#[[3]]}&/@data


CompareDataWithFourOrientations[data_,plotRange_:{0,25}]:=
	GraphicsRow[
		Join[
			{PlotSplittingData[data,0,plotRange]},
			PlotSplittingData[#,0,plotRange]&/@Table[FakeData[$MachineEpsilon,ExtractPositionArray[data],n,Sequence@@$typicalSolution],{n,4}]
		],
		ImageSize->1500
	]


(* ::Text:: *)
(*The following functions can be derived by going to the lab and looking at which ways the motors go. When facing the setup with your back to the computer, we want the lab frame to be the same as NVHamiltonian: z is up, x is to the right, and y is back towards the wall (right handed coordinates system).*)


MagToLab[{x_,y_,z_}]:={y,-x,-z};
LabToMag[{x_,y_,z_}]:={-y,x,-z};


(* ::Text:: *)
(*First get the Hamiltonians for each orientation in units of MHz, where the magnetic field is entered in the LabFrame. This is stored in MagHam.dat because the simplification it involves is lengthy.*)


Get[Evaluate@FileNameJoin[{"data","MagHam.dat"}]];


(* ::Text:: *)
(*Now the zeeman splitting is simply the difference of the two biggest eigenvalues of MagHam*)


NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},n_]:=Abs@First@Differences@Eigenvalues[MagHam[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},n],2]


(* ::Text:: *)
(*Define for orientations 5-8 too...*)


NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},5]:=NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},1]
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},6]:=NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},2]
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},7]:=NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},3]
NVZeemanSplitting[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},8]:=NVZeemanSplitting[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},4]


(* ::Text:: *)
(*Similarly figure out the centre frequencies:*)


NVCenterFrequency[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},n_]:=With[{e=Sort@Eigenvalues[MagHam[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},n]]},((e[[2]]-e[[1]])+(e[[3]]-e[[1]]))/2]
NVCenterFrequency[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},5]:=NVCenterFrequency[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},1]
NVCenterFrequency[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},6]:=NVCenterFrequency[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},2]
NVCenterFrequency[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},7]:=NVCenterFrequency[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},3]
NVCenterFrequency[{Bx_,By_,Bz_},{\[Theta]z1_,\[Theta]y_,\[Theta]z2_},8]:=NVCenterFrequency[{Bx,By,Bz},{\[Theta]z1,\[Theta]y,\[Theta]z2},4]


FakeData[\[Sigma]_,points_,nvOrientation_Integer,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},{B1_,B2_,B3_},m_]:=
	With[{x=#[[1]],y=#[[2]],z=#[[3]]},
		{
			x,y,z,
			RandomVariate[NormalDistribution[
				NVZeemanSplitting[Field[MagToLab[{z1,z2,z3}-{x,y,z}],m]+{B1,B2,B3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},nvOrientation],\[Sigma]]],
			nvOrientation
		}
	]&/@points
FakeData[\[Sigma]_,points_,nvOrientations_List,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},{B1_,B2_,B3_},m_]:=
	With[{x=#[[1]],y=#[[2]],z=#[[3]],nv=#[[4]]},
		{
			x,y,z,
			RandomVariate[NormalDistribution[
				NVZeemanSplitting[Field[MagToLab[{z1,z2,z3}-{x,y,z}],m]+{B1,B2,B3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},nv],\[Sigma]]],
			nv
		}
	]&/@(Join[points\[Transpose],{nvOrientations}]\[Transpose])
FakeData[\[Sigma]_,points_,nvOrientation_,solution_]:=FakeData[\[Sigma],points,nvOrientation,Sequence@@solution]


PlotSplittingData[data_,minsplit_:0,plotRange_:{0,25}]:=
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
			PlotRange->plotRange,
			PlotLabel->"Max Splitting: "<>ToString[max]<>"MHz",
			Epilog-> Inset[Framed[Column[Table[Style["Orientation "<>ToString[n], 12,colours[n]],{n,4}]], Background->LightYellow], {Right, Bottom}, {Right, Bottom}]]
	]


(* ::Text:: *)
(*Here is the cost function for the fitting. WARNING: it is duplicated in the CalibrationField for efficiency reasons. If you change one, be sure to change the other.*)


MagnetFittingCostFunction[{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},{B1_,B2_,B3_},m_,data_]:=
	Sum[
			(NVZeemanSplitting[
				Field[MagToLab[{z1,z2,z3}-data[[k,{1,2,3}]]],m]+{B1,B2,B3},
				{\[Alpha]1,\[Alpha]2,\[Alpha]3},
				data[[k,5]]
			]-data[[k,4]])^2
		,
		{k,1,Length[data]}
	]/Length[data];
MagnetFittingCostFunction[solution_,data_]:=MagnetFittingCostFunction[Sequence@@solution,data]


MostInformativeSpot[{z11_,z12_,z13_},{\[Alpha]11_,\[Alpha]12_,\[Alpha]13_},{B11_,B12_,B13_},m1_,{z21_,z22_,z23_},{\[Alpha]21_,\[Alpha]22_,\[Alpha]23_},{B21_,B22_,B23_},m2_,nvOrientation_,minsplitting_:15,maxsplitting_:280,minz_:6]:=
	Block[{r1,r2,r3,S1,S2,cost,result,values},
		S1[rr1_Real,rr2_Real,rr3_Real]:=NVZeemanSplitting[Field[MagToLab[{z11-rr1,z12-rr2,z13-rr3}],m1]+{B11,B12,B13},{\[Alpha]11,\[Alpha]12,\[Alpha]13},nvOrientation];
		S2[rr1_Real,rr2_Real,rr3_Real]:=NVZeemanSplitting[Field[MagToLab[{z21-rr1,z22-rr2,z23-rr3}],m2]+{B21,B22,B23},{\[Alpha]21,\[Alpha]22,\[Alpha]23},nvOrientation];
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


CalibrateField[data_,m_,constraints_:{{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]},{-\[Infinity],\[Infinity]}},method_:"NelderMead",minsplit_:13]:=
	Block[{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3,B1,B2,B3,cost,sdata,cons,result,values},
		sdata=Select[data,(#[[4]]>=minsplit)&];
		cost[zz1_,zz2_,zz3_,\[Alpha]\[Alpha]1_,\[Alpha]\[Alpha]2_,\[Alpha]\[Alpha]3_,BB1_,BB2_,BB3_]:=
			Sum[
				(
					NVZeemanSplitting[
						Field[MagToLab[{zz1,zz2,zz3}-sdata[[k,{1,2,3}]]],m]+{BB1,BB2,BB3},
						{\[Alpha]\[Alpha]1,\[Alpha]\[Alpha]2,\[Alpha]\[Alpha]3},sdata[[k,5]]]-sdata[[k,4]]
				)^2,
				{k,1,Length[sdata]}
			]/Length[sdata];
		{result,values}=NMinimize[
			{
				Hold[cost[z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3,B1,B2,B3]],
				constraints[[1,1]]<=z1<=constraints[[1,2]],
				constraints[[2,1]]<=z2<=constraints[[2,2]],
				constraints[[3,1]]<=z3<=constraints[[3,2]],
				constraints[[4,1]]<=\[Alpha]1<=constraints[[4,2]],
				constraints[[5,1]]<=\[Alpha]2<=constraints[[5,2]],
				constraints[[6,1]]<=\[Alpha]3<=constraints[[6,2]],
				constraints[[7,1]]<=B1<=constraints[[7,2]],
				constraints[[7,1]]<=B2<=constraints[[7,2]],
				constraints[[7,1]]<=B3<=constraints[[7,2]]
			},
			{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3,B1,B2,B3},
			Method->method,
			MaxIterations->3000
		];
		{z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3,B1,B2,B3}={z1,z2,z3,\[Alpha]1,\[Alpha]2,\[Alpha]3,B1,B2,B3}/.values;
		{result,{{z1,z2,z3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},{B1,B2,B3},m}}
	]


PositionGivenField[{z1_,z2_,z3_},{B1_,B2_,B3_},m_,Bvec_,minz_:6]:=
	Block[{cost,Bcart,result,values,r1,r2,r3},
		Bcart=Value@Cartesian@Bvec;
		cost[rr1_Real,rr2_Real,rr3_Real]:=Norm[Bcart-(Field[MagToLab[{z1,z2,z3}-{rr1,rr2,rr3}],m]+{B1,B2,B3})]^2;
		{result,values}=NMinimize[
			{
				cost[r1,r2,r3],
				0<=r1<=25,0<=r2<=25,minz<=r3<=25
			},
			{r1,r2,r3},
			MaxIterations->3000,
			Method->"NelderMead"
		];
		{Sqrt[result],{r1,r2,r3}/.values}
	]
PositionGivenField[solution_,B_,BCoords_:Cartesian,minz_:6]:=PositionGivenField[solution[[1]],solution[[3]],solution[[4]],B,BCoords,minz]


AlignField[{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},{B1_,B2_,B3_},m_,nvOrientation_,absB_,minz_:6]:=
	PositionGivenField[
		{z1,z2,z3},
		{B1,B2,B3},
		m,
		Vector[FrameChange[{0,0,absB},IdentityFrame,FrameCompose[EulerAngles[\[Alpha]1,\[Alpha]2,\[Alpha]3],NVOrientationToFrame[nvOrientation]]],Cartesian],
		minz
	];
AlignField[solution_,nvOrientation_,absB_,minz_:6]:=AlignField[Sequence@@solution,nvOrientation,absB,minz]


PredictSplitting[{r1_,r2_,r3_},nvOrientation_,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},{B1_,B2_,B3_},m_]:=
	NVZeemanSplitting[Field[MagToLab[{z1,z2,z3}-{r1,r2,r3}],m]+{B1,B2,B3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},nvOrientation]
PredictSplitting[{r1_,r2_,r3_},nvOrientation_,solution_]:=PredictSplitting[{r1,r2,r3},nvOrientation,Sequence@@solution]


PredictCenterFreq[{r1_,r2_,r3_},nvOrientation_,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},{B1_,B2_,B3_},m_]:=
	NVCenterFrequency[Field[MagToLab[{z1,z2,z3}-{r1,r2,r3}],m]+{B1,B2,B3},{\[Alpha]1,\[Alpha]2,\[Alpha]3},nvOrientation]
PredictCenterFreq[{r1_,r2_,r3_},nvOrientation_,solution_]:=PredictCenterFreq[{r1,r2,r3},nvOrientation,Sequence@@solution]


PredictNVPASVector[{r1_,r2_,r3_},nvOrientation_,{z1_,z2_,z3_},{\[Alpha]1_,\[Alpha]2_,\[Alpha]3_},{B1_,B2_,B3_},m_]:=
	Vector[FrameChange[Field[MagToLab[{z1,z2,z3}-{r1,r2,r3}],m]+{B1,B2,B3},IdentityFrame,FrameCompose[EulerAngles[\[Alpha]1,\[Alpha]2,\[Alpha]3],NVOrientationToFrame[nvOrientation]]],Cartesian]
PredictNVPASVector[{r1_,r2_,r3_},nvOrientation_,solution_]:=PredictNVPASVector[{r1,r2,r3},nvOrientation,Sequence@@solution]


End[];


(* ::Section:: *)
(*Carbon Hyperfine Estimation*)


(* ::Subsubsection:: *)
(*Usage Declarations*)


EnhancedLarmourFormula::usage = "EnhancedLarmourFormula[B_,\[Theta]_,\[Phi]_,\[Alpha]_,A1_,A2_,A3_] returns the enhanced Larmour precession frequency in Hz of a carbon coupled to an NV to second order, where {B,\[Theta],\[Phi]} describes the static magnetic field in the PAS, and the hyperfine tensor is {{A1,0,\[Alpha]},{0,A2,0},{\[Alpha],0,A3}} in the same frame."


EnhancedLarmourVisibility::usage = "EnhancedLarmourVisibility[staticField,\[Alpha]_,A1_,A2_,A3_,tp_,\[CapitalOmega]_] does a full simulation of an enhanced larmour experiment to return both the enhanced larmour frequency, and the visibility of the frequency.";


EnhancedLarmourSimulation::usage = "EnhancedLarmourSimulation[staticField,A,T,\[CapitalOmega],doplot:True]";


ELVHyperfineReconstruction::usage = "ELVHyperfineReconstruction[enhancedLarmourFreqs,staticFields,A\[Phi],method:\"NelderMead\"] Enhanced Larmour with Visibility";


EnhancedLarmourHyperfineReconstruction::usage = "EnhancedLarmourHyperfineReconstruction[enhancedLarmourFreqs,staticFields,A\[Phi],method:\"NelderMead\"]";
EnhancedLarmourHyperfineErrorBars::usage = "EnhancedLarmourHyperfineErrorBars[enhancedLarmourFreqs_,enhancedLarmourStds_,staticFields_,staticFieldStds_,A\[Phi]_,N_,opt:OptionsPattern[EnhancedLarmourHyperfineReconstruction]]";


(* ::Subsubsection:: *)
(*Implementations*)


Begin["`Private`"];


(* ::Text:: *)
(*Load the enhanced Larmour formula from file to save time. It can be calculated from scratch by calling:*)
(*	A={{A1,0,\[Alpha]},{0,A2,0},{\[Alpha],0,A3}}*)
(*	H=NVAverageHamiltonian[2,\[CapitalDelta],Carbon[A],StaticField->Vector[{B,\[Theta],\[Phi]},Spherical],Numerical->False];*)
(*	With[{HH=H[[3;;4,3;;4]]},Sqrt[Tr[H.X]^2+Tr[H.Y]^2+Tr[H.Z]^2]]*)
(*With subsequent value insertions and simplifications...*)


Get[Evaluate@FileNameJoin[{"data","EnhancedLarmourFormula.dat"}]];


(* ::Text:: *)
(*This private function takes the output of EvalPulse and tries to fit a cosine to the first observable.*)


FitCosine[data_]:=
	Module[{model,A,B,\[Omega],t,soln},
		model[t_]=A+B*Cos[2\[Pi] \[Omega] t];
		soln=FindFit[First@Observables[data,TimeVector->True],model[t],{A,B,\[Omega]},{t},Method->NMinimize];
		{model,soln}
	]


EnhancedLarmourVisibility[staticField_,\[Alpha]_,A1_,A2_,A3_,tp_,\[CapitalOmega]_]:=Module[{\[Rho],M,H,D,A,UP,dU,ndata=100,dt,T,data,\[Omega]elGuess,model,t,offset,amp,\[Omega],soln},
	\[Omega]elGuess=EnhancedLarmourFormula[Sequence@@Value@Spherical@staticField,\[Alpha],A1,A2,A3];
	T=4/\[Omega]elGuess;
	A=10^6*{{A1,0,\[Alpha]},{0,A2,0},{\[Alpha],0,A3}};
	D=10^6*NVCenterFrequency[Value@Cartesian[staticField],{0,0,0},1]-(10^6*NVZeemanSplitting[Value@Cartesian[staticField],{0,0,0},1]+Norm[A[[3]]])/2;
	H=NVAverageHamiltonian[2,D,Carbon[A],Felton09Nucleus["14N"],StaticField->staticField,Numerical->True]/(10^6);
	UP=MatrixExp[-I*tp*(H+2*\[Pi]*\[CapitalOmega]*{{0,1,0},{1,0,1},{0,1,0}}\[CircleTimes]IdentityMatrix[6]/Sqrt[2])];
	dt=Round[T*D/ndata]*(10^6/D);
	dU=MatrixExp[-I*dt*H];
	M=UP\[ConjugateTranspose].({{0,0,0},{0,1,0},{0,0,0}}\[CircleTimes]IdentityMatrix[6]).UP;
	\[Rho]={{0,0,0},{0,1,0},{0,0,0}}\[CircleTimes]IdentityMatrix[6]/6;
	\[Rho]=UP.\[Rho].UP\[ConjugateTranspose];
	\[Rho]=\[Rho]*(IdentityMatrix[3]\[CircleTimes]ConstantArray[1,{6,6}]);
	data=Table[\[Rho]=dU.\[Rho].dU\[ConjugateTranspose];{n*dt,Re@Tr[\[Rho].M]},{n,ndata}];
	model[t_]=offset+amp*Cos[2\[Pi] \[Omega] t];
	soln=FindFit[data,model[t],{{offset,.75},{amp,0.5},{\[Omega],\[Omega]elGuess/10^6}},{t}];
	{offset,amp,\[Omega]}/.soln
]


EnhancedLarmourSimulation[staticField_,A_,T_,\[CapitalOmega]_,doplot_:True,order_:-1]:=
	Module[
		{Hsim,stepsize,D,P0,rabidata,imperfectdata,perfectdata,pmod,psol,imod,isol,t\[Pi],UP,\[Rho]0,M0,fit,\[Omega]el,fig,U},
		(* Set the synthesizer to be at the leftmost of the four peaks *)
		D=10^6*NVCenterFrequency[Value@Cartesian[staticField],{0,0,0},1]-(10^6*NVZeemanSplitting[Value@Cartesian[staticField],{0,0,0},1]+Norm[A[[3]]])/2;
		stepsize=(10^6/D)/30;
		(* Define the Hamiltonian *)
		U=Module[{\[Theta],\[Phi]},
			\[Theta]=ArcTan[A[[3,1]],A[[3,2]]];
			\[Phi]=ArcTan[A[[3,3]],Norm[A[[3,{1,2}]]]];
			U=Subscript[\[DoubleStruckOne], 9]\[CircleTimes](MatrixExp[-I (\[Theta]/2) Y].MatrixExp[-I (\[Phi]/2) Z]);
		];
		If[order>=0,
			Hsim=U.Simplify[NVAverageHamiltonian[order,D,Carbon[A],Felton09Nucleus["14N"],StaticField->staticField,Numerical->True]/(10^6)].U\[ConjugateTranspose];,
			Hsim[t_]=Simplify[NVAverageHamiltonian[-1,D,Carbon[A],Felton09Nucleus["14N"],StaticField->staticField,Numerical->True][t/10^6]/(10^6)];
		];
		P0=S0\[CircleTimes]Subscript[\[DoubleStruckOne], 6];
		(* Run a rabi experiment to see where the CNOT is *)
		PrintTemporary["Determining CNOT time..."];
		rabidata=EvalPulse[Hsim,{{{.5/\[CapitalOmega],\[CapitalOmega]}},2\[Pi]{Sx\[CircleTimes]Subscript[\[DoubleStruckOne], 6]}},InitialState->P0/6,Observables->{P0},PollingInterval->30*stepsize,StepSize->stepsize];
		t\[Pi]=First@SelMin[First@N@Observables[rabidata,TimeVector->True],{2}];
		(* The imperfect CNOT *)
		PrintTemporary["Computing nonideal CNOT..."];
		UP=Last@Unitaries@EvalPulse[Hsim,{{{t\[Pi],\[CapitalOmega]}},2\[Pi]{Sx\[CircleTimes]Subscript[\[DoubleStruckOne], 6]}},StepSize->stepsize];
		\[Rho]0=UP.P0.UP\[ConjugateTranspose]/6;
		\[Rho]0=\[Rho]0*(IdentityMatrix[3]\[CircleTimes]ConstantArray[1,{6,6}]);
		M0=UP\[ConjugateTranspose].P0.UP;
		(* Gather the imperfect data *)
		PrintTemporary["Performing experiment with nonideal CNOT..."];
		imperfectdata=Module[{\[Rho]=\[Rho]0,U=Last@Unitaries@EvalPulse[Hsim,10^6/D,StepSize->stepsize],nsteps,dstep},
			nsteps=Round[10^-6*D*T];
			dstep=Round[nsteps/200];
			{{Observables,{Table[\[Rho]=U.\[Rho].U\[ConjugateTranspose];Re[Tr[\[Rho].M0]],{n,0,nsteps}][[1;;-1;;dstep]]}\[Transpose]},{TimeVector,Table[(10^6/D)*n,{n,0,nsteps}][[1;;-1;;dstep]]}}
		];
		{imod,isol}=NVCalibrationTools`Private`FitCosine[imperfectdata];
		(* The perfect CNOT *)
		UP=MatrixExp[-I (\[Pi]/2){{0,0,0},{0,0,1},{0,1,0}}\[CircleTimes]Subscript[\[DoubleStruckOne], 3]\[CircleTimes]{{1,0},{0,0}}];
		\[Rho]0=UP.P0.UP\[ConjugateTranspose]/6;
		\[Rho]0=\[Rho]0*(IdentityMatrix[3]\[CircleTimes]ConstantArray[1,{6,6}]);
		M0=UP\[ConjugateTranspose].P0.UP;
		(* Gather the perfect data *)
		PrintTemporary["Performing experiment with ideal CNOT..."];
		perfectdata=Module[{\[Rho]=\[Rho]0,U=Last@Unitaries@EvalPulse[Hsim,10^6/D,StepSize->stepsize],nsteps,dstep},
			nsteps=Round[10^-6*D*T];
			dstep=Round[nsteps/1000];
			{{Observables,{Table[\[Rho]=U.\[Rho].U\[ConjugateTranspose];Re[Tr[\[Rho].M0]],{n,0,nsteps}][[1;;-1;;dstep]]}\[Transpose]},{TimeVector,Table[(10^6/D)*n,{n,0,nsteps}][[1;;-1;;dstep]]}}
		];
		{pmod,psol}=NVCalibrationTools`Private`FitCosine[perfectdata];
		(* Use the formula to figure it out *)
		\[Omega]el=EnhancedLarmourFormula[Sequence@@(Value@Spherical[staticField]),A[[3,1]]/10^6,A[[1,1]]/10^6,A[[2,2]]/10^6,A[[3,3]]/10^6]/10^3;
		(* Plot or not *)
		fig=If[doplot,
			Column[{
				GraphicsRow[{
					ListPlot[Observables[rabidata,TimeVector->True],Joined->True,PlotLabel->"Rabi experiment -- chose time "<>ToString[1000*t\[Pi]]<>"ns"],
					Show[{
						ListPlot[Observables[imperfectdata,TimeVector->True],Joined->True,PlotRange->{0,1},PlotLabel->"Actual CNOT. Fit frequency: "<>ToString[isol[[3,2]]]<>"MHz\n(A,B)=("<>ToString[isol[[1,2]]]<>","<>ToString[isol[[2,2]]]<>")"],
						Plot[imod[t]/.isol,{t,0,Max@TimeVector@imperfectdata}, PlotRange->{0,1},PlotStyle->{Thick,Red}]
					}],
					Show[{
						ListPlot[Observables[perfectdata,TimeVector->True],Joined->True,PlotRange->{0,1},PlotLabel->"Ideal CNOT. Fit frequency: "<>ToString[psol[[3,2]]]<>"MHz\n(A,B)=("<>ToString[psol[[1,2]]]<>","<>ToString[psol[[2,2]]]<>")"],
						Plot[pmod[t]/.psol,{t,0,Max@TimeVector@perfectdata}, PlotRange->{0,1},PlotStyle->{Thick,Red}]
					}]
				},ImageSize->1000],
				TableForm[{Style[#,Bold]&/@{"Non-ideal Simulation","Ideal Simulation","\!\(\*SubscriptBox[\(\[Omega]\), \(L\)]\) formula"},{isol[[3,2]]*10^3,psol[[3,2]]*10^3,\[Omega]el}}]
			}]
		];
		{isol,psol,fig}
	]


Options[ELVHyperfineReconstruction]={"weights"->Automatic,Method->"NelderMead","\[Alpha]Constraints"->{-5,0},"A1Constraints"->{12,18},"A2Constraints"->{10,15}};
ELVHyperfineReconstruction[enhancedLarmourFreqs_,staticFields_,Asplit_,A\[Phi]_,OptionsPattern[]]:=
	Module[
		{\[Alpha],A1,A2,A3,cost,result,values,rotStaticFields,\[Alpha]c,A1c,A2c,A3c,B,\[Theta],\[Phi],data,weights},

		weights=OptionValue["weights"];
		If[weights===Automatic,weights=ConstantArray[1/Length[enhancedLarmourFreqs],Length[enhancedLarmourFreqs]]];

		rotStaticFields=Value[Spherical[#]]&/@staticFields;
		rotStaticFields[[All,2]]=rotStaticFields[[All,2]]-A\[Phi];

		(* We can measure the hyperfine splitting, which gives A3^2+a^2 *)
		A3=Sign[Asplit]*Sqrt[Asplit^2-\[Alpha]^2];

		(* Compiling the cost function gives a ten-fold speedup *)
		cost=Module[{expr,\[Alpha]\[Alpha],AA1,AA2,AA3},
			expr=Total[weights*(EnhancedLarmourFormula[Sequence@@#,\[Alpha]\[Alpha],AA1,AA2,A3/.\[Alpha]->\[Alpha]\[Alpha]]&/@rotStaticFields - enhancedLarmourFreqs)^2]/10^6;
			Compile[{\[Alpha],A1,A2},Evaluate[expr/.{\[Alpha]\[Alpha]->\[Alpha],AA1->A1,AA2->A2}]]
		];

		\[Alpha]c=OptionValue["\[Alpha]Constraints"];
		A1c=OptionValue["A1Constraints"];
		A2c=OptionValue["A2Constraints"];

		{result,values}=NMinimize[
			{
				Hold[cost[\[Alpha],A1,A2]],
				\[Alpha]c[[1]]<=\[Alpha]<=\[Alpha]c[[2]],
				A1c[[1]]<=A1<=A1c[[2]],
				A2c[[1]]<=A2<=A2c[[2]]
			},
			{\[Alpha],A1,A2},
			Method->OptionValue[Method],
			MaxIterations->3000
		];
		{Sqrt@result, {{A1,0,\[Alpha]},{0,A2,0},{\[Alpha],0,A3}}/.values}
	]


Options[EnhancedLarmourHyperfineReconstruction]={"weights"->Automatic,Method->"NelderMead","\[Alpha]Constraints"->{-5,0},"A1Constraints"->{12,18},"A2Constraints"->{10,15}};
EnhancedLarmourHyperfineReconstruction[enhancedLarmourFreqs_,staticFields_,Asplit_,A\[Phi]_,OptionsPattern[]]:=
	Module[
		{\[Alpha],A1,A2,A3,cost,result,values,rotStaticFields,\[Alpha]c,A1c,A2c,A3c,B,\[Theta],\[Phi],data,weights},

		weights=OptionValue["weights"];
		If[weights===Automatic,weights=ConstantArray[1/Length[enhancedLarmourFreqs],Length[enhancedLarmourFreqs]]];

		rotStaticFields=Value[Spherical[#]]&/@staticFields;
		rotStaticFields[[All,2]]=rotStaticFields[[All,2]]-A\[Phi];

		(* We can measure the hyperfine splitting, which gives A3^2+a^2 *)
		A3=Sign[Asplit]*Sqrt[Asplit^2-\[Alpha]^2];

		(* Compiling the cost function gives a ten-fold speedup *)
		cost=Module[{expr,\[Alpha]\[Alpha],AA1,AA2,AA3},
			expr=Total[weights*(EnhancedLarmourFormula[Sequence@@#,\[Alpha]\[Alpha],AA1,AA2,A3/.\[Alpha]->\[Alpha]\[Alpha]]&/@rotStaticFields - enhancedLarmourFreqs)^2]/10^6;
			Compile[{\[Alpha],A1,A2},Evaluate[expr/.{\[Alpha]\[Alpha]->\[Alpha],AA1->A1,AA2->A2}]]
		];

		\[Alpha]c=OptionValue["\[Alpha]Constraints"];
		A1c=OptionValue["A1Constraints"];
		A2c=OptionValue["A2Constraints"];

		{result,values}=NMinimize[
			{
				Hold[cost[\[Alpha],A1,A2]],
				\[Alpha]c[[1]]<=\[Alpha]<=\[Alpha]c[[2]],
				A1c[[1]]<=A1<=A1c[[2]],
				A2c[[1]]<=A2<=A2c[[2]]
			},
			{\[Alpha],A1,A2},
			Method->OptionValue[Method],
			MaxIterations->3000
		];
		{Sqrt@result, {{A1,0,\[Alpha]},{0,A2,0},{\[Alpha],0,A3}}/.values}
	]


EnhancedLarmourHyperfineErrorBars[enhancedLarmourFreqs_,enhancedLarmourStds_,staticFields_,staticFieldStds_,Asplit_,A\[Phi]_,N_,opt:OptionsPattern[EnhancedLarmourHyperfineReconstruction]]:=
	Module[
		{data,pp,entry},

		DistributeDefinitions[staticFields,staticFieldStds,enhancedLarmourFreqs,enhancedLarmourStds,Asplit,A\[Phi],opt,EnhancedLarmourHyperfineReconstruction];
		DistributeDefinitions["NVHamiltonian`"];

		data=ParallelTable[
			Module[{vectorNoise,elfNoise,noisyStaticFields,noisyELFs},
				vectorNoise[v_,\[Sigma]_]:=Vector[((RandomVariate[NormalDistribution[#,\[Sigma]]])&/@Value[Cartesian[v]]),Cartesian];
				elfNoise[\[Omega]_,\[Sigma]_]:=RandomVariate[NormalDistribution[\[Omega],\[Sigma]]];
				noisyStaticFields=Inner[vectorNoise,staticFields,staticFieldStds,List];
				noisyELFs=Inner[elfNoise,enhancedLarmourFreqs,enhancedLarmourStds,List];
				EnhancedLarmourHyperfineReconstruction[noisyELFs,noisyStaticFields,Asplit,A\[Phi],opt]
			],
			{n,N}
		];

		Print["Min Cost Value:  "<>ToString[Min[data[[All,1]]]]<>"kHz"<>"\n"<>
			"Max Cost Value:  "<>ToString[Max[data[[All,1]]]]<>"kHz"<>"\n"<>
			"Mean Cost Value: "<>ToString[Mean[data[[All,1]]]]<>"kHz"];

		pp[mean_,std_]:=ToString[mean]<>"\[PlusMinus]"<>ToString[std];
		entry[i_,j_]:=pp[Mean[data[[All,2,i,j]]],StandardDeviation[data[[All,2,i,j]]]];

		Print[MatrixForm[Array[entry,{3,3}]]]
		
	]


End[];


(* ::Section::Closed:: *)
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
