-----------------------------------------------------------
Installation
-----------------------------------------------------------

To install this library, just open up `Install.m` in 
Mathematica and goto

    Evaluation > Evaluate Notebook

Explanation of installation for those interested:
This installer simply moves the contents of the 
`src` folder into the search path of Mathematica,
which in this case, is the Applications subfolder of 

    $UserBaseDirectory

On Mac and Linux, simlinks are used for the convenience 
of developing code in `src`.
On Windows, we do an actual copy.

-----------------------------------------------------------
Usage
-----------------------------------------------------------

First of all, in every document that you want to 
use NVSim, the command
	Needs["NVSim`"]
needs to be included. (NVSim automatically imports 
quantum-utils if you were wondering.)
Until a better help system is included, Example.nb is 
the place to start.

A start has been made at a better documentation system. 
This documentation is installed when you run the 
Install.nb notebook. You can search for NVSim in the
documentation to get to it, or if you type in 
NVHamiltonian anywhere in any note book and press F1 
while the cursor is in the word, the documentation 
should popup.

This documentation is probably the best reference for 
the NVHamiltonian parameters.

-----------------------------------------------------------
Folder layout
-----------------------------------------------------------

The NVSim packages are kept in the **src** folder.

A couple of old files that some people may still want 
pieces of are in the **old** folder. This folder will 
be removed eventually.

Install.nb and README.txt will be the only files in the 
base directory eventually, all documentation will be
put into new folders.

The user-end documentation files are in the **nvdoc** 
folder.

The pre-deployed documentation source files are in the
**nvdocsrc** folder. I use Wolfram Workbench to deal 
with the documentation system.
