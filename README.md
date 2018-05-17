-----------------------------------------------------------
INSTALLATION
-----------------------------------------------------------

If you are looking at this file, you have probably 
figured out GIT. You will also need to grab the 
latest quantum-utils version off SVN and install it.

To install NVSim, just open up Install.m in Mathematica 
and goto
	Evaluation > Evaluate Notebook

Explanation of Installation for those interested:
All "installation" does is put the packages where
Mathematica expects to find them. When the Needs[] 
command is called, Mathematica looks in two places:
	$UserBaseDirectory
	$BaseDirectory
Since $BaseDirectory is often a protected folder needing 
root permissions, the installation puts NVSim into
$UserBaseDirectory. On Mac and Linux, simlinks are used.
On Windows, we need to literally copy the file there
(this is why Install.nb needs to be run everytime a
Windows user grabs the newest version off SVN).

-----------------------------------------------------------
WHERE TO BEGIN
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
FOLDER LAYOUT
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
