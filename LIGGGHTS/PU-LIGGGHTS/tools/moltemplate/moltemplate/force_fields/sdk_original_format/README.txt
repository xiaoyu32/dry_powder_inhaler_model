The files in this directory are used to create the "sdk.lt" file
(containing SDK force field parameters for moltemplate).
These .PRM files are distributed with "EMC" written by Pieter J. in 't Veld.

Conversion from EMC (.PRM) format to moltemplate (.LT) format was 
done using the "emcprm2lt.py" script written by David Stelter.
Here is an example how to use the emcprm2lt.py script:

emcprm2lt.py --pair-style=lj/sdk/coul/long --bond-style=harmonic --angle-style=sdk sdk_lipids.prm sdk_cholesterol.prm --name=SDK_lipid+chol --units 

This will generate a file named "sdk.lt" which (in this example)
only includes the force field parameters for lipids and cholestrol.
Later you can define new molecules in moltemplate using:

import "sdk.lt"
NewMolecule inherits SDK {
  write("Data Atoms") {...atom coordinates and types go here...}
  write("Data Bond List") {...list of bonds goes here...}
}

This is only part of the SDK force field and is to be used for lipids
and cholesterol only.  Lipid parameters were taken from:

Shinoda et al. J. Phys. Chem. B, Vol. 114, No. 20, 2010
http://dx.doi.org/10.1021/jp9107206

Cholesterol parameters were taken from:
MacDermaid et al. J. Chem. Phys, 143(24), 243144, 2015
http://dx.doi.org/10.1063/1.4937153

You can define lipids with any topology built using the following types:

Name    Structure         Charge
NC      -CH2CH2-N-(CH3)3  +1
NH      -CH2CH2-NH3       +1
PH      -PO4-             -1
PHE     -PO4- (PE lipid)  -1
GL      -CH2CH-CH2-
EST1    -CH2CO2-
EST2    -H2CO2-
CMD2    -HC=CH- (cis)
CM      -CH2CH2CH2-
CT      CH3CH2CH2-
CT2     CH3CH2-
W       (H2O)3

This coarse-grainng allows for design of a wide variety of lipids. 

The following types are defined specifically for cholesterol:

Name    Structure       Location
C2T     -CH-(CH3)2      Tail
CM2     -CH2-CH2-       Tail
CM2R    -CH2-CH2-       Ring A
CMDB    -CH2-C=CH-      Ring A/B
CMB     -CH2-CH-CH-     Ring B/C
CMR     -CH-CH2-CH2-    Ring B/C
CMR5    -CH2-CH2-CH-    Ring D
CTB     -CH2-CH-CH3-    Tail
CTBA    -C-CH3          Ring A/B
CTBB    -C-CH3          Ring C/D
OAB     -CH-OH          Ring A

See the provided reference for details on the CG cholesterol topology. 
A 5.0-10.0 timestep is used when using cholesterol.

Several limiations, due to missing parameters:
-use of cholesterol with type "NH" is not possible.
-use of cholesterol with type "PHE" is not possible.

---- Credits: ----

emcprm2lt.py was written by David Stelter
EMC was written by Pieter J. in 't Veld
SDK was created by Shinoda, DeVane, Klein, J.Phys.Chem.B, Vol. 114, No. 20, 2010

---- additional citation request ----

Since we borrowed force field parameters from files distributed with EMC,
if you use files generated by "emcprm2lt.py", please also cite the EMC paper:
P. J. in ???t Veld and G. C. Rutledge, Macromolecules 2003, 36, 7358



