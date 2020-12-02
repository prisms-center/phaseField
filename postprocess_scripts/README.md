

## Python postprocessing scripts to work in conjunction with the VisIt CLI



### Installation of the visit CLI (MacOS)

1) Once VisIt is installed ([see installation instructions](https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/gui_manual/Intro/Installing_VisIt.html)) you need to add the path to the VisIt installation directory to you path, e.g.,

\$ export PATH=/Applications/VisIt.app/Contents/Resources/bin:$PATH

You can test if the cli works by typing (in a terminal window)

\$ visit -cli

This should open the cli environment and a VisIt window. 

2) [Sourcing scripts](https://www.visitusers.org/index.php?title=Using_CLI#Starting_the_CLI)

To run in batch mode, type:

\$ visit -cli -s scriptname.py

Or to run entirely in batch mode (with no window apperaring)

\$ visit -cli -nowin -s scriptname.py

Note that for full batch mode to work, the following lines need to be added to the script:

import sys
sys.exit()

#### phase_fraction.py

This script calculates the phase fraction (volume fraction for 2D, area fraction for 3D) based on the field variable "n". 

The phase fraction is calculated as the numerical integrals $\phi=\frac{1}{V}\int_V n\,d^3r $ for a 3D system and $\phi=\frac{1}{A}\int_A n\,d^2r $ for a 2D system, where V (A) is the total volume (area) of the system.

The results for all time states are outputted in the file **phi_vs_t.txt**.

For each row, the first column corresponds to the frame number, the second to the time and the third one to the phase fraction.

#### domain_count.py

This script counts the number of separate regions based on the values of field variable "n". Each interconnected region for which n<-0.5 is counted as a separate domain. The results for all time states are outputted in the file **nodoms_vs_t.txt**.

For each row, the first column corresponds to the frame number, the second to the time and the third one to the number of domains.

#### interface_area.py

This script calculates the total area (or length, in a 2D system) of the interface between two phases defined by the field variable "n".   

The numerical integral $A_{int}=2\int_{\delta\Omega} n\,d^2r $ is used to calculate the *area* of a 2D interface in a 3D system, where  $\delta\Omega$ is the contour surface at $n=0.5$.

The numerical integral  $L_{int}=2\int_{\delta l} n\,dr $ is used to the *length* of the 1D interface in a 2D system, where  $\delta l$  is the contour line at $n=0.5$.

The results for all time states are outputted in the file **iarea_vs_t.txt**. For each row, the first column corresponds to the frame number, the second to the time and the third one to the total interface area (or length).