## Python postprocessing scripts to work in conjunction with the VisIt CLI

### Installation of the visit CLI (MacOS)

1) Once VisIt is installed ([see installation instructions](https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/gui_manual/Intro/Installing_VisIt.html)) you need to add the path to the VisIt installation directory to you path, e.g.,

\$ export PATH=/Applications/VisIt.app/Contents/Resources/bin:$PATH

You can test if the cli works by typing (in a terminal window)

\$ visit -cli

This should open the cli environment and a VisIt window. 

2) [Sourcing scripts](https://www.visitusers.org/index.php?title=Using_CLI#Starting_the_CLI)

To run in batch mode, type:

\$ visit -cli -s <scriptname.py>

Or to run entirely in batch mode (with no window apperaring)

\$ visit -cli -nowin -s <scriptname.py>

Note that for full batch mode to work, the following lines need to be added to the script (the first line at the beginning and the second at the end):

import sys
sys.exit()

### Postprocessing scripts

#### plot_and_save.py

This script creates a pseudocolor (in 2D) or contour (in 3D) plot for each time states and saves the serieas of plots as png images. The default plotting variable is "n".

#### phase_fraction.py

This script calculates the phase fraction (volume fraction for 2D, area fraction for 3D) based on the field variable "n". 

The phase fraction is calculated as the numerical integrals ![formula](https://render.githubusercontent.com/render/math?math=\phi=\frac{1}{V}\int_Vnd^3r) for 
a 3D system and ![formula](https://render.githubusercontent.com/render/math?math=\phi=\frac{1}{A}\int_An\,d^2r) for a 2D system, where *V* (*A*) is the total volume (area) of the system.

The results for all time states are outputted in the file **phi_vs_t.txt**.

For each row, the first column corresponds to the frame number, the second to the time and the third one to the phase fraction.

#### interface_area.py

This script calculates the total area (or length, in a 2D system) of the interface between two phases defined by the field variable "n".   

The numerical integral $A_{int}=2\int_{\delta\Omega}nd^2r$ ![formula](https://render.githubusercontent.com/render/math?math=A_{int}=2\int_{\delta\Omega}nd^2r) is used to calculate the *area* of a 2D interface 
in a 3D system, where ![formula](https://render.githubusercontent.com/render/math?math=\delta\Omega) is the contour surface at *n=0.5*.

The numerical integral ![formula](https://render.githubusercontent.com/render/math?math=L_{int}=2\int_{\delta%20l}ndr) is used to the *length* of the 1D 
interface in a 2D system, where ![formula](https://render.githubusercontent.com/render/math?math=\delta%20l) is the contour line at *n=0.5*.

The results for all time states are outputted in the file **iarea_vs_t.txt**. 
For each row, the first column corresponds to the frame number, the second to the time and the third one to the total interface area (or length).

#### domain_stats.py

This script calculates the following quantities for each of the time states: 
1) Number of domains, defined as the number or separate regions based on the values of field variable, "n". Each interconnected region for which n>0.5 is counted as a separate domain. 
2) Average domain size (area in 2D and volume in 3D)
3)  Standard deviation of the domain sizes. 

Results for all time states are outputted in the file **domain_stats_vs_t.txt**.

For each row, the first column corresponds to the frame number, the second to the time,  the third column to the number of domains, the fourth column to the average size, and the fifth column to the standard deviation of the size.

#### splitvtufiles.py

This script splits the output .vtu files into several files, allowing for parallel visualization in VisIt or ParaView. To run it simply 
1. Open a terminal

2. Go to the directory where the series of output files from a single simulation are located (these files must be generated using the single file per process option, see the **Output** section from [this page](https://prisms-center.github.io/phaseField/doxygen_files/input_file.html))

3. Type:

  \$ python splitvtufiles.py <no. files>

  where <no. files> is the number of files into which you want to split each solution\*.vtu file.
