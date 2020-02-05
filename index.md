---
layout: default
---
[![PRISMS-PF Logo](assets/logo.png)](https://prisms-center.github.io/phaseField/)

***
## Overview
PRISMS-PF is a powerful, massively parallel finite element code for conducting phase field and other related simulations of microstructural evolution. The phase field method is commonly used for predicting the evolution if microstructures under a wide range of conditions and material systems. PRISMS-PF provides a simple interface for solving customizable systems of partial differential equations of the type commonly found in phase field models, and has 24 pre-built application modules, including for precipitate evolution, grain growth, dendritic solidification, and spinodal decomposition.

![PRISMS-PF Example Results](assets/example_bar.png)

### Features and Capabilities
- Matrix-free finite element framework for improved performance over traditional finite element approaches
- Parallelization at the inter-node, intra-node, and intra-core levels (MPI, threads, vectorization), with near ideal scaling beyond 1,000 cores
- Adaptive meshing to greatly reduce problem sizes
- Support for high order elements, with up to 5th order spatial accuracy
- Support for explicit nucleus placement to enable simulations that include nucleating phases
- Grain-remapping algorithm to facilitate simulations of polycrystals with thousands of grains
- Simple interface to solve an arbitrary number of coupled PDEs
- Straightforward Docker-based installation

### Acknowledgements
This code is developed by the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan which is supported by the U.S. Department of Energy (DOE), Office of Basic Energy Sciences, Division of Materials Sciences and Engineering under Award #DE-SC0008637.

### License
PRISMS-PF is released under the GNU Lesser General Public License (LGPL).

***
## Announcements
- 1/14/2020: The PRISMS Center will hold a workshop at the TMS 2020 149th Annual Meeting and Exhibition, where participants will become familiar with PRISMS software tools, including PRISMS-PF, and will gain the knowledge about when to use these tools. For more information and registration [visit this site](https://www.tms.org/TMS2020/Professional_Development/Workshops_and_Courses/PRISMS__Training_Workshop/TMS2020/pd/PRISMS_Workshop.aspx?hkey=7238a67d-383c-4223-880d-8ea3d87f4729&_zs=DiBXC1&_zl=xeeV5).

- 12/2/2019: The PRISMS Center at the University of Michigan is looking to hire a new research faculty member to lead the development of PRISMS-PF and participate in phase field research related to the center's goals. If you are interested in applying, [please visit this page](https://careers.umich.edu/job_detail/180053/research_investigatorasst_research_scientist).

- 6/11/2019: The registration for the 2019 PRISMS Center Workshop is now open. The workshop is August 5-9 in Ann Arbor, MI. The first three days are dedicated to hands-on training with the PRISMS codes (including a day for PRISMS-PF). The final two days are a series of technical presentations by members of the PRISMS Center and distinguished outside guests. You can register [here](http://www.prisms-center.org/#/workshop).

- 6/10/2019: A [new GUI-based nanoHUB tool](https://nanohub.org/tools/prismspfmisfit) has been released that is powered by PRISMS-PF. It calculates the equilibrium shape of a precipitate due to anisotropic interfacial energy and misfit strain.

- 9/11/2018: Please take a moment to fill out the short [2018 PRISMS-PF User Survey](https://goo.gl/forms/rAp8cJAeBjqsi5ep1). The survey results will guide future PRISMS-PF development and help us justify continued funding.

- 8/21/2018: The tentative dates for the 2019 PRISMS Workshop are August 5-9 in Ann Arbor, MI. The first three days are dedicated to hands-on training with the PRISMS codes (including a day for PRISMS-PF). The final two days are a series of technical presentations by members of the PRISMS Center and distinguished outside guests. For more information, please send a message to [prismsphasefield.dev@umich.edu](mailto:prismsphasefield.dev@umich.edu).

- 8/21/2018: [Version 2.1 released.](https://github.com/prisms-center/phaseField/releases/tag/v2.1) This is a moderate-level update to v2.0.2. The structure of equations and ICs/BCs files has been updated for improved legibility and flexibility. A new grain-remapper has been introduced to handling polycrystalline simulations with many grains and a new hybrid Newton/Picard nonlinear solver has been added.

([older announcements](pages/announcements.html))

***
## General Links
- [PRISMS Center homepage](http://www.prisms-center.org/#/home) <br>
- [PRISMS-PF code repository](https://github.com/prisms-center/phaseField) <br>
- [User registration link](http://goo.gl/forms/GXo7Im8p2Y)

***
## Getting Started
- [Quick Start Guide](https://github.com/prisms-center/phaseField#quick-start-guide) <br>
- [Installation instructions](doxygen_files/install.html) <br>
- [User manual](doxygen_files/manual.html) <br>
- [Repository of training slides and exercises](https://goo.gl/BBTkJ8)

***
## Getting Help
- [PRISMS-PF forum](https://groups.google.com/forum/#!forum/prisms-pf-users) <br>
- [Code documentation](doxygen_files/index.html) <br>
- [Email the developers](mailto:prismsphasefield.dev@umich.edu)

***
## Publications Citing PRISMS-PF

[Aagesen, Adams, Allison, et al., PRISMS: An Integrated, Open-Source Framework for Accelerating Predictive Structural Materials Science, JOM, 71 (2018)](https://link.springer.com/article/10.1007%2Fs11837-018-3079-6)

[DeWitt and Thornton, Phase Field Modeling of Microstructural Evolution in Computational Materials System Design, Shin and Saal, Eds., Springer Nature, London (2018)](https://link.springer.com/chapter/10.1007/978-3-319-68280-8_4)

[DeWitt, Solomon, Natarajan, Araullo-Peters, Rudraraju, Aagesen, Puchala, Marquis, Van der Ven, Thornton, and Allison, Misfit-driven β′′′ precipitate composition and morphology in Mg-Nd alloys, Acta Materialia, 136 (2017)](https://www.sciencedirect.com/science/article/pii/S1359645417305281)


(Note: If you use PRISMS-PF in one of your publications, please send the publication information to [prismsphasefield.dev@umich.edu](mailto:prismsphasefield.dev@umich.edu) to help us demonstrate our impact to our funding agency.)
