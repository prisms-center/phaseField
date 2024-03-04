---
layout: default
---
[![PRISMS-PF Logo](assets/logo.png)](https://prisms-center.github.io/phaseField/)

***
## Overview
PRISMS-PF is a powerful, massively parallel finite element code for conducting phase field and other related simulations of microstructural evolution. The phase field method is commonly used for predicting the evolution if microstructures under a wide range of conditions and material systems. PRISMS-PF provides a simple interface for solving customizable systems of partial differential equations of the type commonly found in phase field models, and has 29 pre-built application modules, including for precipitate evolution, grain growth, dendritic solidification, corrosion and spinodal decomposition.

![PRISMS-PF Example Results](assets/example_bar.png)

***
## Announcements
- 7/23/2021: [Version 2.3 released.](https://github.com/prisms-center/phaseField/releases/tag/v2.3) This version contains moderate changes from version 2.2, including new applications, new postprocessing scripts, improvements in performance, bug fixes and comparibility with the latest version of deal.II. For a detailed discussion of the new features, bug fixes and other changes in this version, please consult the [version_changes.md](https://github.com/prisms-center/phaseField/blob/master/version_changes.md) document.
- 8/13/2023: We are pleased to announce that PRISMS-PF is now compatible with deal.II 9.4.2 or older (up to 9.2.0). We are in the process of ensuring compatibility with the latest deal.II version (9.5.1). Note that from now on the oldest deal.II version compatible with PRISMS-PF is 9.2.0. You can find all of the deal.II available versions for different platforms [here](https://www.dealii.org/downloads/). Please report any issues [here](https://github.com/prisms-center/phaseField/issues) or [contact us](mailto:prisms-pf@umich.edu) directly.
- 7/16/2023: We are pleased to announce the release of a new application, alloySolidification_uniform, to simulate solidification of a binary alloy at uniform temperature. In the model implemented [A. Karma, Phys. Rev. Lett. 87, 115701 (2001)] latent heat is assumed to diffuse much faster than impurities and, therefore, the temperature field is considered to be fixed by external conditions. In contrast to alloySolidification, this application considers solidification under uniform temperature and no solute diffusion in the solid. The implementation of the model in this application has been validated by comparing them the results to those reported in Figure 1 from the paper cited above.
- 7/16/2023: 2023 PRISMS Workshop and Training Sessions. Everyone interested in PRISMS-PF (or in materials modeling, in general) is encouraged to attend the 2023 PRISMS Center Workshop (August 21-25), returning this year as an in person event. The training sessions will be on 21-23 of August and the symposia will be on 24-25 of August. The training part will include two sessions dedicated to a hands-on tutorial for PRISMS-PF on Monday Aug. 21:\
10AM - 12PM (EDT)\
1PM - 3PM (EDT)\
For more information and to register please visit [this page](http://www.prisms-center.org/#/workshop) (all attendees must register). There will also be training sessions for other PRISMS frameworks (see registration page for details).
- 6/30/2022: New Microgalvanic Corrosion (corrosion_microgalvanic) application release! This application simulates the evolution of the metal-electrolyte interface during free immersion due to the microgalvanic coupling between the anodic and cathodic metals. A manuscript detailing the implementation of the model, the scaling performance of the code, and simulation examples in 2D and 3D has been submitted to MRS Communications (currently under revision). 
- 7/26/2021: 2021 PRISMS Virtual Workshop and Training Sessions. Everyone interested in PRISMS-PF (or in materials modeling, in general) is invited to attend the 2021 PRISMS Center Workshop (August 3-6) and training sessions (August 9-13). Both of these events will be held virtually this year. The training part will include two sessions dedicated to a hands-on tutorial for PRISMS-PF on the following dates:
Monday, Aug. 9, 10 AM -12 PM (EDT)
Wednesday, Aug. 11, 10 AM -12 PM (EDT). 
Registration for both the workshop and the training is free! Please register on [this page](https://docs.google.com/forms/d/e/1FAIpQLSfLY5TxU768H2XW0FJuXQtqarpoW0NWLJ9p0ucqAEgXQCaNQA/viewform) (all attendees must register). There will also be training sessions for other PRISMS frameworks (see registration page for details).
- 7/23/2021: PRISMS-PF Docker Image Update. We have updated the PRISMS-PF Docker image to use the current version of PRISMS-PF (2.2).  This new image solves the discrepancy between the previous input parameters filename extension (parameters.in) and the current filename for each application (parameters.prm).  The image is availabel in the PRISMS-PF [Docker repository](https://hub.docker.com/r/prismspf/prismspf). Please follow the updated [instructions for the Docker installation](https://prisms-center.github.io/phaseField/doxygen_files/install_prismspf.html) to download and launch containers using the new image.
- 7/23/2021: [Version 2.2 released.](https://github.com/prisms-center/phaseField/releases/tag/v2.2) This version contains moderate changes from version 2.1.2, the main one being the release of new applications. For a detailed discussion of the new features, bug fixes and other changes in this version, please consult the [version_changes.md](https://github.com/prisms-center/phaseField/blob/master/version_changes.md) document.
- 7/6/2021: New Spinodal Decomposition (spinodalDecomposition) application release! This application simulates domain coarsening in time following spinodal decomposition, which is goverened by Cahn-Hilliard dynamics. In contrast to the Cahn-Hilliard application, the initial condition for the concentration is set to uniformly distributed random values around a concentration within the spinodal region. 
- 3/15/2021: New Alloy Solidification (alloySolidification) application release! This application simulates the directional solidification of a binary alloy  in the dilute limit, taking into account solute segregation, which leads to dendritic growth. The model employed [Echebarria et al., Phys. Rev. E 70, 061604 (2004)] includes an anti-trapping solute current required to correct for spurious effects that arise from considering an interface thickness much larger than the physical solid-liquid interface.
- 3/8/2021: We have uploaded new video tutorial covering the installation of prerequisited for PRISMS-PF. This video can be watched on the new [PRISMS Center YouTube channel](https://www.youtube.com/channel/UCZXc3007JuBCGKDcneD_umA) and is also available in the [Tutorials](https://github.com/prisms-center/phaseField/blob/gh-pages/pages/tutorial.md) section.
- 12/1/2020: We have uploaded new directory to the PRISMS-PF GitHub repository which includes scripts to perform processing tasks on output data files generated by PRISMS-PF. These scripts are written in Python and work in combination with the [VisIt CLI](https://www.visitusers.org/index.php?title=Using_CLI#Starting_the_CLI). Currently, we have three scripts to calculate: (1) phase fraction, (2) number of domains, and (3) area of an interface but more will be added soon. A description of each of these scripts as well as installation instructions for the VisIt CLI can be found [here](https://github.com/prisms-center/phaseField/tree/master/postprocess_scripts). 
- 11/24/2020: We have uploaded new video tutorial covering a simple example on how to set up and run a simulation for explicit nucleation and growth in PRISMS-PF. This video can be watched on the new [PRISMS Center YouTube channel](https://www.youtube.com/channel/UCZXc3007JuBCGKDcneD_umA) and is also available in the [Tutorials](https://github.com/prisms-center/phaseField/blob/gh-pages/pages/tutorial.md) section.
- 10/13/2020: We have updated PRISMS-PF to make it compatible with the latest [deal.II](https://www.dealii.org/) version (9.2.0). The main change from the previous deal.II version (9.1.2) is that now only input files with the ".prm" extension are accepted. In order to comply with that, the default name for the input file in all PRISMS-PF applications has been changed to "parameters.prm" ("parameters.in" will no longer work). After cloning the latest PRISMS-PF repository, make sure to use the correct extension for input files when creating new applications. Please report any issues [here](https://github.com/prisms-center/phaseField/issues).
- 9/25/2020: The PRISMS-PF [GUI-based nanoHUB tool](https://nanohub.org/tools/prismspfmisfit), which calculates the equilibrium shape of a precipitate particle due to anisotropic interfacial energy and misfit strain, has been updated for improved performance.
- 8/16/2020: A new video tutorial covering PRISMS-PF installation, running of a sample application and results visualization has been uploaded to the new [PRISMS Center YouTube channel](https://www.youtube.com/channel/UCZXc3007JuBCGKDcneD_umA). It is also visible in the [Tutorials](https://github.com/prisms-center/phaseField/blob/gh-pages/pages/tutorial.md) section. More PRISMS-PF tutorial videos will be available soon. 
- 8/1/2020: We are working on resolving compatibility issues between PRISMS-PF and the latest version of deal.ii (9.2.0). Until we have a working version, please use [deal.ii version 9.1.1](https://www.dealii.org/9.1.1/index.html). We will send an update when these issues are resolved.
- 6/23/2020: PRISMS-PF Docker Image Update. We have updated the PRISMS-PF Docker image to use the current version of PRISMS-PF (2.1.2). This update should solve the compatibility issues between the core library and the applications from the previous image. We have also created a [new repository](https://hub.docker.com/r/prismspf/prismspf) on Docker Hub where the new image, as well as all upcoming images, will be stored. Please follow the updated [instructions for the Docker installation](https://prisms-center.github.io/phaseField/doxygen_files/install_prismspf.html) to download and launch containers using the new image.
- 6/13/2020: New corrosion application release! This application simulates the evolution of the metal-electrolyte interface during the anodic corrosion reaction. The model employed [Chadwick et al., J. Electrochem. Soc.,10, C633-C646 (2018)] uses the phase-field and smoothed-boundary methods to track the moving metal/electrolyte interface and to couple it to mass transport (diffusion and migration) within the electrolyte and Butler-Volmer electrochemical kinetics. 
- 3/26/2020: A new open-access article titled [PRISMS-PF: A general framework for phase-field modeling with a matrix-free finite element method](https://www.nature.com/articles/s41524-020-0298-5) has been published on *npj Computational Materials*. In this article we introduce PRISMS-PF as a new open-source framework for phase-field modeling, emphasizing ease of use, flexibility and adaptability to a wide range of applications, and computational performance. We recommend this article to everyone interested in an overview of the PRISMS-PF framework and its applications. Please use this publication as the standard reference when citing PRISMS-PF. We'll keep an eye on references to that paper, but if you use PRISMS-PF in your work, please send us the citation so that we can advertise your work on the PRISMS-PF website and can pass usage information to our funding agency.

([older announcements](pages/announcements.md))

***

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

### Citing PRISMS-PF

Please cite [the following reference](https://www.nature.com/articles/s41524-020-0298-5) when discussing PRISMS-PF in a publication:

S. DeWitt, S. Rudraraju, D. Montiel, W.B. Andrews, and K. Thornton. PRISMS-PF: A general framework for phase-field modeling with a matrix-free finite element method. _npj Computational Materials_ __6__, 29 (2020).

If additionally you would like to cite a specific release of PRISMS-PF, please use the following format:

PRISMS-PF, v2.1.2 (2019). Available from https://github.com/prisms-center/phaseField. DOI: 10.5281/zenodo.3357005.

For DOI information for other releases, please refer to [this site](https://zenodo.org/record/3357005).

***

## General Links
- [PRISMS Center homepage](http://www.prisms-center.org/#/home) <br>
- [PRISMS-PF code repository](https://github.com/prisms-center/phaseField) <br>
- [User registration link](http://goo.gl/forms/GXo7Im8p2Y)
- [YouTube channel](https://www.youtube.com/channel/UCZXc3007JuBCGKDcneD_umA/playlists)
- [Virtual Machine](http://www.prisms-center.org/#/ctools/software) (contains most of the PRISMS Center tools, including PRISMS-PF)

***
## Getting Started
- [Quick start guide](https://github.com/prisms-center/phaseField#quick-start-guide) <br>
- [Installation instructions](doxygen_files/install.html) <br>
- [User manual](doxygen_files/manual.html) <br>
- [Repository of training slides and exercises](https://goo.gl/BBTkJ8)
- [Tutorials](https://github.com/prisms-center/phaseField/blob/gh-pages/pages/tutorial.md)

***
## Getting Help
- [PRISMS-PF forum](https://groups.google.com/forum/#!forum/prisms-pf-users) <br>
- [Code documentation](doxygen_files/index.html) <br>
- [Email the developers](mailto:prismsphasefield.dev@umich.edu)

***
## Publications That Use PRISMS-PF

[Wang and Liu, Multiscale thermo-kinetic characterization for β′ and β1 precipitation in Mg-Sm alloys, Acta Mater. 254, 119011 (2023)](https://doi.org/10.1016/j.actamat.2023.119011)

[Pendl and Hochrainer, Coupling stress fields and vacancy diffusion in phase-field models of voids as pure vacancy phase, Comput. Mater. Sci. 224, 112157 (2023)](https://doi.org/10.1016/j.commatsci.2023.112157)

[Goel, Lyu, DeWitt, Montiel, and Thornton, Simulating microgalvanic corrosion in alloys using the PRISMS phase-field framework, MRS Communications 12, 1050–1059 (2022)](https://doi.org/10.1557/s43579-022-00266-6)

[Bhagat and Rudraraju, Modeling of dendritic solidification and numerical analysis of the
phase-field approach to model complex morphologies in alloys, arXiv preprint (2022)](https://doi.org/10.48550/arXiv.2210.14449)

[Kinzer and Chandran, A Phase-Field Study on the Effects of Nanoparticles on Solidification and Grain Growth, arXiv preprint (2022)](https://doi.org/10.48550/arXiv.2207.07153)

[Gao, Wang, Li, et al,, Cerium-alloyed ultra-high strength maraging steel with good ductility: Experiments, first-principles calculations and phase-field simulations, Materials Science and Engineering: A 846, 14330 (2022)](https://doi.org/10.1016/j.msea.2022.143306)

[Yao, Montiel, and Allison, Investigating the Effects of Dendrite Evolution on Microsegregation in Al–Cu Alloys by Coupling Experiments, Micro-modeling, and Phase-Field Simulations. Metall Mater Trans A 53, 3341–3356 (2022)](https://doi.org/10.1007/s11661-022-06748-5)

[Cao, Zhang, Meng, and Zhang, Analyzing effects of temperature gradient and scan rate on metal additive manufacturing microstructure by using phase field-finite element method, Modelling Simul. Mater. Sci. Eng. 30, 034003 (2022)](https://iopscience.iop.org/article/10.1088/1361-651X/ac4f3a)

[Brewick, Simulating Pitting Corrosion in AM 316L Microstructures Through Phase Field Methods and Computational Modeling, J. Electrochem. Soc. 169, 011503 (2022)](https://iopscience.iop.org/article/10.1149/1945-7111/ac4935)

\*[DeWitt, Rudraraju, Montiel, Andrews and Thornton, PRISMS-PF: A general framework for phase-field modeling with a matrix-free finite element method, npj Comput Mater 6, 29 (2020)](https://www.nature.com/articles/s41524-020-0298-5)

[Aagesen, Adams, Allison, et al., PRISMS: An Integrated, Open-Source Framework for Accelerating Predictive Structural Materials Science,  JOM 70, 2298–2314 (2018)](https://link.springer.com/article/10.1007%2Fs11837-018-3079-6)

[DeWitt, Solomon, Natarajan, Araullo-Peters, Rudraraju, Aagesen, Puchala, Marquis, Van der Ven, Thornton, and Allison, Misfit-driven β′′′ precipitate composition and morphology in Mg-Nd alloys, Acta Materialia, 136, 378-389 (2017)](https://www.sciencedirect.com/science/article/pii/S1359645417305281)

\* Please use this publication as the standard reference when citing PRISMS-PF.

Note: If you use PRISMS-PF in one of your publications, please send the publication information to [prismsphasefield.dev@umich.edu](mailto:prismsphasefield.dev@umich.edu) to help us demonstrate our impact to our funding agency.

***
## Publications That Cite PRISMS-PF

[Gao, Peng, Zhang, et al., Profound strengthening and toughening effect of reinforcement aspect ratio in composite with network architecture, Journal of Alloys and Compounds, 167444 (2022)](https://doi.org/10.1016/j.jallcom.2022.167444)

[Endo, Matsuda, Tanaka, et al., A phase-field model by an Ising machine and its application to the phase-separation structure of a diblock polymer, Sci Rep 12, 10794 (2022)](https://doi.org/10.1038/s41598-022-14735-4)

[Stewart, Recent progress on the mesoscale modeling of architected thin-films via phase-field formulations of physical vapor deposition, Computational Materials Science 211, 111503 (2022)](https://doi.org/10.1016/j.commatsci.2022.111503)

[Hong and Viswanathan, Open-Sourcing Phase-Field Simulations for Accelerating Energy Materials Design and Optimization, ACS Energy Letters 2020 5, 3254-3259 (10)](https://pubs.acs.org/doi/abs/10.1021/acsenergylett.0c01904)

[DeWitt and Thornton, Phase Field Modeling of Microstructural Evolution in Computational Materials System Design, Shin and Saal, Eds., Springer Nature, London (2018)](https://link.springer.com/chapter/10.1007/978-3-319-68280-8_4)

[Yaghoobi, Ganesan, Sundar, Lakshmanan, Rudraraju, Allison, and Sundararaghavan, PRISMS-Plasticity: An open-source crystal plasticity finite element software, Computational Materials Science 169, 109078 (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0927025619303696)

[Wheeler, Keller, DeWitt, Jokisaari, Schwen, Guyer, Aagesen, Heinonen, Tonks, Voorhees, and Warren, 2019. PFHub: The Phase-Field Community Hub. Journal of Open Research Software, 7(1), 29 (2019)](https://openresearchsoftware.metajnl.com/article/10.5334/jors.276/)

[Hötzer, Reiter, Hierl, Steinmetz, Selzer, and Nestler, The parallel multi-physics phase-field framework Pace3D, Journal of Computational Science 26, 1-12 (2018)](https://www.sciencedirect.com/science/article/abs/pii/S1877750317310116)

[Tonks and Aagesen, The Phase Field Method: Mesoscale Simulation Aiding Material DiscoveryAnnual Review of Materials Research 49, 79–102 (2019)](https://www.annualreviews.org/doi/abs/10.1146/annurev-matsci-070218-010151)
