What is phaseField_PRISMS?

  phaseField_PRISMS is a C++ code implementing various phase field and
  mechanics models and its applications to some physical problems
  (spinodal decomposition, precipitate growth, coupled
  transport+mechanics, etc)

  It is a Finite Element Method (FEM) code and is build on top of the
  deal.II open source finite element library (http://www.dealii.org)
  
  This code is developed by the PRedictive Integrated Structural
  Materials Science (PRISMS) Center at University of Michigan which is
  supported by the U.S. Department of Energy, Office of Basic Energy
  Sciences, Division of Materials Sciences and Engineering under Award
  #DE-SC0008637

Installation and Use:

  1) Configure, compile and install the deal.II library. No additional
     libraries are required besides the installation of deal.ii with
     the following configuration flags:
     -DDEAL_II_WITH_MPI=ON  
     -DDEAL_II_WITH_P4EST=ON 
     -DDEAL_II_WITH_LAPACK=ON

     Download: http://www.dealii.org/download.html
     Detailed instructions: http://www.dealii.org/8.1.0/readme.html
      
  2) Fork the repo https://github.com/prisms-center/phaseField on
  GitHub (How to fork: https://help.github.com/articles/fork-a-repo/)
  and clone the phaseField repository using your GitHub username:

  $ git clone git@github.com:username/phaseField.git 
  (OR)
  $ git clone https://github.com/username/phaseField.git  

  [Note: phaseField is currently a private repository on GitHub and
  hence you need to be authorised to access the repository. Contact
  the developers/Mailing-List for requesting access]
   
  3) Running phaseField applications, for example cahnHilliard 
  $ cd applications/cahnHilliard
  $ cmake CMakeLists.txt
  For debug mode [default mode]:
  $ make debug
  For optimized mode:
  $ make release 
  Execution (serial runs):
  $ make run
  Execution (parallel runs):
  $ mpiexec -np nprocs ./main
  [here nprocs denotes the number of processors]
  
  4) Updates: Since phaseField code is still under active development,
  regular code and documentation updates are pushed to the upstream
  repo (https://github.com/prisms-center/phaseField) and we strongly
  recommend users to synchronize their respective forks at regular
  intervals or when requested by the developers through the
  announcements on the mailing list. 
  (How to sync: https://help.github.com/articles/syncing-a-fork/)

  5) Visualization: Results are written out in standard vtk/vtu/vts
  formats which can be visualized with the following open source
  applications:
  1. VisIt (https://wci.llnl.gov/simulation/computer-codes/visit/downloads)
  2. Paraview (http://www.paraview.org/download/)

Getting started:

  The applications are located under /applications and each
  application has a detailed README with instructions and related
  documentation.

  Applications are intended to serve as (1) Demonstration of the
  capabilities of this library, (2) Provide a framework for
  further development of specialized/advanced applications by
  users. Currently the following applications are provided:

  cahnHilliard
  allenCahn
  mechanics	
  coupledCahnHilliardAllenCahn
  betaPrimePrecipitate


  Applications under development/testing are preceded by an
  underscore, such as:
  _fickianFlux
  _coupledFickianFluxMechanics
  _coupledCahnHilliardMechanics


License:

  GNU Lesser General Public License (LGPL). Please see the file
  LICENSE for details.

Mailing List:
  
  PRISMS.phaseField@umich.edu

Further information, questions, issues and bugs:

  Contact the phaseField mailing list at PRISMS.phaseField@umich.edu  



