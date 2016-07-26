//implementation of various models of anisotropy for
//st.venant-kirchoff material model of elasticity

#ifndef ANISOTROPY_MECHANICS_H
#define ANISOTROPY_MECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include <deal.II/base/table.h>
#include <deal.II/base/conditional_ostream.h>

//Each material model is characterized by the number of independent
//constants required to characterize its elasticity tensor in the Voigt
//notation

//3D models:
//ISOTROPIC - 2 constants [E, nu], where E-modulus and nu-poisson's ratio
//TRANSVERSE- 5 constants [C11 C33 C44 C12 C13]
//ORTHOTROPIC- 9 constants [C11 C22 C33 C44 C55 C66 C12 C13 C23]
//ANISOTROPIC- 21 constants [C11 C22 C33 C44 C55 C66 C12 C13 C14 C15
//C16 C23 C24 C25 C26 C34 C35 C36 C45 C46 C56]

//2D models:
//ISOTROPIC- (Plane Strain) 2 constants [E, nu]
//ANISOTROPIC- 6 constants [C11 C22 C33 C12 C13 C23]

//1D models:
//ISOTROPIC- 1 constant [E]

enum elasticityModel {ISOTROPIC, TRANSVERSE, ORTHOTROPIC, ANISOTROPIC, ANISOTROPIC2D};

template <int dim>
void getCIJMatrix(elasticityModel model, double constants[], dealii::Table<2, double>& CIJ, dealii::ConditionalOStream& pcout){
  CIJ.fill(0.0);
  pcout << "Reading material model:";
  switch (dim){
  case 1:{
    pcout << " 1D ";
    //1D models
    switch (model){
    case ISOTROPIC:{  
      pcout << " ISOTROPIC \n"; 
      CIJ[0][0]=constants[0];
      break;
    }
    default:{
      std::cout << "\nelasticityModels: Supported models in 1D - ISOTROPIC\n"; 
      std::cout << "See /src/elasticityModels.h\n";       
      exit(-1);
    }
    }
    break;
  }
  case 2:{
   pcout << " 2D ";
    //2D models
    switch (model){
    case ISOTROPIC:{
      pcout << " ISOTROPIC \n"; 
      double E=constants[0], nu=constants[1];
      double  mu=E/(2*(1+nu)), lambda= nu*E/((1+nu)*(1-2*nu)); 
      CIJ[0][0]=lambda+2*mu;
      CIJ[1][1]=lambda+2*mu;
      CIJ[2][2]=mu;
      CIJ[0][1]=CIJ[1][0]=lambda;
      break;
    }
    case ANISOTROPIC:{
      pcout << " ANISOTROPIC \n"; 
      CIJ[0][0]=constants[0]; //C11
      CIJ[1][1]=constants[1]; //C22
      CIJ[2][2]=constants[2]; //C33
      CIJ[0][1]=CIJ[1][0]=constants[3]; //C12
      CIJ[0][2]=CIJ[2][0]=constants[4]; //C13
      CIJ[1][2]=CIJ[2][1]=constants[5]; //C23
      break;
    }
    default:{
      std::cout << "\nelasticityModels: Supported models in 2D - ISOTROPIC/ANISOTROPIC\n"; 
      std::cout << "See /src/elasticityModels.h\n";       
      exit(-1);
    }
    }
    break;
  }
  case 3:{
    pcout << " 3D ";
    //3D models
    switch (model){
    case ISOTROPIC:{
      pcout << " ISOTROPIC \n"; 
      double E=constants[0], nu=constants[1];
      double  mu=E/(2*(1+nu)), lambda= nu*E/((1+nu)*(1-2*nu)); 
      CIJ[0][0]=lambda+2*mu;
      CIJ[1][1]=lambda+2*mu;
      CIJ[2][2]=lambda+2*mu;
      CIJ[3][3]=mu;
      CIJ[4][4]=mu;
      CIJ[5][5]=mu;
      CIJ[0][1]=CIJ[1][0]=lambda;
      CIJ[0][2]=CIJ[2][0]=lambda;
      CIJ[1][2]=CIJ[2][1]=lambda;
      break;
    }
    case TRANSVERSE:{
      pcout << " TRANSVERSE \n"; 
      CIJ[0][0]=constants[0]; //C11
      CIJ[1][1]=constants[0]; //C11
      CIJ[2][2]=constants[1]; //C33
      CIJ[3][3]=constants[2]; //C44
      CIJ[4][4]=constants[2]; //C44
      CIJ[5][5]=(constants[0]-constants[3])/2.0; //(C11-C12)/2
      CIJ[0][1]=CIJ[1][0]=constants[3]; //C12 
      CIJ[0][2]=CIJ[2][0]=constants[4]; //C13 
      CIJ[1][2]=CIJ[2][1]=constants[4]; //C13 
      break;
    }
    case ORTHOTROPIC:{
      pcout << " ORTHOTROPIC \n"; 
      CIJ[0][0]=constants[0]; //C11
      CIJ[1][1]=constants[1]; //C22
      CIJ[2][2]=constants[2]; //C33
      CIJ[3][3]=constants[3]; //C44
      CIJ[4][4]=constants[4]; //C55
      CIJ[5][5]=constants[5]; //C66
      CIJ[0][1]=CIJ[1][0]=constants[6]; //C12 
      CIJ[0][2]=CIJ[2][0]=constants[7]; //C13 
      CIJ[1][2]=CIJ[2][1]=constants[8]; //C23 
      break;
    }
    case ANISOTROPIC:{
      pcout << " ANISOTROPIC \n"; 
      CIJ[0][0]=constants[0]; //C11
      CIJ[1][1]=constants[1]; //C22
      CIJ[2][2]=constants[2]; //C33
      CIJ[3][3]=constants[3]; //C44
      CIJ[4][4]=constants[4]; //C55
      CIJ[5][5]=constants[5]; //C66
      CIJ[0][1]=CIJ[1][0]=constants[6]; //C12 
      CIJ[0][2]=CIJ[2][0]=constants[7]; //C13 
      CIJ[0][3]=CIJ[3][0]=constants[8]; //C14 
      CIJ[0][4]=CIJ[4][0]=constants[9]; //C15 
      CIJ[0][5]=CIJ[5][0]=constants[10]; //C16
      CIJ[1][2]=CIJ[2][1]=constants[11]; //C23 
      CIJ[1][3]=CIJ[3][1]=constants[12]; //C24 
      CIJ[1][4]=CIJ[4][1]=constants[13]; //C25 
      CIJ[1][5]=CIJ[5][1]=constants[14]; //C26
      CIJ[2][3]=CIJ[3][2]=constants[15]; //C34 
      CIJ[2][4]=CIJ[4][2]=constants[16]; //C35 
      CIJ[2][5]=CIJ[5][2]=constants[17]; //C36
      CIJ[3][4]=CIJ[4][3]=constants[18]; //C45 
      CIJ[3][5]=CIJ[5][3]=constants[19]; //C46
      CIJ[4][5]=CIJ[5][4]=constants[20]; //C56
      break;
    }
    default:{
      std::cout << "\nelasticityModels: Supported models in 3D - ISOTROPIC/TRANSVERSE/ORTHOTROPIC/ANISOTROPIC\n"; 
      std::cout << "See /src/elasticityModels.h\n";       
      exit(-1);
    }
    }
    break;
  }
  default:{
    std::cout << "\nelasticityModels: DIM is not 1/2/3\n"; 
    exit(-1);
  }
  }
  //print CIJ to terminal
  pcout << "Elasticity matrix (Voigt notation):\n";
  char buffer[100];
  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
    for (unsigned int j=0; j<2*dim-1+dim/3; j++){
      sprintf(buffer, "%8.3e ", CIJ[i][j]);
      pcout << buffer;
    }
    pcout << "\n";
  }
  pcout << "\n";
}

template <int dim>
void getCIJMatrix(elasticityModel model, std::vector<double> constants, dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> >& CIJ, dealii::ConditionalOStream& pcout){
  //CIJ.fill(0.0);
  pcout << "Reading material model:";
  switch (dim){
  case 1:{
    pcout << " 1D ";
    //1D models
    switch (model){
    case ISOTROPIC:{
      pcout << " ISOTROPIC \n";
      CIJ[0][0]=constants[0];
      break;
    }
    default:{
      std::cout << "\nelasticityModels: Supported models in 1D - ISOTROPIC\n";
      std::cout << "See /src/elasticityModels.h\n";
      exit(-1);
    }
    }
    break;
  }
  case 2:{
   pcout << " 2D ";
    //2D models
    switch (model){
    case ISOTROPIC:{
      pcout << " ISOTROPIC \n";
      double E=constants[0], nu=constants[1];
      double  mu=E/(2*(1+nu)), lambda= nu*E/((1+nu)*(1-2*nu));
      CIJ[0][0]=lambda+2*mu;
      CIJ[1][1]=lambda+2*mu;
      CIJ[2][2]=mu;
      CIJ[0][1]=CIJ[1][0]=lambda;
      break;
    }
    case ANISOTROPIC:{
      pcout << " ANISOTROPIC \n";
      CIJ[0][0]=constants[0]; //C11
      CIJ[1][1]=constants[1]; //C22
      CIJ[2][2]=constants[2]; //C33
      CIJ[0][1]=CIJ[1][0]=constants[3]; //C12
      CIJ[0][2]=CIJ[2][0]=constants[4]; //C13
      CIJ[1][2]=CIJ[2][1]=constants[5]; //C23
      break;
    }
    default:{
      std::cout << "\nelasticityModels: Supported models in 2D - ISOTROPIC/ANISOTROPIC\n";
      std::cout << "See /src/elasticityModels.h\n";
      exit(-1);
    }
    }
    break;
  }
  case 3:{
    pcout << " 3D ";
    //3D models
    switch (model){
    case ISOTROPIC:{
      pcout << " ISOTROPIC \n";
      double E=constants[0], nu=constants[1];
      double  mu=E/(2*(1+nu)), lambda= nu*E/((1+nu)*(1-2*nu));
      CIJ[0][0]=lambda+2*mu;
      CIJ[1][1]=lambda+2*mu;
      CIJ[2][2]=lambda+2*mu;
      CIJ[3][3]=mu;
      CIJ[4][4]=mu;
      CIJ[5][5]=mu;
      CIJ[0][1]=CIJ[1][0]=lambda;
      CIJ[0][2]=CIJ[2][0]=lambda;
      CIJ[1][2]=CIJ[2][1]=lambda;
      break;
    }
    case TRANSVERSE:{
      pcout << " TRANSVERSE \n";
      CIJ[0][0]=constants[0]; //C11
      CIJ[1][1]=constants[0]; //C11
      CIJ[2][2]=constants[1]; //C33
      CIJ[3][3]=constants[2]; //C44
      CIJ[4][4]=constants[2]; //C44
      CIJ[5][5]=(constants[0]-constants[3])/2.0; //(C11-C12)/2
      CIJ[0][1]=CIJ[1][0]=constants[3]; //C12
      CIJ[0][2]=CIJ[2][0]=constants[4]; //C13
      CIJ[1][2]=CIJ[2][1]=constants[4]; //C13
      break;
    }
    case ORTHOTROPIC:{
      pcout << " ORTHOTROPIC \n";
      CIJ[0][0]=constants[0]; //C11
      CIJ[1][1]=constants[1]; //C22
      CIJ[2][2]=constants[2]; //C33
      CIJ[3][3]=constants[3]; //C44
      CIJ[4][4]=constants[4]; //C55
      CIJ[5][5]=constants[5]; //C66
      CIJ[0][1]=CIJ[1][0]=constants[6]; //C12
      CIJ[0][2]=CIJ[2][0]=constants[7]; //C13
      CIJ[1][2]=CIJ[2][1]=constants[8]; //C23
      break;
    }
    case ANISOTROPIC:{
      pcout << " ANISOTROPIC \n";
      CIJ[0][0]=constants[0]; //C11
      CIJ[1][1]=constants[1]; //C22
      CIJ[2][2]=constants[2]; //C33
      CIJ[3][3]=constants[3]; //C44
      CIJ[4][4]=constants[4]; //C55
      CIJ[5][5]=constants[5]; //C66
      CIJ[0][1]=CIJ[1][0]=constants[6]; //C12
      CIJ[0][2]=CIJ[2][0]=constants[7]; //C13
      CIJ[0][3]=CIJ[3][0]=constants[8]; //C14
      CIJ[0][4]=CIJ[4][0]=constants[9]; //C15
      CIJ[0][5]=CIJ[5][0]=constants[10]; //C16
      CIJ[1][2]=CIJ[2][1]=constants[11]; //C23
      CIJ[1][3]=CIJ[3][1]=constants[12]; //C24
      CIJ[1][4]=CIJ[4][1]=constants[13]; //C25
      CIJ[1][5]=CIJ[5][1]=constants[14]; //C26
      CIJ[2][3]=CIJ[3][2]=constants[15]; //C34
      CIJ[2][4]=CIJ[4][2]=constants[16]; //C35
      CIJ[2][5]=CIJ[5][2]=constants[17]; //C36
      CIJ[3][4]=CIJ[4][3]=constants[18]; //C45
      CIJ[3][5]=CIJ[5][3]=constants[19]; //C46
      CIJ[4][5]=CIJ[5][4]=constants[20]; //C56
      break;
    }
    default:{
      std::cout << "\nelasticityModels: Supported models in 3D - ISOTROPIC/TRANSVERSE/ORTHOTROPIC/ANISOTROPIC\n";
      std::cout << "See /src/elasticityModels.h\n";
      exit(-1);
    }
    }
    break;
  }
  default:{
    std::cout << "\nelasticityModels: DIM is not 1/2/3\n";
    exit(-1);
  }
  }
  //print CIJ to terminal
  pcout << "Elasticity matrix (Voigt notation):\n";
  char buffer[100];
  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
    for (unsigned int j=0; j<2*dim-1+dim/3; j++){
      sprintf(buffer, "%8.3e ", CIJ[i][j][0]);
      pcout << buffer;
    }
    pcout << "\n";
  }
  pcout << "\n";
}

#endif
