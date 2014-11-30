#ifndef STIFFELEMENT_HPP
#define STIFFELEMENT_HPP

#include <iostream>
#include "element.hpp"
#include "mesh.hpp"
#include "material.hpp"
#include "cblas.h"

using namespace std;


class EStiffness{
private:
  const Element* element;
  const Material* material;
  const size_t K_size;            // stiffness matrix dimension
  double **K, *K_data;            // element stiffness matrix
  double *P, *Q;                  // global equation number
public:
  EStiffness(const Material*,const Element*);
  ~EStiffness();
  void Compute_Element_Stiffness(double const&);
  void Compute_Equation_Number(const vector<int>&);
};



// Functions


EStiffness :: EStiffness(const Material* m,const Element* e)
  : element(e),material(m), K_size(2*e->Quad->Qpoints())
{
  assert(K_size != 0);
  P = new double [K_size];
  Q = new double [K_size];
  K = new double* [K_size];
  K_data = new double [K_size*K_size];
  for(size_t i = 0; i < K_size; i++){
    K[i] = &K_data[i*K_size];
  }
}


EStiffness :: ~EStiffness(){
  delete [] K_data;
  delete [] K;
  delete [] P;
  delete [] Q;
}



void EStiffness :: Compute_Element_Stiffness(double const& thickness){
  double B[3][K_size];
  double** C = material->Get_Element_Stiffness();
  double* QW = element->Quad->QWeights();
  double result[8][3];
  double beta[4] = {0.0,1.0,1.0,1.0};

  assert(K_size == 8);
  for(size_t z = 0; z < K_size/2; z++){
    // compute B
    B[0][0] = element->dN1_dxi[z]*element->dxi_dx[z] + element->dN1_deta[z]*element->deta_dx[z];
    B[0][1] = 0.0;
    B[1][0] = 0.0;
    B[1][1] = element->dN1_dxi[z]*element->dxi_dy[z] + element->dN1_deta[z]*element->deta_dy[z];
    B[2][0] = B[1][1];
    B[2][1] = B[0][0];

    B[0][2] = element->dN2_dxi[z]*element->dxi_dx[z] + element->dN2_deta[z]*element->deta_dx[z];
    B[0][3] = 0.0;
    B[1][2] = 0.0;
    B[1][3] = element->dN2_dxi[z]*element->dxi_dy[z] + element->dN2_deta[z]*element->deta_dy[z];
    B[2][2] = B[1][3];
    B[2][3] = B[0][2];

    B[0][4] = element->dN3_dxi[z]*element->dxi_dx[z] + element->dN3_deta[z]*element->deta_dx[z];
    B[0][5] = 0.0;
    B[1][4] = 0.0;
    B[1][5] = element->dN3_dxi[z]*element->dxi_dy[z] + element->dN3_deta[z]*element->deta_dy[z];
    B[2][4] = B[1][5];
    B[2][5] = B[0][4];

    B[0][6] = element->dN4_dxi[z]*element->dxi_dx[z] + element->dN4_deta[z]*element->deta_dx[z];
    B[0][7] = 0.0;
    B[1][6] = 0.0;
    B[1][7] = element->dN4_dxi[z]*element->dxi_dy[z] + element->dN4_deta[z]*element->deta_dy[z];
    B[2][6] = B[1][7];
    B[2][7] = B[0][6];

    //  for(int i = 0; i < 3; i++){
    //    for(int j = 0; j < 8; j++){
    //      cout << setw(12) << B[i][j];
    //    }
    //    cout << endl;
    //  }
    //  cout << endl;

    double alpha = element->J[z]*thickness*QW[z];
    // matrix multiplication --> result = B(transpose)*C*alpha
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,K_size,3,3,alpha,&B[0][0],K_size,&C[0][0],3,0.0,&result[0][0],3);

    //  for(int i = 0; i < 8; i++){
    //    for(int j = 0; j < 3; j++){
    //      cout << setw(15) << result[i][j];
    //    }
    //    cout << endl;
    //  }

    // matrix multiplication --> K = result*B
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,K_size,K_size,3,1.0,&result[0][0],3,&B[0][0],K_size,beta[z],&K[0][0],K_size);
  }

//  for(int i = 0; i < 8; i++){
//    for(int j = 0; j < 8; j++){
//      cout << setw(15) << K[i][j];
//    }
//    cout << endl;
//  }
//  cout << endl;

}


void EStiffness :: Compute_Equation_Number(const vector<int>& node){
  for(size_t i = 0; i < K_size/2; i++){
    P[i*2] = (node[i]-1)*2;
    Q[i*2] = P[i*2];
    P[i*2+1] = P[i*2] + 1;
    Q[i*2+1] = Q[i*2] + 1;
  }

//  for(size_t i = 0; i < K_size; i++){
//    cout << setw(12) << P[i];
//  }
//  cout << endl;

}





#endif // STIFFELEMENT_HPP
