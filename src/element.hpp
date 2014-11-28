#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <iostream>
#include <iomanip>
#include "quadrature.hpp"
#include "mesh.hpp"

extern "C" {
void dgesv_(int *n, int *nrhs,  double *a,  int  *lda,
            int *ipivot, double *b, int *ldb, int *info);
}

using namespace std;


class Element{

protected:
  double *alpha,*beta;
  double **Na, *Na_data;
  double *J;
  double *dx_dxi, *dx_deta;
  double *dy_dxi, *dy_deta;
  double *dxi_dx, *dxi_dy;
  double *deta_dx, *deta_dy;

public:

  Element();

  virtual ~Element();

  virtual void Element_setup() = 0;
  virtual void Compute_mapping_coeff() = 0;
  virtual void Compute_Shape_Function() = 0;
  //virtual void Compute_Jacobian() = 0;
  virtual void Compute_dx_dxi() = 0;
//  virtual void Compute_dx_deta(double const&,double const&) = 0;
//  virtual void Compute_dy_dxi(double const&,double const&) = 0;
//  virtual void Compute_dy_deta(double const&,double const&) = 0;
//  virtual void Compute_dxi_dx(double const&,double const&) = 0;
//  virtual void Compute_dxi_dy(double const&,double const&) = 0;
//  virtual void Compute_deta_dx(double const&,double const&) = 0;
//  virtual void Compute_deta_dy(double const&,double const&) = 0;

};



class Quad4 : public Element{

protected:
  const Quadrature *Quad;
  const vector<Node> node;
  const Face face;
public:
  Quad4(Quadrature const*,const vector<Node>&, const Face&);
  ~Quad4();
  virtual void Element_setup();
  virtual void Compute_mapping_coeff();
  virtual void Compute_Shape_Function();
  //virtual void Compute_Jacobian();
  virtual void Compute_dx_dxi();
//  virtual void Compute_dx_deta(double const&,double const&);
//  virtual void Compute_dy_dxi(double const&,double const&);
//  virtual void Compute_dy_deta(double const&,double const&);
//  virtual void Compute_dxi_dx(double const&,double const&);
//  virtual void Compute_dxi_dy(double const&,double const&);
//  virtual void Compute_deta_dx(double const&,double const&);
//  virtual void Compute_deta_dy(double const&,double const&);

};





// Functions

Element::Element(){
  alpha = NULL;
  beta = NULL;
  Na = NULL;
  Na_data = NULL;
  J = NULL;
  dx_dxi = NULL;
  dx_deta = NULL;
  dy_dxi = NULL;
  dy_deta = NULL;
  dxi_dx = NULL;
  dxi_dy = NULL;
  deta_dx = NULL;
  deta_dy = NULL;
}



Element :: ~Element(){
  //free memory
  if(alpha != NULL){
    delete [] alpha;
  }
  if(beta != NULL){
    delete [] beta;
  }
  if(Na != NULL){
    delete [] Na;
  }
  if(Na_data != NULL){
    delete [] Na_data;
  }
}

Quad4::Quad4(Quadrature const* quad, const vector<Node>& n, const Face& f)
  :Quad(quad), node(n), face(f)
{
  alpha = NULL;
  beta = NULL;
  Na = NULL;
  Na_data = NULL;
  J = NULL;
  dx_dxi = NULL;
  dx_deta = NULL;
  dy_dxi = NULL;
  dy_deta = NULL;
  dxi_dx = NULL;
  dxi_dy = NULL;
  deta_dx = NULL;
  deta_dy = NULL;
}


Quad4::~Quad4(){

}



void Quad4 ::Element_setup(){
  assert(Quad != NULL);
  Compute_mapping_coeff();
  Compute_Shape_Function();
}



void Quad4 :: Compute_mapping_coeff(){

  const int n = Quad->Qpoints(), nrhs = 2;
  int lda = n, ldb = n;
  int info, ipiv[n];
  double solution[2][n];
  double **mat = Quad->QMapping();
  alpha = new double [n];
  beta = new double [n];

  // get x and y coordinates of face nodes
  for(int i = 0; i < n; i++){
    solution[0][i] = node[face.nodes[i]-1].x;
    solution[1][i] = node[face.nodes[i]-1].y;
  }

  /* solve AX = B using lapack solver */
  dgesv_(&n, &nrhs, &mat[0][0], &lda, ipiv, &solution[0][0], &ldb, &info);
  assert(info == 0);

//  cout << "Solution: "<< info << endl;
//  for(int i = 0; i < n; i++){
//    alpha[i] = solution[0][i];
//    beta[i] = solution[1][i];
//    cout << setw(10) << alpha[i] << setw(10) << beta[i] << endl;
//  }

}


void Quad4 :: Compute_Shape_Function(){

  const int n = Quad->Qpoints();
  double* QXi  = Quad->QXipoints();
  double* QEta = Quad->QEtapoints();
  Na = new double* [n];
  Na_data = new double [n*n];
  for(int i = 0; i < n; i++){
    Na[i] = &Na_data[i*n];
  }

  for(int i = 0; i < n; i++){
    //compute shape functions at quadrature points
    Na[0][i] = 0.25*(1.0-QXi[i])*(1.0-QEta[i]);
    Na[1][i] = 0.25*(1.0+QXi[i])*(1.0-QEta[i]);
    Na[2][i] = 0.25*(1.0+QXi[i])*(1.0+QEta[i]);
    Na[3][i] = 0.25*(1.0-QXi[i])*(1.0+QEta[i]);
  }

//  cout << "Shape Functions:" << endl;
//  for(int i = 0; i < n; i++){
//    for(int j = 0; j < n; j++){
//      cout << setw(12) << Na[i][j];
//    }
//    cout << endl;
//  }

}


void Quad4 :: Compute_dx_dxi(){
  const int n = Quad->Qpoints();
  dx_dxi = new double [n];

}




#endif // ELEMENT_HPP
