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
//  virtual void Compute_Jacobian(double const&,double const&) = 0;
//  virtual void Compute_shape_function(double const&,double const&) = 0;
//  virtual void Compute_dx_dxi(double const&,double const&) = 0;
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
//  virtual void Compute_Jacobian(double const&,double const&);
//  virtual void Compute_shape_function(double const&,double const&);
//  virtual void Compute_dx_dxi(double const&,double const&);
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
}

Quad4::Quad4(Quadrature const* quad, const vector<Node>& n, const Face& f)
  :Quad(quad), node(n), face(f)
{
  alpha = NULL;
  beta = NULL;
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
}



void Quad4 :: Compute_mapping_coeff(){

  int n = 4, nrhs = 2;
  int lda = 4, ldb = 4;
  int info, ipiv[4];
  double solution[2][4];
  double **mat = Quad->QMapping();
  alpha = new double [4];
  beta = new double [4];

  // get x and y coordinates of face nodes
  for(int i = 0; i < 4; i++){
    solution[0][i] = node[face.nodes[i]-1].x;
    solution[1][i] = node[face.nodes[i]-1].y;
  }

  /* solve AX = B using lapack solver */
  dgesv_(&n, &nrhs, &mat[0][0], &lda, ipiv, &solution[0][0], &ldb, &info);
  assert(info == 0);

  cout << "Solution: "<< info << endl;
  for(int i = 0; i < 4; i++){
    alpha[i] = solution[0][i];
    beta[i] = solution[1][i];
    cout << setw(10) << alpha[i] << setw(10) << beta[i] << endl;
  }

}




#endif // ELEMENT_HPP
