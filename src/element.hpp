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
  friend class EStiffness;
protected:
  const Quadrature *Quad;         // quadrature info
  const Face& face;               // face associated with element
  double *alpha,*beta;            // mapping coeff to master element
  double **Na, *Na_data;          // shape function evaluated at quadrature points
  double *J;                      // jacobian evaluated at quadrature points
  double *dx_dxi, *dx_deta;       // derivative of x wrt xi and eta evaluated at quadrature points
  double *dy_dxi, *dy_deta;       // derivative of y wrt xi and eta evaluated at quadrature points
  double *dxi_dx, *dxi_dy;        // derivative of xi wrt x and y evaluated at quadrature points
  double *deta_dx, *deta_dy;      // derivative of eta wrt x and y evaluated at quadrature points

public:

  Element(Quadrature const* quad, const Face& f);
  virtual ~Element();

  virtual void Element_setup(vector<Node> const&) = 0;
  virtual void Compute_mapping_coeff(vector<Node> const&) = 0;
  virtual void Compute_Shape_Function() = 0;
  virtual void Compute_dX_dXI() = 0;
  virtual void Compute_Jacobian() = 0;
  virtual void Compute_dXI_dX() = 0;

};



class Quad4 : public Element{
  friend class EStiffness;
protected:

public:
  Quad4(Quadrature const*, const Face&);
  ~Quad4();
  virtual void Element_setup(vector<Node> const&);
  virtual void Compute_mapping_coeff(vector<Node> const&);
  virtual void Compute_Shape_Function();
  virtual void Compute_dX_dXI();
  virtual void Compute_Jacobian();
  virtual void Compute_dXI_dX();

};





// Functions

Element::Element(Quadrature const* quad, const Face& f)
  :Quad(quad), face(f)
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


//free memory
Element :: ~Element(){
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
  if(dx_dxi != NULL){
    delete [] dx_dxi;
  }
  if(dx_deta != NULL){
    delete [] dx_deta;
  }
  if(dy_dxi != NULL){
    delete [] dy_dxi;
  }
  if(dy_deta != NULL){
    delete [] dy_deta;
  }
  if(dxi_dx != NULL){
    delete [] dxi_dx;
  }
  if(dxi_dy != NULL){
    delete [] dxi_dy;
  }
  if(deta_dx != NULL){
    delete [] deta_dx;
  }
  if(deta_dy != NULL){
    delete [] deta_dy;
  }
  if(J != NULL){
    delete [] J;
  }
}

Quad4::Quad4(Quadrature const* quad, const Face& f)
  :Element(quad,f)
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



void Quad4 :: Element_setup(vector<Node> const& node){
  assert(Quad != NULL);
  Compute_mapping_coeff(node);
  Compute_Shape_Function();
  Compute_dX_dXI();
  Compute_Jacobian();
  Compute_dXI_dX();
}



void Quad4 :: Compute_mapping_coeff(vector<Node> const& node){

  int n = Quad->Qpoints(), nrhs = 2;
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

  cout << "Solution: "<< info << endl;
  for(int i = 0; i < n; i++){
    alpha[i] = solution[0][i];
    beta[i] = solution[1][i];
    cout << setw(15) << alpha[i] << setw(15) << beta[i] << endl;
  }
  cout << endl;

}


void Quad4 :: Compute_Shape_Function(){

  const int n = Quad->Qpoints();
  const double* QXi  = Quad->QXipoints();
  const double* QEta = Quad->QEtapoints();
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

  cout << "Shape Functions:" << endl;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      cout << setw(12) << Na[i][j];
    }
    cout << endl;
  }
  cout << endl;

}


void Quad4 :: Compute_dX_dXI(){

  assert(alpha != NULL || beta != NULL);
  const int n = Quad->Qpoints();
  const double* QXi = Quad->QXipoints();
  const double* QEta = Quad->QEtapoints();
  dx_dxi = new double [n];
  dx_deta = new double [n];
  dy_dxi = new double [n];
  dy_deta = new double [n];

  cout << setw(15) << "dx_dxi" << setw(15) << "dx_deta" << setw(15) << "dy_dxi" << setw(15) << "dy_deta" << endl;
  for(int i = 0; i < n; i++){
    dx_dxi[i] = alpha[1] + alpha[3]*QEta[i];
    dx_deta[i] = alpha[2] + alpha[3]*QXi[i];
    dy_dxi[i] = beta[1] + beta[3]*QEta[i];
    dy_deta[i] = beta[2] + beta[3]*QXi[i];
    cout << setw(15) << dx_dxi[i] << setw(15) << dx_deta[i]
         << setw(15) << dy_dxi[i] << setw(15) << dy_deta[i] << endl;
  }
  cout << endl;
}


void Quad4 :: Compute_Jacobian(){
  assert(dx_dxi != NULL || dx_deta != NULL ||dy_dxi != NULL ||dy_deta != NULL);
  const int n = Quad->Qpoints();
  J = new double [n];

  cout << "Jacobian :" << endl;
  for(int i = 0; i < n; i++){
    J[i] = dx_dxi[i]*dy_deta[i] - dy_dxi[i]*dx_deta[i];
    cout << setw(15) << J[i];
  }
  cout << endl << endl;

  /* assert positive jacobian */
  assert(*J > 0);
}


void Quad4 :: Compute_dXI_dX(){
  assert(dx_dxi != NULL || dx_deta != NULL ||dy_dxi != NULL ||dy_deta != NULL);
  assert(*J < 1e8);
  const int n = 4;
  dxi_dx = new double [n];
  dxi_dy = new double [n];
  deta_dx = new double [n];
  deta_dy = new double [n];

  cout << setw(15) << "dxi_dx" << setw(15) << "dxi_dy" << setw(15) << "deta_dx" << setw(15) << "deta_dy" << endl;
  for(int i = 0; i < n; i++){
    dxi_dx[i]  =  dy_deta[i]/J[i];
    dxi_dy[i]  = -dx_deta[i]/J[i];
    deta_dx[i] = -dy_dxi[i]/J[i];
    deta_dy[i] =  dx_dxi[i]/J[i];
    cout << setw(15) << dxi_dx[i] << setw(15) << dxi_dy[i]
         << setw(15) << deta_dx[i] << setw(15) << deta_dy[i] << endl;
  }

}



#endif // ELEMENT_HPP
