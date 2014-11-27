#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <iostream>
#include "quadrature.hpp"

using namespace std;


class Element{

protected:
  Quadrature *Quad;
  double *alpha,*beta;
  double *J;
  double *dx_dxi, *dx_deta;
  double *dy_dxi, *dy_deta;
  double *dxi_dx, *dxi_dy;
  double *deta_dx, *deta_dy;

public:

  Element();
  virtual ~Element();

  void Set_Quadrature(Quadrature*);

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

public:
  ~Quad4();
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

void Element :: Set_Quadrature(Quadrature* quad){
  Quad = quad;
}


Element :: ~Element(){
  if(J == NULL){
    cout << "J doesnt exist" << endl;
  }
}


void Quad4 :: Compute_mapping_coeff(){

}

Quad4::~Quad4(){

}





#endif // ELEMENT_HPP
