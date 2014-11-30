#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <iostream>
#include <cassert>
#include <iomanip>

using namespace std;

class Material{
  friend class PreProcessor;
private:
  double E; // young's modulus
  double nu; // poisson's ratio

  double **Estiff, *Estiff_data;    // Element stiffness

public:

  Material();
  ~Material();
  Material(double const&,double const&);
  void set_YoungsModulus(double const&);
  void set_PoissonsRatio(double const&);
  void Compute_Elastic_Stiffness();
  void Print_Elastic_Stiffness();
  double** Get_Element_Stiffness() const {return Estiff;}
  void Allocate_Estiff();

};


/******************* Functions **************************/

Material :: Material(){
  E = 0.0;
  nu = 0.0;
  Allocate_Estiff();
}

Material :: Material(double const& YoungsMod, double const& PoissonRatio){
  set_YoungsModulus(YoungsMod);
  set_PoissonsRatio(PoissonRatio);
  Allocate_Estiff();
}

void Material :: Allocate_Estiff(){
  Estiff = new double* [3];
  Estiff_data = new double [9];
  for(int i = 0; i < 3; i++){
    Estiff[i] = &Estiff_data[i*3];
  }
}

void Material :: set_YoungsModulus(const double & YoungsMod){
  E = YoungsMod;
}

void Material :: set_PoissonsRatio(const double & PoissonRatio){
  nu = PoissonRatio;
}

void Material :: Compute_Elastic_Stiffness(){
  assert(E != 0 || nu != 0);

  /* Refer: Introduction to Finite Element Method, J. N. Reddy, pg. 609, Eq. 11.2.11  */
  Estiff[0][0] = E/(1-nu*nu);
  Estiff[0][1] = E*nu/(1-nu*nu);
  Estiff[0][2] = 0.0;
  Estiff[1][0] = E*nu/(1-nu*nu);
  Estiff[1][1] = E/(1-nu*nu);
  Estiff[1][2] = 0.0;
  Estiff[2][0] = 0.0;
  Estiff[2][1] = 0.0;
  Estiff[2][2] = E/(2*(1+nu));
}

void Material :: Print_Elastic_Stiffness(){
  cout << "Elastic Stiffness :"<< endl;
  cout << setw(15) << Estiff[0][0] << setw(15) << Estiff[0][1] << setw(15) << Estiff[0][2] << endl;
  cout << setw(15) << Estiff[1][0] << setw(15) << Estiff[1][1] << setw(15) << Estiff[1][2] << endl;
  cout << setw(15) << Estiff[2][0] << setw(15) << Estiff[2][1] << setw(15) << Estiff[2][2] << endl;
  cout << endl;
}



Material :: ~Material(){
  delete [] Estiff;
  delete [] Estiff_data;
}


#endif // MATERIAL_HPP
