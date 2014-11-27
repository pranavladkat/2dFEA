#include <iostream>
#include "mesh.hpp"
#include "material.hpp"
#include "preprocessor.hpp"

using namespace std;

int main(int argc, char* argv[]){

  Mesh mesh("L-mesh.dat");
  mesh.ReadMeshFile();
  mesh.ValidateMesh();
  //mesh.WriteMesh(Mesh::MATLAB);

  Material steel(30E+6,0.25);
  steel.Compute_Elastic_Stiffness();
  steel.Print_Elastic_Stiffness();

  PreProcessor pre(&mesh,&steel);
  pre.Set_quadrature_rule(Q2D_2point);
  pre.Create_Quadrature_Objects();

  pre.temp();

  cout << "Hello Pranav!" << endl;

  return 0;
}
