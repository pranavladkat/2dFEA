#include <iostream>
#include "mesh.hpp"
#include "material.hpp"
#include "cblas.h"

using namespace std;

int main(int argc, char* argv[]){

  Mesh mesh("L-mesh.dat");
  mesh.ReadMeshFile();
  mesh.ValidateMesh();
  //mesh.WriteMesh(Mesh::MATLAB);

  Material steel(30E+6,0.25);
  steel.Compute_Elastic_Stiffness();
  steel.Print_Elastic_Stiffness();


  cout << "Hello Pranav!" << endl;

  return 0;
}
