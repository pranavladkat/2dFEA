#include <iostream>
#include "mesh.hpp"
#include "material.hpp"
#include "preprocessor.hpp"
#include "solver.hpp"

using namespace std;

int main(int argc, char* argv[]){

  PetscInitialize(&argc,&argv,(char*)0,NULL);

  Mesh mesh("4x4Quad.dat");
  mesh.ReadMeshFile();
  mesh.Set_Thickness(0.1);
  mesh.ValidateMesh();
  mesh.WriteMesh(Mesh::MATLAB);

  Material steel(3.0E+7,0.3);
  steel.Compute_Elastic_Stiffness();
  steel.Print_Elastic_Stiffness();

  PreProcessor pre(&mesh,&steel);
  pre.Set_quadrature_rule(Q2D_2point);
  pre.Create_Quadrature_Objects();

  pre.Compute_Element_properties();
  pre.Compute_Element_stiffness();
  pre.Assemble_Stiffness_Matrix();
  pre.set_pointload(-1000.0);
  pre.Apply_BC();

  FEA_Solver solver(&pre);
  solver.solve_disp();
  solver.write_sol_disp();

  cout << "Program Finished!" << endl;

  // call destructor to free PETSc objects before PetscFinalize()
  pre.~PreProcessor();
  solver.~FEA_Solver();
  PetscFinalize();

  return 0;
}
