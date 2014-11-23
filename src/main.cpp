#include <iostream>
#include "mesh.hpp"

using namespace std;

int main(int argc, char* argv[]){

  Mesh mesh("L-mesh.dat");
  mesh.ReadMeshFile();


  cout << "Hello Pranav!" << endl;

  return 0;
}
