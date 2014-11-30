#ifndef MESH_HPP
#define MESH_HPP

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <iomanip>

using namespace std;


/*
 * CLASS NODE -> contains mesh coordinate details
 */
class Node{
  friend class Mesh;
  friend class PreProcessor;
  friend class Quad4;
private:
  int NodeID;
  double x,y,z;
};


/*
 * CLASS FACE -> contains mesh faces and its connected nodes
 */
class Face{
  friend class Mesh;
  friend class PreProcessor;
  friend class Quad4;
private:
  typedef enum {TRI, QUAD} FaceType;
  FaceType Ftype;
  int FaceID;
  vector<int> nodes;
};


/*
 * CLASS BOUNDARY -> identifies boundary nodes
 */
class Boundary{
  friend class Mesh;
  friend class PreProcessor;
private:
  typedef enum {NODE, ELEMENT} BoundaryType;

  string name;
  BoundaryType BType;
  vector<int> nodes;
};


/*
 * CLASS MESH -> Reads the mesh file and populates mesh data
 */
class Mesh{
  friend class PreProcessor;
private:
  vector<Node> node;
  vector<Face> face;
  vector<Boundary> boundary;
  bool set_filename;
  string filename;
  bool isQuadPresent, isTriPresent;
  double thickness;

public:

  typedef enum {MATLAB,PLOT3D} OUTPUT_MESH_FORMAT;
  Mesh();
  Mesh(string const&);
  void SetMeshFilename(string const&);
  void ReadMeshFile();
  void ValidateMesh();
  void WriteMesh(OUTPUT_MESH_FORMAT const&);
  void Set_Thickness(double const&);
  double Get_Thickness() const ;
};


/************* Function Definitions for Class Mesh ******************/

Mesh::Mesh(){
  set_filename = false;
  isQuadPresent = false;
  isTriPresent = false;
}

Mesh::Mesh(const string &a){
  SetMeshFilename(a);
  isQuadPresent = false;
  isTriPresent = false;
}

void Mesh::SetMeshFilename(const string &a){
  filename = a;
  set_filename = true;
  //cout << "Mesh Filename = " << filename << endl;
}


void Mesh::ReadMeshFile(){

  /* assert if input filename is set */
  assert(set_filename);

  ifstream mfile(filename);
  assert(mfile.is_open());

  string parameter;

  while(!mfile.eof()){

    mfile >> parameter;
    //cout << parameter << endl;
    if(parameter.compare("#Nodes") == 0){
      for(;;){
        Node read_node;
        mfile >> read_node.NodeID;
        if(read_node.NodeID == -1){
          break;
        }
        mfile >> read_node.x >> read_node.y >> read_node.z;
        node.push_back(read_node);
        //cout << setw(12) << read_node.x << setw(12) << read_node.y << setw(12) << read_node.z << endl;
      }
    }

    if(parameter.compare("#Elements") == 0){

      for(;;){
        Face read_face;
        int check, temp;
        mfile >> check;
        if(check == -1){
          break;
        }
        mfile >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp;
        mfile >> read_face.FaceID;
        mfile >> temp; read_face.nodes.push_back(temp);
        mfile >> temp; read_face.nodes.push_back(temp);
        mfile >> temp; read_face.nodes.push_back(temp);
        mfile >> temp; read_face.nodes.push_back(temp);

        assert(read_face.nodes[2] != read_face.nodes[3]);
        read_face.Ftype = Face::QUAD;
        if(read_face.Ftype == Face::QUAD && !isQuadPresent){
          isQuadPresent = true;
          cout << "Quad Face is present" << endl;
        }
        if(read_face.Ftype == Face::TRI && !isTriPresent){
          isTriPresent = true;
          cout << "Tri Face is present" << endl;
        }

        face.push_back(read_face);

        //cout << setw(12) << read_face.FaceID << setw(12) << read_face.nodes[0] << setw(12) << read_face.nodes[1]
        //     << setw(12) << read_face.nodes[2] << setw(12) << read_face.nodes[3] << endl;
      }
    }

    if(parameter.compare("#NamedSelection") == 0){

      double size;
      Boundary read_b;
      mfile >> read_b.name;
      string btype;
      mfile >> btype;

      if(btype.compare("NODE")==0){
        read_b.BType = Boundary::NODE;
      }else{
        cerr << "ERROR in Boundary "<< read_b.name << endl;
      }

      mfile >> size;
      for(int i = 0; i < size; i++){
        int buf;
        mfile >> buf;
        read_b.nodes.push_back(buf);
        //cout << read_b.nodes[i] << "\t";
      }
      //cout << endl;
      boundary.push_back(read_b);
    }

    if(parameter.compare("#End") == 0 || parameter.compare("#end") == 0){
      break;
    }

  } // end while
  mfile.close();
} // end mesh read function


void Mesh::ValidateMesh(){
  cout << "Number of nodes = " << node.size() << endl;
  cout << "Number of faces = " << face.size() << endl;
  cout << "Number of boundaries = " << boundary.size() << endl;
  cout << "size of nodes: " << sizeof(Node)*node.size() << " bytes" << endl;
  cout << "size of faces: " << sizeof(Face)*face.size() << " bytes" << endl;
  cout << endl;
}


void Mesh::WriteMesh(OUTPUT_MESH_FORMAT const& output_mesh_format){

  if(output_mesh_format == MATLAB){

    ofstream mfile("mesh.m");
    assert(mfile.is_open());

    // write x component
    mfile << "x = [" << endl;
    for(size_t i = 0; i < node.size(); i++){
      mfile << node[i].x << endl;
    }
    mfile << "];" << endl;
    // write y component
    mfile << "y = [" << endl;
    for(size_t i = 0; i < node.size(); i++){
      mfile << node[i].y << endl;
    }
    mfile << "];" << endl;
    // write z component
    mfile << "z = [" << endl;
    for(size_t i = 0; i < node.size(); i++){
      mfile << node[i].z << endl;
    }
    mfile << "];" << endl;

  } // end matlab mesh format write

}  // end writemesh function


void Mesh :: Set_Thickness(double const& Thickness){
  thickness = Thickness;
}

double Mesh :: Get_Thickness() const {
  return thickness;
}


#endif // MESH_HPP
