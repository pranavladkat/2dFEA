#ifndef MESH_HPP
#define MESH_HPP

#include <iostream>
#include <vector>

using namespace std;


class Node{
  friend class Mesh;
private:
  int NodeID;
  double x,y,z;
};


class Face{
  friend class Mesh;
private:
  typedef enum {TRI, QUAD} FaceType;
  FaceType Ftype;
  int FaceID;
  int N1, N2, N3, N4;
};


class Boundary{
  friend class Mesh;
private:
  typedef enum {NODE, ELEMENT} BoundaryType;

  string name;
  BoundaryType BType;
  int size;
  int* BoundaryElement;

};


class Mesh{



};

#endif // MESH_HPP
