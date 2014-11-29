#ifndef STIFFELEMENT_HPP
#define STIFFELEMENT_HPP

#include <iostream>
#include "element.hpp"
#include "mesh.hpp"

using namespace std;


class EStiffness{
private:

  const Element* element;
  const size_t K_size;            // stiffness matrix dimension
  double **K, *K_data;            // element stiffness matrix
  double *p, *q;                  // global equation number
public:
  EStiffness(const Element*);
  ~EStiffness();

};



// Functions


EStiffness ::EStiffness(const Element* e)
  : element(e), K_size(2*e->Quad->Qpoints())
{
  K = NULL;
  K_data = NULL;
  p = NULL;
  q = NULL;
}


EStiffness :: ~EStiffness(){
  if(K != NULL)
    delete [] K;
  if(K_data != NULL)
    delete [] K;
  if(p != NULL)
    delete [] p;
  if(q != NULL)
    delete [] q;
}









#endif // STIFFELEMENT_HPP
