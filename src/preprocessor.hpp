#ifndef PREPROCESSOR_HPP
#define PREPROCESSOR_HPP

#include <iostream>
#include <vector>
#include "mesh.hpp"
#include "material.hpp"
#include "element.hpp"
#include "quadrature.hpp"

using namespace std;

class PreProcessor{

private:
  const Mesh *mesh;
  const Material *material;
  Quadrature *Quad_Quad, *Quad_Tri;
  vector<Element*> element;
  Element *elem;

  Quadrature_Rule QRule;

public:

  PreProcessor(Mesh const*, Material const*);
  ~PreProcessor();

  void Set_quadrature_rule(Quadrature_Rule const &);
  void Create_Quadrature_Objects();

  void temp();


};


/********************* functions ************************/

PreProcessor::PreProcessor(Mesh const* msh, Material const* matrl)
  :mesh(msh), material(matrl){
  Quad_Quad = NULL;
  Quad_Tri = NULL;
}


void PreProcessor :: Set_quadrature_rule(const Quadrature_Rule &qrule){
  QRule = qrule;
}

PreProcessor :: ~PreProcessor(){
  if(Quad_Quad != NULL){
    delete Quad_Quad;
  }
  if(Quad_Tri != NULL){
    delete Quad_Tri;
  }
  for(size_t i = 0; i < element.size(); i++){
    delete element[i];
  }
}


void PreProcessor :: Create_Quadrature_Objects(){
  if(QRule == Q2D_2point && mesh->isQuadPresent){
    Quad_Quad = new Quadrature_2PQuad4;
    Quad_Quad->Setup_Quadrature();
    Quad_Quad->Print_Quadrature_Info();
  }
  if(QRule == Q2D_3point && mesh->isQuadPresent){
    Quad_Quad = new Quadrature_3PQuad4;
    Quad_Quad->Setup_Quadrature();
    Quad_Quad->Print_Quadrature_Info();
  }
}



void PreProcessor :: temp(){
  cout << element.size() << endl;
  elem = new Quad4;
  element.push_back(elem);
  cout << element.size() << endl;

}



#endif // PREPROCESSOR_HPP
