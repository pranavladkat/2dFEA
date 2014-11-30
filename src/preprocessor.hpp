#ifndef PREPROCESSOR_HPP
#define PREPROCESSOR_HPP

#include <iostream>
#include <vector>
#include "petscksp.h"
#include "mesh.hpp"
#include "material.hpp"
#include "element.hpp"
#include "quadrature.hpp"
#include "stiffelement.hpp"
#include "functions.h"

using namespace std;

class PreProcessor{

private:
  const Mesh *mesh;
  const Material *material;
  Quadrature_Rule QRule;
  Quadrature *Quad_Quad, *Quad_Tri;
  vector<Element*> element;
  vector<EStiffness*> stiffness;
  Mat KMat;
  Vec RHS, Solution;

public:

  PreProcessor(Mesh const*, Material const*);
  ~PreProcessor();

  void Set_quadrature_rule(Quadrature_Rule const &);
  void Create_Quadrature_Objects();
  void Compute_Element_properties();
  void Compute_Element_stiffness();
  void Assemble_Stiffness_Matrix();

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
    delete stiffness[i];
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



void PreProcessor :: Compute_Element_properties(){

  Element *elem;
  for(size_t i = 0; i < mesh->face.size(); i++){
    if(mesh->face[i].Ftype == Face::QUAD){
      elem = new Quad4(Quad_Quad,mesh->face[i]);
    }
    elem->Element_setup(mesh->node);
    element.push_back(elem);
  }

}



void PreProcessor :: Compute_Element_stiffness(){
  assert(mesh->Get_Thickness() != 0);
  EStiffness *estiff;
  for(size_t i = 0; i < mesh->face.size(); i++){
    estiff= new EStiffness(material,element[i]);
    estiff->Compute_Element_Stiffness(mesh->Get_Thickness());
    estiff->Compute_Equation_Number(mesh->face[i].nodes);
    stiffness.push_back(estiff);
  }
}



void PreProcessor :: Assemble_Stiffness_Matrix(){
  assert(stiffness.size() != 0);
  int GDof = 2*mesh->node.size();
  cout << "GDof = " << GDof << endl;

  /* initialize K matrix */
  MatCreate(PETSC_COMM_WORLD,&KMat);
  MatSetSizes(KMat,PETSC_DECIDE,PETSC_DECIDE,GDof,GDof);
  MatSetFromOptions(KMat);
  MatSetUp(KMat);
  MatZeroEntries(KMat);
  MatAssemblyBegin(KMat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(KMat,MAT_FINAL_ASSEMBLY);

  for(size_t e = 0; e < element.size(); e++){

    int*  P = stiffness[e]->Get_P();
    int*  Q = stiffness[e]->Get_Q();
    double** K = stiffness[e]->Get_K();
    int K_size = stiffness[e]->Get_K_size();

    for(int z = 0; z < K_size; z++){
      MatSetValues(KMat,1,&Q[z],K_size,P,K[z],ADD_VALUES);
    }

  }

  MatAssemblyBegin(KMat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(KMat,MAT_FINAL_ASSEMBLY);

  //WriteMat(KMat,"KMat");

}




#endif // PREPROCESSOR_HPP
