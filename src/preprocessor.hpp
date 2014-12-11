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
  void Compute_RHS();

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

  /* initialize K matrix */
  MatCreate(PETSC_COMM_WORLD,&KMat);
  MatSetSizes(KMat,PETSC_DECIDE,PETSC_DECIDE,GDof,GDof);
  MatSetFromOptions(KMat);
  MatSetUp(KMat);
  MatZeroEntries(KMat);
  MatAssemblyBegin(KMat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(KMat,MAT_FINAL_ASSEMBLY);

  for(size_t e = 0; e < element.size(); e++){
    const int*  P = stiffness[e]->Get_P();
    double** K = stiffness[e]->Get_K();
    int K_size = stiffness[e]->Get_K_size();

    for(int z = 0; z < K_size; z++){
      MatSetValues(KMat,1,&P[z],K_size,P,K[z],ADD_VALUES);
    }
  }

  MatAssemblyBegin(KMat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(KMat,MAT_FINAL_ASSEMBLY);

  //WriteMat(KMat,"KMat");

}


void PreProcessor :: Compute_RHS(){

  int GDof = 2*mesh->node.size();

  VecCreate(PETSC_COMM_WORLD,&RHS);
  VecSetSizes(RHS,PETSC_DECIDE,GDof);
  VecSetFromOptions(RHS);
  VecSet(RHS,0.0);
  VecDuplicate(RHS,&Solution);

  PetscReal *_RHS;
  VecGetArray(RHS,&_RHS);

  int BCnode = mesh->boundary[1].nodes[0];
  int BCP = (BCnode-1)*2 + 1;
  _RHS[BCP] = -1000.0;

  VecRestoreArray(RHS,&_RHS);

  int bottom_node[mesh->boundary[0].nodes.size()];
  for(size_t i = 0; i < mesh->boundary[0].nodes.size(); i++){
    bottom_node[i] = mesh->boundary[0].nodes[i];
  }

  int rows[2*mesh->boundary[0].nodes.size()];
  for(size_t i = 0; i < mesh->boundary[0].nodes.size(); i++){
    rows[i*2] = (bottom_node[i]-1)*2;
    rows[i*2+1] = rows[i*2] + 1;
    //cout << bottom_node[i] << "\t" << rows[i*2] << "\t" << rows[i*2+1] << endl;
  }


  MatZeroRows(KMat,2*mesh->boundary[0].nodes.size(),rows,1.0,NULL,NULL);

  WriteMat(KMat,"KMat");
  WriteVec(RHS,"RHS");

  KSP ksp;
  int itn;
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,KMat,KMat);
  KSPSetTolerances(ksp,1.e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
  KSPSolve(ksp,RHS,Solution);
  KSPGetIterationNumber(ksp,&itn);
  PetscPrintf(PETSC_COMM_WORLD,"KSP: %d\n",itn);

  ofstream sol("sol.dat");
  PetscReal *_sol;
  VecGetArray(Solution,&_sol);
  for(int i = 0; i < GDof; i+=2){
    sol << sqrt(pow(_sol[i],2)+pow(_sol[i+1],2)) << endl;
  }
  sol.close();
  VecRestoreArray(Solution,&_sol);



  WriteVec(Solution,"solution");

}




#endif // PREPROCESSOR_HPP
