#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <iostream>
#include "preprocessor.hpp"

using namespace std;

class FEA_Solver{

private:
  const PreProcessor* prep;
  Vec Solution;
  KSP ksp;
public:
  FEA_Solver(const PreProcessor* pre)
    : prep(pre)
  {
    VecDuplicate(prep->RHS,&Solution);
  }

  void solve_disp(double tol = 1e-12){

    int itn;
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,prep->KMat,prep->KMat);
    KSPSetTolerances(ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    KSPSolve(ksp,prep->RHS,Solution);
    KSPGetIterationNumber(ksp,&itn);
    PetscPrintf(PETSC_COMM_WORLD,"Iterations taken by KSP: %d\n",itn);

  }

  void write_sol_disp(){

    ofstream disp_total("disp_total.dat");
    ofstream disp_u("disp_u.dat");
    ofstream disp_v("disp_v.dat");
    PetscReal *_sol;
    VecGetArray(Solution,&_sol);
    for(size_t i = 0; i < prep->GDof; i+=2){
      disp_total << sqrt(pow(_sol[i],2)+pow(_sol[i+1],2)) << endl;
      disp_u << _sol[i] << endl;
      disp_v << _sol[i+1] << endl;
    }
    disp_total.close();
    disp_u.close();
    disp_v.close();
    VecRestoreArray(Solution,&_sol);

    //WriteVec(Solution,"solution");

  }



  ~FEA_Solver(){
    VecDestroy(&Solution);
    KSPDestroy(&ksp);
  }

};

#endif // SOLVER_HPP
