#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

typedef enum {Q2D_2point, Q2D_3point} Quadrature_Rule;

class Quadrature{

protected:
  int Quadrature_points;            // number of quadrature points
  double *QW, *QXi, *QEta;          // quadrature weights and points
  double **mapping, *mapping_data;  // mapping matrix to compute mapping coeff for each element

public:
  virtual ~Quadrature();

  virtual void Setup_Quadrature() = 0;
  virtual void Print_Quadrature_Info() = 0;
  int Qpoints() const {return Quadrature_points;}
  double* QWeights() const {return QW;}
  double* QXipoints() const {return QXi;}
  double* QEtapoints() const {return QEta;}
  double** QMapping() const {return mapping;}

};


class Quadrature_2PQuad4 : public Quadrature{
protected:
public:
  virtual void Setup_Quadrature();
  virtual void Print_Quadrature_Info();
};


class Quadrature_3PQuad4 : public Quadrature{
public:
  virtual void Setup_Quadrature();
  virtual void Print_Quadrature_Info();
};



// functions


Quadrature :: ~Quadrature(){
  delete [] QW;
  delete [] QXi;
  delete [] QEta;
  delete [] mapping_data;
  delete [] mapping;
}


void Quadrature_2PQuad4 :: Setup_Quadrature(){
  Quadrature_points = 4;
  QW = new double [4];
  QXi = new double [4];
  QEta = new double [4];
  mapping_data = new double [16];
  mapping = new double* [4];
  for(int i = 0; i < 4; i++){
    mapping[i] = &mapping_data[i*4];
  }

  QW[0] = 1.0;
  QW[1] = 1.0;
  QW[2] = 1.0;
  QW[3] = 1.0;

  QXi[0] = -1.0/sqrt(3);
  QXi[1] =  1.0/sqrt(3);
  QXi[2] =  1.0/sqrt(3);
  QXi[3] = -1.0/sqrt(3);

  QEta[0] = -1.0/sqrt(3);
  QEta[1] = -1.0/sqrt(3);
  QEta[2] =  1.0/sqrt(3);
  QEta[3] =  1.0/sqrt(3);

  // transposed form because of lapack solver
  mapping[0][0] =  1.0;
  mapping[1][0] = -1.0;
  mapping[2][0] = -1.0;
  mapping[3][0] =  1.0;
  mapping[0][1] =  1.0;
  mapping[1][1] =  1.0;
  mapping[2][1] = -1.0;
  mapping[3][1] = -1.0;
  mapping[0][2] =  1.0;
  mapping[1][2] =  1.0;
  mapping[2][2] =  1.0;
  mapping[3][2] =  1.0;
  mapping[0][3] =  1.0;
  mapping[1][3] = -1.0;
  mapping[2][3] =  1.0;
  mapping[3][3] = -1.0;
}

void Quadrature_2PQuad4 :: Print_Quadrature_Info(){
  cout << "Quadrature = 2D Gaussian Quadrature : 4 Points" <<endl;
  cout << "Quadrature points : " << Qpoints() << endl;
  cout << "Weights:" << endl;
  cout << setw(10) << QW[0] << setw(10) << QW[1] << setw(10) << QW[2] << setw(10) << QW[3] << endl;
  cout << "Xi:" << endl;
  cout << setw(10) << QXi[0] << setw(10) << QXi[1] << setw(10) << QXi[2] << setw(10) << QXi[3] << endl;
  cout << "Eta:" << endl;
  cout << setw(10) << QEta[0] << setw(10) << QEta[1] << setw(10) << QEta[2] << setw(10) << QEta[3] << endl;
  cout << "Mapping matrix(Transposed):" << endl;
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      cout << setw(10) << mapping[i][j];
    }
    cout << endl;
  }
  cout << endl;
}



void Quadrature_3PQuad4 :: Setup_Quadrature(){
  Quadrature_points = 9;
  QW = new double [9];
  QXi = new double [9];
  QEta = new double [9];

  QW[0] = (5.0/9.0)*(5.0/9.0);
  QW[1] = (8.0/9.0)*(5.0/9.0);
  QW[2] = (5.0/9.0)*(5.0/9.0);
  QW[3] = (5.0/9.0)*(8.0/9.0);
  QW[4] = (8.0/9.0)*(8.0/9.0);
  QW[5] = (5.0/9.0)*(8.0/9.0);
  QW[6] = (5.0/9.0)*(5.0/9.0);
  QW[7] = (8.0/9.0)*(5.0/9.0);
  QW[8] = (5.0/9.0)*(5.0/9.0);

  QXi[0] = -sqrt(3.0/5.0);
  QXi[1] =  0.0;
  QXi[2] =  sqrt(3.0/5.0);
  QXi[3] = -sqrt(3.0/5.0);
  QXi[4] =  0.0;
  QXi[5] = -sqrt(3.0/5.0);
  QXi[6] =  sqrt(3.0/5.0);
  QXi[7] =  0.0;
  QXi[8] = -sqrt(3.0/5.0);

  QEta[0] = -sqrt(3.0/5.0);
  QEta[1] =  0.0;
  QEta[2] =  sqrt(3.0/5.0);
  QEta[3] = -sqrt(3.0/5.0);
  QEta[4] =  0.0;
  QEta[5] = -sqrt(3.0/5.0);
  QEta[6] =  sqrt(3.0/5.0);
  QEta[7] =  0.0;
  QEta[8] = -sqrt(3.0/5.0);
}

void Quadrature_3PQuad4 :: Print_Quadrature_Info(){
  cout << "Quadrature = 2D Gaussian Quadrature : 9 Points" <<endl;
  cout << "Quadrature points : " << Qpoints() << endl;
  cout << "Weights:" << endl;
  for(int i = 0; i < Qpoints(); i++)
    cout << setw(10) << QW[i];
  cout << endl;
  cout << "Xi:" << endl;
  for(int i = 0; i < Qpoints(); i++)
    cout << setw(10) << QXi[i];
  cout << endl;
  cout << "Eta:" << endl;
  for(int i = 0; i < Qpoints(); i++)
    cout << setw(10) << QEta[i];
  cout << endl;
  cout << endl;
}



#endif // QUADRATURE_HPP
