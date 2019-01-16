#include <iostream>
#include <fstream>
#include <chrono>

#include "eulerpoisson.h"

using namespace std;

int main(int argc, char** argv)
{
  double xmin = 0.0;
  double xmax = 1.0;

  int Nx = 100;
  double hx = (xmax-xmin)/(Nx+1.));

  double tfinal = 10.;
  double dt = 0.001;

  double gamma = 1.7;
  double g = 10;

  double CI_rho = 10.;
  double CI_u = 0.1;
  double CI_E = CI_rho/((gamma-1.)*CI_rho) + CI_u*CI_u/2.;

  int schema = 1;
  //1 -> Rusanov
  //2 -> Relaxation

  Schema_VF_1D *VF;

  if (Schema_VF_1D == 1){
    VF = new Rusanov();
    VF->Initialize(xmin,xmax,Nx,hx,dt,CI_rho,CI_u,CI_E,gamma,g);
    auto start = chrono::high_resolution_clock::now();
    VF->TimeScheme(tfinal);
    auto finish = chrono::high_resolution_clock::now();
  }

  if (Schema_VF_1D == 2){
    VF = new Relaxation();
    VF->Initialize(xmin,xmax,Nx,hx,dt,CI_rho,CI_u,CI_E,gamma,g);
    auto start = chrono::high_resolution_clock::now();
    VF->TimeScheme(tfinal);
    auto finish = chrono::high_resolution_clock::now();
  }

  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();

  // Affichage du chrono
  cout << "CPU time = "<< t << " seconds" << endl;

  delete VF;

  return 0;
}
