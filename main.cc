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
  double hx = (xmax-xmin)/(Nx+1.);

  double tfinal = 50.;
  double dt = 0.0001;

  double gamma = 1.4;
  double g = 10.;

  double a = 10.;

  double CI_rho = 10.;
  double CI_u = 0.01;
  double CI_E = 1/(gamma-1.) + CI_u*CI_u/2.;

  int schema = 2;
  //1 -> Rusanov
  //2 -> Relaxation

  Schema_VF_1D *VF;

  if (schema == 1){
    VF = new Rusanov();
    VF->Initialize(xmin,xmax,Nx,hx,dt,CI_rho,CI_u,CI_E,gamma,g,a);
    VF->TimeScheme(tfinal);
  }

  if (schema == 2){
    VF = new Relaxation();
    VF->Initialize(xmin,xmax,Nx,hx,dt,CI_rho,CI_u,CI_E,gamma,g,a);
    VF->TimeScheme(tfinal);
  }

  delete VF;

  return 0;
}
