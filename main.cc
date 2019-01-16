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

  int nb_iterations = int(ceil(tfinal / deltaT));

  double CI_rho = 10.;
  double CI_u = 1.;
  double CI_E = 100.;

  double CL_u_g = 0.;
  double CL_u_d = 0.;
  double CL_rho_g = 0.;
  double CL_rho_d = 0.;
  double CL_E_g = 0.;
  double CL_E_d = 0.;
  double CL_phi_g = 0.;
  double CL_phi_d = 0.;

  int schema = 1;
  //1 -> Rusanov
  //2 -> Relaxation

  Schema_VF_1D *VF;

  if (Schema_VF_1D == 1){
    VF = new Rusanov();
    VF->Initialize(xmin,xmax,Nx,hx,dt,CI_rho,CI_u,CI_E,
    CL_u_g,CL_u_d,CL_rho_g,CL_rho_d,CL_E_g,CL_E_d,CL_phi_g,CL_phi_d);
    auto start = chrono::high_resolution_clock::now();
    VF->TimeScheme(tfinal);
    auto finish = chrono::high_resolution_clock::now();
  }

  if (Schema_VF_1D == 2){
    VF = new Relaxation();
    VF->Initialize(xmin,xmax,Nx,hx,dt,CI_rho,CI_u,CI_E,
    CL_u_g,CL_u_d,CL_rho_g,CL_rho_d,CL_E_g,CL_E_d,CL_phi_g,CL_phi_d);
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
