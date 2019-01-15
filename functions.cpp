#include "functions.h"

using namespace std;

//Constructeur :
Schema_VF_1D::Schema_VF_1D()
{
}
//Destructeur :
Schema_VF_1D::~Schema_VF_1D()
{
}

void Initialize(double xmin, double xmax, int Nx, double hx,
   double dt, double CI_rho, double CI_u, double CI_E) {

    _xmin = xmin; _xmax = xmax; _Nx = Nx; _hx = hx;
    _dt = dt; nb_iterations


    _Gravity.resize(Nx);
    _Wsol[].resize(3);
    _Wsol_moins[].resize(3);

    for (int i = 0; i < 3; ++i)
    {
      _Wsol[i].resize(Nx);
      _Wsol_moins[i].resize(Nx);
    }

    for (int j = 0; j < _Nx; ++j){
      _Wsol[0][j] = CI_rho;
      _Wsol[1][j] = CI_rho*CI_u;
      _Wsol[2][j] = CI_rho*CI_E;
    }

}

void SaveSol(const std::string& name_file) {}

void Poisson() {}

//####################################################################
//Rusanov
//####################################################################

void Rusanov::Initialize(DataFile data_file)
{
  Schema_VF_1D::Initialize(double xmin, double xmax, int Nx, double hx,
     double dt, double CI_rho, double CI_u, double CI_E);

  _Fl[].resize(3);
  _Fr[].resize(3);

  for (int i = 0; i < 3; ++i)
   {
     _Fl[i].resize(Nx);
     _Fr[i].resize(Nx);
   }

}

//####################################################################
//Relaxation
//####################################################################

void Relaxation::Initialize(DataFile data_file)
{
  Schema_VF_1D::Initialize(double xmin, double xmax, int Nx, double hx,
     double dt, double CI_rho, double CI_u, double CI_E);

  _Fdl[].resize(5);
  _Fdr[].resize(5);
  _Wd[].resize(5);

  for (int i = 0; i < 5; ++i)
   {
     _Fdl[i].resize(Nx);
     _Fdr[i].resize(Nx);
     _Wd[i].resize(Nx);
   }

}
