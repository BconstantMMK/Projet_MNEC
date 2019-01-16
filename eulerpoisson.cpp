#include "eulerpoisson.h"

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
   double dt, double CI_rho, double CI_u, double CI_E,
      double CL_u_g, double CL_u_d, double CL_rho_g, double CL_rho_d,
          double CL_E_g, double CL_E_d, double CL_phi_g, double CL_phi_d) {

    _xmin = xmin; _xmax = xmax; _Nx = Nx; _hx = hx;
    _dt = dt;
    _CL_u_g = CL_u_g;
    _CL_u_d = CL_u_d;
    _CL_rho_g = CL_rho_g;
    _CL_rho_d = CL_rho_d;
    _CL_E_g = CL_E_g;
    _CL_E_d = CL_E_d;
    _CL_phi_g = CL_phi_g;
    _CL_phi_d = CL_phi_d;

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

    _LapMat1D.resize(3);
    double alpha = -2./(_hx * _hx);
    double beta = 1/(_hx * _hx);

    for (int i = 0; i < 3; i++)
    {
      _LapMat1D[i].resize(Nx);
    }
    for (int j = 0; j < _Nx; ++j){
      _LapMat1D[0][j] = beta;
      _LapMat1D[1][j] = alpha;
      _LapMat1D[2][j] = beta;
    }
    _LapMat1D[0][0] = 0;
    _LapMat1D[0][_Nx-1] = 0;
    _LapMat1D[2][0] = 0;
    _LapMat1D[2][_Nx-1] = 0;

}

void SaveSol(const std::string& name_file) {}

void Poisson() {}

//####################################################################
//Rusanov
//####################################################################

void Rusanov::Initialize(DataFile data_file)
{
  Schema_VF_1D::Initialize(double xmin, double xmax, int Nx, double hx,
     double dt, double CI_rho, double CI_u, double CI_E,
        double CL_u_g, double CL_u_d, double CL_rho_g, double CL_rho_d,
            double CL_E_g, double CL_E_d, double CL_phi_g, double CL_phi_d);

  _Fg[].resize(3);
  _Fd[].resize(3);

  for (int i = 0; i < 3; ++i)
   {
     _Fg[i].resize(Nx);
     _Fd[i].resize(Nx);
   }

}

//####################################################################
//Relaxation
//####################################################################

void Relaxation::Initialize(DataFile data_file)
{
  Schema_VF_1D::Initialize(double xmin, double xmax, int Nx, double hx,
     double dt, double CI_rho, double CI_u, double CI_E);

  _Fdg[].resize(5);
  _Fdd[].resize(5);
  _Wdelta[].resize(5);

  for (int i = 0; i < 5; ++i)
   {
     _Fdg[i].resize(Nx);
     _Fdd[i].resize(Nx);
     _Wdelta[i].resize(Nx);
   }

}
