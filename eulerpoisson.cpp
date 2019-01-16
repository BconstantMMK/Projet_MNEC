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
          double CL_E_g, double CL_E_d, double CL_phi_g, double CL_phi_d)
{

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

void Poisson()
{
  double pi = acos(-1);
  double G = 6.67*pow(10,-11);
  int kmax = _Nx + 100;
  std::vector<double> b;
  b.resize(_Nx);

  for (int i = 0; i < _Nx; ++i){
    b[i] = 4*pi*G*_Wsol[0][i];
  }

  _Gravity = CG(_LapMat1D, b, _Gravity, 0.000001, kmax, _Nx);
}

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

  Schema_VF_1D::Poisson();

}

void Rusanov::Flux() {}

void Rusanov::Euler()
{
  _Wsol_moins = _Wsol;
  Flux();
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < _Nx; ++j){
      _Wsol[i][j] = _Wsol_moins[i][j] - (_dt/_hx)*(_Fg[i][j] - _Fd[i][j]);
    }
  }
}

void Rusanov::Source()
{
  _Wsol_moins = _Wsol;
  for (int j = 0; j < _Nx; ++j){
    _Wsol[1][j] = _Wsol_moins[1][j] - _dt*(_Wsol_moins[0][j]*((_Gravity[j]-_Gravity[j-1])/_hx));
    _Wsol[2][j] = _Wsol_moins[2][j] - _dt*(_Wsol_moins[1][j]*((_Gravity[j]-_Gravity[j-1])/_hx));
  }

}

void Rusanov::TimeScheme(double tfinal) {

  int nbiter = int(ceil(tfinal / _dt));
  for (int n=0; n<nbiter; ++n){
    Euler();
    Source();
    Schema_VF_1D::Poisson();
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

void Rusanov::Flux() {}

void Rusanov::TimeScheme(double tfinal) {}
