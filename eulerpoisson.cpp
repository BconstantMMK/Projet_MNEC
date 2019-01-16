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
   double dt, double CI_rho, double CI_u, double gamma, double g);
{
    _xmin = xmin; _xmax = xmax; _Nx = Nx; _hx = hx;
    _dt = dt;
    _gamma = gamma;
    _g = g;

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

void SaveSol(const std::string& name_file, const int iter)
{
  ofstream flux_sol;
  flux_sol.open(name_file + ".txt", ios::out);
  for (int j = 0; j < _Nx; ++j) {
    flux_sol << iter << " " << (j*_hx/2) << " " << _Wsol[0][j] << " " << _Wsol[1][j]/_Wsol[0][j] << " " _Wsol[2][j]/_Wsol[0][j] << endl;
  }
  flux_sol.close();
}

void Poisson()
{
  double CL_phi_g = 0;
  double CL_phi_d = g;
  double pi = acos(-1);
  double G = 6.67*pow(10,-11);
  int kmax = _Nx + 100;
  std::vector<double> b;
  b.resize(_Nx);

  b[0] -= CL_phi_g/(_hx*_hx);
  b[_Nx-1] -= Cm_phi_d/(_hx*_hx);

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
     double dt, double CI_rho, double CI_u, double CI_E, double gamma, double g);

  Schema_VF_1D::Poisson();
}

void Rusanov::Euler()
{
  double bi=0.;
  _Wsol_moins = _Wsol;
  Flux();
  for (int j = 1; j < _Nx-1; ++j){
    bi = max({abs(_Wsol_moins[1][j]/_Wsol_moins[0][j]-1),abs(_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1]-1),abs(_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1]-1),
    abs(_Wsol_moins[1][j]/_Wsol_moins[0][j]+1),abs(_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1]+1),abs(_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1]+1),});
    _Wsol[0][j] = _Wsol_moins[0][j] - (_dt/_hx)*(0.5*(_Wsol_moins[1][j+1]-_Wsol_moins[1][j-1])
                                                  - bi*0.5*(_Wsol_moins[0][j+1] -2*_Wsol_moins[0][j] + _Wsol_moins[0][j-1]));
    _Wsol[1][j] = _Wsol_moins[1][j] - (_dt/_hx)*(0.5*((_Wsol_moins[1][j+1]*_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1] + _Wsol_moins[0][j+1])-(_Wsol_moins[1][j-1]*_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1] + _Wsol_moins[0][j-1]))
                                                  - bi*0.5*(_Wsol_moins[1][j+1] -2*_Wsol_moins[1][j] + _Wsol_moins[1][j-1]));
    _Wsol[2][j] = _Wsol_moins[2][j] - (_dt/_hx)*(0.5*((_Wsol_moins[2][j+1]+_Wsol_moins[0][j+1])*_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1])-(_Wsol_moins[2][j-1]+_Wsol_moins[0][j-1])*_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1]))
                                                  - bi*0.5*(_Wsol_moins[2][j+1] -2*_Wsol_moins[2][j] + _Wsol_moins[2][j-1]));
  }

  //CL gauche : u=0 / -dx(rho) = -rho*g / -dx(E) = gE - g/(gamma-1)
  double CL_rho_g = _hx*_Wsol_moins[0][0]*_g + _Wsol_moins[0][0];
  double CL_rhoE_g = _hx*_Wsol_moins[2][0]*_g -_g/(_gamma-1) + _Wsol_moins[2][0];

  bi = max({abs(_Wsol_moins[1][0]/_Wsol_moins[0][0]-1),abs(0-1),abs(_Wsol_moins[1][1]/_Wsol_moins[0][1]-1),
  abs(_Wsol_moins[1][0]/_Wsol_moins[0][0]+1),abs(0+1),abs(_Wsol_moins[1][1]/_Wsol_moins[0][1]+1),});

  _Wsol[0][0] = _Wsol_moins[0][0] - (_dt/_hx)*(0.5*(_Wsol_moins[1][1]-0)
                                                - bi*0.5*(_Wsol_moins[0][1] -2*_Wsol_moins[0][0] + CL_rho_g));
  _Wsol[1][0] = _Wsol_moins[1][0] - (_dt/_hx)*(0.5*((_Wsol_moins[1][1]*_Wsol_moins[1][1]/_Wsol_moins[0][1] + _Wsol_moins[0][1])-(0 + CL_rho_g))
                                                - bi*0.5*(_Wsol_moins[1][1] -2*_Wsol_moins[1][0] + 0));
  _Wsol[2][0] = _Wsol_moins[2][0] - (_dt/_hx)*(0.5*(((_Wsol_moins[2][1]+_Wsol_moins[0][1])*_Wsol_moins[1][1]/_Wsol_moins[0][1])-(0))
                                                - bi*0.5*(_Wsol_moins[2][1] -2*_Wsol_moins[2][0] + CL_rhoE_g));


  //CL droite : u=0 / dx(rho) = -rho*g / dx(E) = gE - g/(gamma-1)
  double CL_rho_d = _hx*_Wsol_moins[0][_Nx-1]*_g + _Wsol_moins[0][_Nx-1];
  double CL_rhoE_d = _hx*_Wsol_moins[2][_Nx-1]*_g -_g/(_gamma-1) + _Wsol_moins[2][_Nx-1];

  bi = max({abs(_Wsol_moins[1][_Nx-1]/_Wsol_moins[0][_Nx-1]-1),abs(0-1),abs(_Wsol_moins[1][_Nx-2]/_Wsol_moins[0][_Nx-2]-1),
  abs(_Wsol_moins[1][_Nx-1]/_Wsol_moins[0][_Nx-1]+1),abs(0+1),abs(_Wsol_moins[1][_Nx-2]/_Wsol_moins[0][_Nx-2]+1),});

  _Wsol[0][_Nx-1] = _Wsol_moins[0][_Nx-2] - (_dt/_hx)*(0.5*(0-_Wsol_moins[1][_Nx-2])
                                                - bi*0.5*(_Wsol_moins[0][_Nx-2] -2*_Wsol_moins[0][_Nx-1] + CL_rho_d));
  _Wsol[1][_Nx-1] = _Wsol_moins[1][_Nx-2] - (_dt/_hx)*(0.5*((0 + CL_rho_d)-(_Wsol_moins[1][_Nx-2]*_Wsol_moins[1][_Nx-2]/_Wsol_moins[0][_Nx-2] + _Wsol_moins[0][_Nx-2]))
                                                - bi*0.5*(_Wsol_moins[1][_Nx-2] -2*_Wsol_moins[1][_Nx-1] + 0));
  _Wsol[2][_Nx-1] = _Wsol_moins[2][_Nx-2] - (_dt/_hx)*(0.5*(0-((_Wsol_moins[2][1]+_Wsol_moins[0][1])*_Wsol_moins[1][1]/_Wsol_moins[0][1]))
                                                - bi*0.5*(_Wsol_moins[2][_Nx-2] -2*_Wsol_moins[2][_Nx-1] + CL_rhoE_d));
}

void Rusanov::Source()
{
  double CL_phi_g = 0;
  _Wsol_moins = _Wsol;
  for (int j = 1; j < _Nx; ++j){
    _Wsol[1][j] = _Wsol_moins[1][j] - _dt*(_Wsol_moins[0][j]*((_Gravity[j]-_Gravity[j-1])/_hx));
    _Wsol[2][j] = _Wsol_moins[2][j] - _dt*(_Wsol_moins[1][j]*((_Gravity[j]-_Gravity[j-1])/_hx));
  }
  _Wsol[1][_Nx-1] = _Wsol_moins[1][_Nx-1] - _dt*(_Wsol_moins[0][_Nx-1]*((_Gravity[j]-CL_phi_g)/_hx));
  _Wsol[2][_Nx-1] = _Wsol_moins[2][_Nx-1] - _dt*(_Wsol_moins[1][_Nx-1]*((_Gravity[j]-CL_phi_g)/_hx));
}

void Rusanov::TimeScheme(double tfinal) {

  int nbiter = int(ceil(tfinal / _dt));
  for (int iter=0; iter<nbiter; ++n){
    Euler();
    Source();
    Schema_VF_1D::Poisson();
    Schema_VF_1D::SaveSol("EquilibriumFlow1D", iter);
    //Barre de Chargement ------------------------------------------
    int i_barre;
    const int p = floor((((double)iter) / ((double)nbiter)) * 100);
    printf("[");
    for (i_barre = 0; i_barre <= p; i_barre += 2)
    printf("*");
    for (; i_barre <= 100; i_barre += 2)
    printf("-");
    printf("] %3d %%", p);
    for (i_barre = 0; i_barre < 59; ++i_barre)
    printf("%c", 8);
    fflush(stdout);
    //--------------------------------------------------------------
  }
}


//####################################################################
//Relaxation
//####################################################################

void Relaxation::Initialize(DataFile data_file)
{
  Schema_VF_1D::Initialize(double xmin, double xmax, int Nx, double hx,
     double dt, double CI_rho, double CI_u, double CI_E, double gamma);

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
