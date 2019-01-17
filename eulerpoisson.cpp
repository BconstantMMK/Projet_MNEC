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

void Schema_VF_1D::Initialize(double xmin, double xmax, int Nx, double hx, double dt, double CI_rho, double CI_u, double CI_E, double gamma, double g, double a)
{
    _xmin = xmin; _xmax = xmax; _Nx = Nx; _hx = hx;
    _dt = dt;
    _gamma = gamma;
    _g = g;
    _SFSG = true;
    _CI_rho = CI_rho;
    _a = a;

    _Gravity.resize(Nx);
    _Wsol.resize(3);
    _Wsol_moins.resize(3);

    for (int i = 0; i < 3; ++i)
    {
      _Wsol[i].resize(Nx+2);
      _Wsol_moins[i].resize(Nx+2);
    }

    for (int j = 0; j < _Nx; ++j){
      _Wsol[0][j+1] = CI_rho;
      _Wsol[1][j+1] = CI_rho*CI_u;
      _Wsol[2][j+1] = CI_rho*CI_E;
      _Gravity[j] = j*_hx*g;
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
    _LapMat1D[0][0] = 0.;
    _LapMat1D[0][_Nx-1] = 0.;
    _LapMat1D[2][0] = 0.;
    _LapMat1D[2][_Nx-1] = 0.;
}

void Schema_VF_1D::UpdateCL()
{
  //CL gauche : u=0 / dx(rho) = -rho_eq*g / dx(E) = gE - g/(gamma-1)
  _Wsol_moins[0][0] = _CI_rho;
  // _Wsol_moins[0][0] = _Wsol_moins[0][1] - _hx*_g*_CI_rho;
  _Wsol_moins[1][0] = 0.;
  _Wsol_moins[2][0] = 1./(1+(_g*_hx))*(_Wsol_moins[2][1]/_Wsol_moins[0][1]+_g*_hx/(_gamma-1));
  _Wsol_moins[2][0] *= _Wsol_moins[0][0];

  //CL droite : u=0 / dx(rho) = -rho_eq*g / dx(E) = gE - g/(gamma-1)
  _Wsol_moins[0][_Nx+1] = _CI_rho*exp(-_g);
  // _Wsol_moins[0][_Nx+1] = _Wsol_moins[0][_Nx] - _hx*_g*_CI_rho*exp(-_g);
  _Wsol_moins[1][_Nx+1] = 0.;
  _Wsol_moins[2][_Nx+1] = 1./(1-(_g*_hx))*(_Wsol_moins[2][_Nx]/_Wsol_moins[0][_Nx]-_g*_hx/(_gamma-1));
  _Wsol_moins[2][_Nx+1] *= _Wsol_moins[0][_Nx+1];
}

void Schema_VF_1D::SaveSol(const std::string& name_file)
{
  ofstream flux_sol;
  flux_sol.open(name_file + ".txt", ios::out);
  for (int j = 0; j < _Nx; ++j) {
    flux_sol << (_hx/2+j*_hx) << " " << _Wsol[0][j+1] << " " << 10*exp(-_g*(_hx/2+j*_hx)) << " "<< _Wsol[1][j+1]/_Wsol[0][j+1] << " " << _Wsol[2][j+1]/_Wsol[0][j+1] << endl;
  }
  flux_sol.close();
}

void Schema_VF_1D::Poisson()
{
  double CL_phi_g = 0;
  double CL_phi_d = _g;
  double pi = acos(-1);
  double G = 6.67*pow(10,-11);
  int kmax = _Nx + 100;
  vector<double> x0(_Nx,0.);
  vector<double> b;
  b.resize(_Nx);

  b[0] -= CL_phi_g/(_hx*_hx);
  b[_Nx-1] -= CL_phi_d/(_hx*_hx);

  for (int i = 0; i < _Nx; ++i){
    b[i] = 4.*pi*G*_Wsol[0][i+1];
  }

  _Gravity = CG(_LapMat1D, b, x0, 0.000001, kmax, _Nx);
}

//####################################################################
//Rusanov
//####################################################################

void Rusanov::Euler()
{
  double bi=0.;
  vector<double> b;
  _Wsol_moins = _Wsol;

  Schema_VF_1D::UpdateCL();

  for (int j = 1; j < _Nx+1; ++j){
    b = {abs(_Wsol_moins[1][j]/_Wsol_moins[0][j]-1),abs(_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1]-1),abs(_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1]-1),
    abs(_Wsol_moins[1][j]/_Wsol_moins[0][j]+1),abs(_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1]+1),abs(_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1]+1)};

    bi = *max_element(begin(b),end(b));

    if (bi*_dt/_hx > 0.5)
    {
      cout << "ATTENTION, LA CONDTION CFL N'EST PAS VERIFIÉE !" << endl;
      cout << "Relancez avec un pas de temps plus petit." << endl;
      _SFSG = false;
      break;
    }

    _Wsol[0][j] = _Wsol_moins[0][j] - (_dt/_hx)*(0.5*(_Wsol_moins[1][j+1]-_Wsol_moins[1][j-1])
                                                  - bi*0.5*(_Wsol_moins[0][j+1] -2*_Wsol_moins[0][j] + _Wsol_moins[0][j-1]));
    _Wsol[1][j] = _Wsol_moins[1][j] - (_dt/_hx)*(0.5*((_Wsol_moins[1][j+1]*_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1] + _Wsol_moins[0][j+1])-(_Wsol_moins[1][j-1]*_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1] + _Wsol_moins[0][j-1]))
                                                  - bi*0.5*(_Wsol_moins[1][j+1] -2*_Wsol_moins[1][j] + _Wsol_moins[1][j-1]));
    _Wsol[2][j] = _Wsol_moins[2][j] - (_dt/_hx)*(0.5*(((_Wsol_moins[2][j+1]+_Wsol_moins[0][j+1])*_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1])-((_Wsol_moins[2][j-1]+_Wsol_moins[0][j-1])*_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1]))
                                                  - bi*0.5*(_Wsol_moins[2][j+1] -2*_Wsol_moins[2][j] + _Wsol_moins[2][j-1]));
  }
}

void Rusanov::Source()
{
  double CL_phi_g = 0;
  _Wsol_moins = _Wsol;
  // for (int j = 1; j < _Nx; ++j){
  //   _Wsol[1][j+1] = _Wsol_moins[1][j+1] - _dt*(_Wsol_moins[0][j+1]*((_Gravity[j]-_Gravity[j-1])/_hx));
  //   _Wsol[2][j+1] = _Wsol_moins[2][j+1] - _dt*(_Wsol_moins[1][j+1]*((_Gravity[j]-_Gravity[j-1])/_hx));
  // }
  // _Wsol[1][1] = _Wsol_moins[1][1] - _dt*(_Wsol_moins[0][1]*((_Gravity[0]-CL_phi_g)/_hx));
  // _Wsol[2][1] = _Wsol_moins[2][1] - _dt*(_Wsol_moins[1][1]*((_Gravity[0]-CL_phi_g)/_hx));

  for (int j = 1; j < _Nx+1; ++j){
    _Wsol[1][j] = _Wsol_moins[1][j] - _dt*(_Wsol_moins[0][j]*_g);
    _Wsol[2][j] = _Wsol_moins[2][j] - _dt*(_Wsol_moins[1][j]*_g);
  }
}

void Rusanov::TimeScheme(double tfinal) {

  // Schema_VF_1D::Poisson();

  int nbiter = int(ceil(tfinal / _dt));
  double normeL2rho;
  double normeL2u;

  ofstream flux_norme;
  flux_norme.open("EquilibriumFlow1D_norme.txt", ios::out);

  for (int iter=0; iter<nbiter; ++iter){
    Euler();
    if (_SFSG == false)
      break;
    Source();
    // Schema_VF_1D::Poisson();

    //Calcul des normes --------------------------------------------
    if (iter % 1000 == 0){
      normeL2rho = 0;
      normeL2u = 0;
      for (int j = 0; j < _Nx; ++j) {
        normeL2rho += (_Wsol[0][j+1]-10*exp(-_g*(_hx/2+j*_hx)))*(_Wsol[0][j+1]-10*exp(-_g*(_hx/2+j*_hx)));
        normeL2u += _Wsol[1][j+1]/_Wsol[0][j+1]*_Wsol[1][j+1]/_Wsol[0][j+1];
      }
      normeL2rho = log(sqrt(_hx*normeL2rho));
      normeL2u = log(sqrt(_hx*normeL2u));
      flux_norme << _dt*iter << " " << normeL2rho << " " << normeL2u << endl;
    }
    //--------------------------------------------------------------

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
  Schema_VF_1D::SaveSol("EquilibriumFlow1D");
  flux_norme.close();
}

//####################################################################
//Relaxation
//####################################################################

void Relaxation::Initialize(double xmin, double xmax, int Nx, double hx, double dt, double CI_rho, double CI_u, double CI_E, double gamma, double g, double a)
{
  Schema_VF_1D::Initialize(xmin, xmax, Nx, hx, dt, CI_rho, CI_u, CI_E, gamma, g, a);
  _Fdg.resize(5);
  _Fdd.resize(5);
  _Wdelta.resize(5);

  _a = a;

  for (int i = 0; i < 5; ++i)
   {
     _Fdg[i].resize(_Nx+2);
     _Fdd[i].resize(_Nx+2);
     _Wdelta[i].resize(_Nx+2);
   }
}

void Relaxation::UpdateFluxCase1(string sens, int j, double sigma)
{

}

void Relaxation::UpdateFluxCase2(string sens, int j, double sigma)
{

}

void Relaxation::UpdateFluxCase3(string sens, int j, double sigma)
{

}

void Relaxation::UpdateFluxCase4(string sens, int j, double sigma)
{

}

void Relaxation::Flux()
{
  _Wsol_moins = _Wsol;
  double sigma_gauche;
  double sigma_droite;

  for (int j=0; j < _Nx+2; ++j)
  {
    _Wdelta[0][j] = 1/_Wsol_moins[0][j];
    _Wdelta[1][j] = _Wsol_moins[1][j]/Wsol_moins[0][j];
    _Wdelta[2][j] = Wsol_moins[2][j]/Wsol_moins[0][j] - pow(Wsol_moins[1][j],2)/Wsol_moins[0][j];
    _Wdelta[3][j] = _Wsol_moins[0][j];
    _Wdelta[4][j] = (_hx/2+j*_hx)*_g;
  }

  Schema_VF_1D::UpdateCL();

  //Gauche
  for (int j=1; j < _Nx+1; ++j)
  {
    sigma_gauche = (_Wdelta[1][j]+_Wdelta[1][j+1])*0.5 - (_Wdelta[3][j+1] - _Wdelta[3][j])/(2*_a);
    if(_Wdelta[1][j]-_a*_Wdelta[0][j] < 0)
      UpdateFluxCase1("g",j,sigma);
    else if(_Wdelta[1][j+1]+_a*_Wdelta[0][j+1] < 0)
      UpdateFluxCase4("g",j,sigma);
    else if(sigma_gauche > 0)
      UpdateFluxCase2("g",j,sigma);
    else
      UpdateFluxCase3("g",j,sigma);

    sigma_droite = (_Wdelta[1][j-1]+_Wdelta[1][j])*0.5 - (_Wdelta[3][j] - _Wdelta[3][j-1])/(2*_a);
    if(_Wdelta[1][j-1]-_a*_Wdelta[0][j-1] < 0)
      UpdateFluxCase1("d",j,sigma);
    else if(_Wdelta[1][j]+_a*_Wdelta[0][j] < 0)
      UpdateFluxCase4("d",j,sigma);
    else if(sigma_droite > 0)
      UpdateFluxCase2("d",j,sigma);
    else
      UpdateFluxCase3("d",j,sigma);
  }
}

void Relaxation::TimeScheme(double tfinal)
{

}
