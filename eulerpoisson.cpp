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
  _Wsol_moins = _Wsol;

  for (int j = 1; j < _Nx+1; ++j){
    _Wsol[1][j] = _Wsol_moins[1][j] - _dt*(_Wsol_moins[0][j]*_g);
    _Wsol[2][j] = _Wsol_moins[2][j] - _dt*(_Wsol_moins[1][j]*_g);
  }
}

void Rusanov::TimeScheme(double tfinal)
{
  int nbiter = int(ceil(tfinal / _dt));
  double normeL2rho;
  double normeL2u;

  ofstream flux_norme;
  flux_norme.open("EF1D_Rusanov_norme.txt", ios::out);

  for (int iter=0; iter<nbiter; ++iter){
    Euler();
    if (_SFSG == false)
      break;
    Source();

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
  Schema_VF_1D::SaveSol("EF1D_Rusanov");
  flux_norme.close();
}

//####################################################################
//Relaxation
//####################################################################

void Relaxation::Initialize(double xmin, double xmax, int Nx, double hx, double dt, double CI_rho, double CI_u, double CI_E, double gamma, double g, double a)
{
  Schema_VF_1D::Initialize(xmin, xmax, Nx, hx, dt, CI_rho, CI_u, CI_E, gamma, g, a);
  _Fdg.resize(3);
  _Fdd.resize(3);
  _Wdelta.resize(5);

  _a = a;

  for (int i = 0; i < 5; ++i)
   {
     _Wdelta[i].resize(_Nx+2);
   }
  for (int i = 0; i < 3; ++i)
  {
    _Fdg[i].resize(_Nx+2);
    _Fdd[i].resize(_Nx+2);
  }
}

void Relaxation::UpdateFluxCase1(string sens, int j, double sigma, double Sl, double Sr)
{
  vector<double> Vl(4,0.),Vll(4,0.),Vr(4,0.);
  if (sens == "g"){
    _Fdg[0][j] = _hx/2.*_dt*_Wsol_moins[0][j] + _hx/2.*_dt*_Wsol_moins[0][j];
    _Fdg[1][j] = _hx/2.*_dt*_Wsol_moins[1][j] + _hx/2.*_dt*_Wsol_moins[1][j];
    _Fdg[2][j] = _hx/2.*_dt*_Wsol_moins[2][j] + _hx/2.*_dt*_Wsol_moins[2][j];
  }

  if (sens == "d"){
    Vl[0] = _Wdelta[0][j-1]*sqrt(1. - 2.*(_Wdelta[4][j] - _Wdelta[4][j-1])/(pow(_Wdelta[1][j-1],2) - pow(_a,2)*pow(_Wdelta[0][j-1],2)));
    Vll[0] = 1./(2.*_a)*(Sr - Sl*Vl[0]/_Wdelta[0][j-1] - (_Wdelta[3][j] + pow(_a,2)*_Wdelta[0][j] - _Wdelta[3][j]+pow(_a,2)*_Wdelta[0][j]));
    Vr[0] = 1./(2.*_a)*(Sr - Sl*Vl[0]/_Wdelta[0][j-1] + (_Wdelta[3][j] + pow(_a,2)*_Wdelta[0][j] - _Wdelta[3][j]+pow(_a,2)*_Wdelta[0][j]));

    Vl[1] = _Wdelta[1][j-1]*Vl[0]/_Wdelta[0][j-1];
    Vll[1] = sigma;
    Vr[1] = sigma;

    Vl[3] = _Wdelta[3][j-1] + pow(_a,2)*_Wdelta[0][j-1] - pow(_a,2)*Vl[0];
    Vll[3] = _Wdelta[3][j-1] + pow(_a,2)*_Wdelta[0][j-1] - pow(_a,2)*Vll[0];
    Vr[3] = _Wdelta[3][j] + pow(_a,2)*_Wdelta[0][j] - pow(_a,2)*Vr[0];

    Vl[2] = _Wdelta[2][j-1] - pow(_Wdelta[3][j-1],2)/(2*pow(_a,2)) + pow(Vl[3],2)/(2*pow(_a,2));
    Vll[2] = _Wdelta[2][j-1] - pow(_Wdelta[3][j-1],2)/(2*pow(_a,2)) + pow(Vll[3],2)/(2*pow(_a,2));
    Vr[2] = _Wdelta[2][j] - pow(_Wdelta[3][j],2)/(2*pow(_a,2)) + pow(Vr[3],2)/(2*pow(_a,2));

    _Fdg[0][j] = _hx/2.*_dt*_Wsol_moins[0][j] + Sl*(1./Vl[0]) + (sigma-Sl)*(1./Vll[0]) + (Sr-sigma)*(1./Vr[0]) + (_hx/2.*_dt-Sr)*_Wsol_moins[0][j];
    _Fdg[1][j] = _hx/2.*_dt*_Wsol_moins[1][j] + Sl*(Vl[1]/Vl[0]) + (sigma-Sl)*(Vll[1]/Vll[0]) + (Sr-sigma)*(Vr[1]/Vr[0]) + (_hx/2.*_dt-Sr)*_Wsol_moins[1][j];
    _Fdg[2][j] = _hx/2.*_dt*_Wsol_moins[2][j] + Sl*(Vl[2]+pow(Vl[1],2)/2.)/Vl[0] + (sigma-Sl)*(Vll[2]+pow(Vll[1],2)/2.)/Vll[0] + (Sr-sigma)*(Vr[2]+pow(Vr[1],2)/2.)/Vr[0] + (_hx/2.*_dt-Sr)*_Wsol_moins[2][j];

  }
}

void Relaxation::UpdateFluxCase2(string sens, int j, double sigma, double Sl, double Sr)
{
  vector<double> Vl(5,0.),Vll(5,0.),Vr(5,0.);
}

void Relaxation::UpdateFluxCase3(string sens, int j, double sigma, double Sl, double Sr)
{
  vector<double> Vl(5,0.),Vll(5,0.),Vr(5,0.);
}

void Relaxation::UpdateFluxCase4(string sens, int j, double sigma, double Sl, double Sr)
{
  vector<double> Vl(5,0.),Vll(5,0.),Vr(5,0.);
}

void Relaxation::Flux()
{
  _Wsol_moins = _Wsol;
  double sigma, Sl, Sr;
  double bi=0.;
  vector<double> b;


  for (int j=0; j < _Nx+2; ++j)
  {
    b = {abs(_Wsol_moins[1][j]/_Wsol_moins[0][j]-1),abs(_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1]-1),abs(_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1]-1),
    abs(_Wsol_moins[1][j]/_Wsol_moins[0][j]+1),abs(_Wsol_moins[1][j-1]/_Wsol_moins[0][j-1]+1),abs(_Wsol_moins[1][j+1]/_Wsol_moins[0][j+1]+1)};

    bi = *max_element(begin(b),end(b));

    if (bi*_dt/_hx > 0.5)
    {
      cout << "ATTENTION, LA CONDTION CFL N'EST PAS VERIFIÉE !" << endl;
      cout << "Relancez avec un pas de temps plus petit." << endl;
      _SFSG = false;
    }

    _Wdelta[0][j] = 1/_Wsol_moins[0][j];
    _Wdelta[1][j] = _Wsol_moins[1][j]/_Wsol_moins[0][j];
    _Wdelta[2][j] = _Wsol_moins[2][j]/_Wsol_moins[0][j] - pow(_Wsol_moins[1][j],2)/_Wsol_moins[0][j];
    _Wdelta[3][j] = _Wsol_moins[0][j];
    _Wdelta[4][j] = (_hx/2+j*_hx)*_g;
  }

  Schema_VF_1D::UpdateCL();

  //Gauche
  for (int j=1; j < _Nx+1; ++j)
  {
    sigma = (_Wdelta[1][j]+_Wdelta[1][j+1])*0.5 - (_Wdelta[3][j+1] - _Wdelta[3][j])/(2*_a);
    Sl = _Wdelta[1][j]-_a*_Wdelta[0][j];
    Sr = _Wdelta[1][j+1]+_a*_Wdelta[0][j+1];
    if(Sl < 0)
      UpdateFluxCase1("g",j,sigma,Sl,Sr);
    else if(Sr < 0)
      UpdateFluxCase4("g",j,sigma,Sl,Sr);
    else if(sigma > 0)
      UpdateFluxCase2("g",j,sigma,Sl,Sr);
    else
      UpdateFluxCase3("g",j,sigma,Sl,Sr);

    sigma = (_Wdelta[1][j-1]+_Wdelta[1][j])*0.5 - (_Wdelta[3][j] - _Wdelta[3][j-1])/(2*_a);
    Sl = _Wdelta[1][j-1]-_a*_Wdelta[0][j-1];
    Sr = _Wdelta[1][j]+_a*_Wdelta[0][j];
    if(Sl < 0)
      UpdateFluxCase1("d",j,sigma,Sl,Sr);
    else if(Sr < 0)
      UpdateFluxCase4("d",j,sigma,Sl,Sr);
    else if(sigma > 0)
      UpdateFluxCase2("d",j,sigma,Sl,Sr);
    else
      UpdateFluxCase3("d",j,sigma,Sl,Sr);
  }
}

void Relaxation::TimeScheme(double tfinal)
{
  int nbiter = int(ceil(tfinal / _dt));
  double normeL2rho;
  double normeL2u;

  ofstream flux_norme;
  flux_norme.open("EF1D_Relaxation_norme.txt", ios::out);

  for (int iter=0; iter<nbiter; ++iter){
    Flux();
    if (_SFSG == false)
      break;

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
  Schema_VF_1D::SaveSol("EF1D_Relaxation");
  flux_norme.close();
}
