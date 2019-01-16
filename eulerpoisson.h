#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>
#include <stdio.h>
#include <math.h>

#include "fonctions1D.h"

class Schema_VF_1D
{
 protected: // Les attributs de la classe

  double _xmin, _xmax, _hx, _dt;
  double _gamma; double _g;
  int _Nx;

  std::vector<std::vector<double> > _Wsol; // (rho, rho*u, rho*E)
  std::vector<std::vector<double> > _Wsol_moins; // vecteur solution à l'itération précédente
  std::vector<double> _Gravity; // (phi)
  std::vector<std::vector<double> > _LapMat1D; // matrice creuse du laplacien

 public: // Méthodes et opérateurs de la classe
  Schema_VF_1D();
  // Constructeur : Initialiser xmin, xmax, Nx, hx, dt, _Wsol_0
  virtual ~Schema_VF_1D();

  void Initialize(double xmin, double xmax, int Nx, double hx, double dt, double CI_rho, double CI_u, double CI_E, double gamma, double g);
  void SaveSol(const std::string& name_file, const int iter);
  void Poisson();
  virtual void TimeScheme(double tfinal) = 0;
  virtual void Euler() = 0;
  virtual void Source() = 0;
  // virtual void Flux() = 0;
};

class Rusanov : public Schema_VF_1D
{
 public:
   void TimeScheme(double tfinal);
   void Euler();
   void Source();
};

// class Relaxation : public Schema_VF_1D
// {
//  protected:
//   std::vector<std::vector<double> > _Wdelta; // (rho, u, epsilon, pi, psi)
//   std::vector<std::vector<double> > _Fdg; // flux gauche pour le problème de relaxation
//   std::vector<std::vector<double> > _Fdd; // flux droit pour le problème de relaxation
//
//  public:
//    void Initialize();
//    void TimeScheme(double tfinal);
//    void Flux();
// };
