#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>
#include <stdio.h>
#include <math.h>

class Schema_VF_1D
{
 protected: // Les attributs de la classe

  double _x_min, _x_max, _h_x, _dt;
  int _Nx;

  std::vector<std::vector<double> > _Wsol; // (rho, rho*u, rho*E)
  std::vector<std::vector<double> > _Wsol_moins;
  std::vector<std::vector<double> > _Gravity; // (phi)

 public: // Méthodes et opérateurs de la classe
  Schema_VF_1D();
  // Constructeur : Initialiser xmin, xmax, Nx, hx, dt, _Wsol_0
  virtual ~Schema_VF_1D();

  void Initialize(xmin,xmax,Nx,hx,dt,CI_rho,CI_u,CI_E);
  void SaveSol(const std::string& name_file);
  void Poisson();
  virtual void TimeScheme(tfinal) = 0;
  virtual void Euler() = 0;
  virtual void Source() = 0;
  virtual void RightFlux() = 0;
  virtual void LeftFlux() = 0;
};

class Rusanov : public Schema_VF_1D
{
 protected:
   std::vector<std::vector<double> > _Fl; // flux gauche
   std::vector<std::vector<double> > _Fr; // flux droit

 public:
   void TimeScheme(tfinal);
   void Euler();
   void Source();
   void RightFlux();
   void LeftFlux();
};

class Relaxation : public Schema_VF_1D
{
 protected:
  std::vector<std::vector<double> > _Wd; // (rho, u, epsilon, pi, psi)
  std::vector<std::vector<double> > _Fdl; // flux gauche pour le problème de relaxation
  std::vector<std::vector<double> > _Fdr; // flux droit pour le problème de relaxation

 public:
   void TimeScheme(tfinal);
   void RightFlux();
   void LeftFlux();
};
