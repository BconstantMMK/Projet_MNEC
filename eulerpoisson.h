#ifndef EULERPOISSON_H_
#define EULERPOISSON_H_

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

  double _x_min, _x_max, _h_x, _dt;
  double _CL_u_g, _CL_u_d, _CL_rho_g, _CL_rho_d, _CL_E_g, _CL_E_d, _CL_phi_g, _CL_phi_d;
  int _Nx;

  std::vector<std::vector<double> > _Wsol; // (rho, rho*u, rho*E)
  std::vector<std::vector<double> > _Wsol_moins; // vecteur solution à l'itération précédente
  std::vector<double> _Gravity; // (phi)
  std::vector<std::vector<double> > _LapMat1D; // matrice creuse du laplacien

 public: // Méthodes et opérateurs de la classe
  Schema_VF_1D();
  // Constructeur : Initialiser xmin, xmax, Nx, hx, dt, _Wsol_0
  virtual ~Schema_VF_1D();

  void Initialize(xmin,xmax,Nx,hx,dt,CI_rho,CI_u,CI_E,CL_u_g,CL_u_d,CL_rho_g,CL_rho_d,CL_E_g,CL_E_d,CL_phi_g,CL_phi_d);
  void SaveSol(const std::string& name_file);
  void Poisson();
  virtual void TimeScheme(tfinal) = 0;
  virtual void Euler() = 0;
  virtual void Source() = 0;
  virtual void Flux() = 0;
};

class Rusanov : public Schema_VF_1D
{
 protected:
   std::vector<std::vector<double> > _Fg; // flux gauche
   std::vector<std::vector<double> > _Fd; // flux droit

 public:
   void TimeScheme(tfinal);
   void Euler();
   void Source();
   void Flux();
};

class Relaxation : public Schema_VF_1D
{
 protected:
  std::vector<std::vector<double> > _Wdelta; // (rho, u, epsilon, pi, psi)
  std::vector<std::vector<double> > _Fdg; // flux gauche pour le problème de relaxation
  std::vector<std::vector<double> > _Fdd; // flux droit pour le problème de relaxation

 public:
   void TimeScheme(tfinal);
   void Flux();
};
