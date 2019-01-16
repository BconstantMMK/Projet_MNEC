#include "fonctions1D.h"

std::vector<double> prodMVC(std::vector<std::vector<double>> A, std::vector<double> x, int Nx)
{
  std::vector<double> y(Nx);

  for (int i = 0; i < Nx; i++)
  {
    y[i] =
        A[0][i] * x[i-1]
      + A[1][i] * x[i]
      + A[2][i] * x[i+1];
  }

  return y;
}

double dot(std::vector<double> u, std::vector<double> v)
{
  //Calcul le produit scalaire entre 2 vecteurs.

  const int N = u.size();
  double y = 0.0;

  for (int i = 0; i < N; i++)
    y += u[i] * v[i];

  return y;
}

std::vector<double> CG(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double> x0, double err, int kmax, int nx)
{
  // Algorithme du gradient conjugué parallèle qui prend en argument uniquement des vecteurs locaux et renvoie un vecteur local.

  int k = 0;

  double norm_r, nr_carre, nr2_carre;
  std::vector<double> w, r, r_1, p, d, x;

  int N = nx;

  w.resize(N);
  r.resize(N);
  p.resize(N);
  d.resize(N);

  x = x0;

  w = prodMVC(A, x, nx);

  for (int i = 0; i < N; i++)
  {
    r[i] = b[i] - w[i];
  }
  p = r;

  nr_carre = dot(r, r);
  norm_r = sqrt(nr_carre); //On stocke ces deux variables puisqu'on s'en sert toutes les deux plusieurs fois dans la suite

  while ((norm_r > err) and (k < kmax))
  {
    d = prodMVC(A, p, nx);

    double alpha = nr_carre / dot(p, d);

    for (int i = 0; i < N; i++)
    {
      x[i] += alpha * p[i];
      r[i] -= alpha * d[i]; //rk devient r(k+1)
    }

    nr2_carre = dot(r, r); //norme de r(k+1)

    double beta = nr2_carre / nr_carre;

    for (int i = 0; i < N; i++)
    {
      p[i] = r[i] + beta * p[i];
    }

    nr_carre = nr2_carre;
    norm_r = sqrt(nr_carre);

    k++;
  }

  return x;
}
