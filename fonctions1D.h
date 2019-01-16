#include<iostream>
#include<mpi.h>
#include<cmath>
#include<vector>

std::vector<double> prodMVC(std::vector<std::vector<double> > A, std::vector<double> x, int nx);

double dot(std::vector<double> u, std::vector<double> v);

std::vector<double> CG(std::vector<std::vector<double> > A, std::vector<double> b, std::vector<double> x0 , double err, int kmax, int nx);
