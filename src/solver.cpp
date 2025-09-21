#include <iostream>
#include <iomanip>
#include "base.hpp"
#include <cmath>
#include <numbers>

/*
      Solov'ev's Solution

      psi(r, z) := 0.5 (b + c0) R^2 z^2 + c0 R zeta(r) z^2 + 0.5 (a - c0) R^2
  zeta(r)^2 zeta(r) := (r^2 - R^2)/(2R)

      we use the boundary given by
      psi(r, z) = k
  with default parameters:
      R = 3, a = 0.6, b = 0.2, c0 = 0.1, k = 3
*/

double e_fix = 0.1;
double e_sor = 0.1;
double omega = 1.5; // relaxation coefficient
double A = 1, mu0 = 4*std::numbers::pi*1e-7, sigma = 0.5;
double Ip = 500e3, psi_bdry = 0.1;

double p_prime(double psi) {
  return 3*A*A/(4*mu0)*(sigma*sigma - sigma)/std::pow(psi, 7);
}

double gg_prime(double psi) {
  return A/(4*std::pow(psi, 3));
}

bool is_converged(Array<double> psi_p, Array<double> psi) {
  bool converged = true;
  for (int i = 0; i < psi.size(); i++) {
    for (int j = 0; j < psi[0].size(); j++) {
      if ((psi[i][j] - psi_p[i][j]) / psi[i][j] > e_fix)
        converged = false;
    }
  }
  return converged;
}

void update_F(const Grid& grid, Array<double>& psi, Array<double>& F) {
  for (int i = 0; i < grid.N_r; i++) {
	for (int j = 0; j < grid.N_z; j++) {
	  F[i][j] = -(mu0*grid.r[i]*grid.r[i]*p_prime(psi[i][j]) + gg_prime(psi[i][j]));
	}
  }
}

void update_sigma(const Grid& grid, const Array<double>& F, double& sigma) {
  double I = 0;

  for (int i = 0; i < grid.N_r; i++) {
	for (int j = 0; j < grid.N_z; j++) {
	  if (grid.boundary[i][j].domain != Domain::INT) continue;
	  I += F[i][j]/(mu0*grid.r[i])*grid.h*grid.h;
	}
  }

  sigma = Ip / (psi_bdry * I);
}

void normalize(const Grid& grid, Array<double>& psi) {
  // normalize psi between 0 at the center and 1 at the boundary
}

void sor(const Grid& grid, Array<double>& psi, Array<double>& F, double& sigma) {
  // find solution to Grad-Shafranov equaiton using
  // successive over relation with fixed F and sigma
}

int main() {
  // May use std::cin for N in the future
  int N = 10 + 1;
  Grid grid(N);
  Array<double> psi(grid.N_r, std::vector<double>(grid.N_z, 1));
  Array<double> psi_p;
  Array<double> F(grid.N_r, std::vector<double>(grid.N_z, 0));
  double sigma;

  do {
	psi_p = psi;
	update_F(grid, psi, F);
	update_sigma(grid, F, sigma);

    sor(grid, psi, F, sigma);
    normalize(grid, psi);
  } while (!is_converged(psi_p, psi));

  std::cout << p_prime(10) << std::endl;
}
