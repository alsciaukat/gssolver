#include <iostream>
#include <iomanip>
#include "base.hpp"

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

int main() {
  // May use std::cin for N
  int N = 10 + 1;
  Grid grid(N);
  Array<double> psi(grid.N_r, std::vector<double>(grid.N_z, 1));
  Array<double> psi_p;
  do {
	psi_p = psi;

  } while (!is_converged(psi_p, psi));
}
