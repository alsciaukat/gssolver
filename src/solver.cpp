#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <ranges>
#include <string>
#include <chrono>
#include <utility>
#include <toml++/toml.hpp>

#include "base.hpp"
#include "dataio.hpp"

void update_F(Field<double> &F, Field<double> &psi, InitialCondition &cond, const Parameters &param) {
  const Grid &grid = psi.grid;
  const std::vector<double> &r = grid.r;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      F[i, j] = -(mu0 * r[i] * r[i] * cond.p_prime(psi[i, j]) +
                  cond.ff_prime(psi[i, j])) + 1e-100;
    }
  }
}

void update_sigma(double &sigma, Field<double> &F, const Parameters &param) {
  double I = 0;
  const Grid &grid = F.grid;
  const std::vector<double> &r = grid.r;

  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      if (grid.boundary[i][j].domain != Domain::INT)
        continue;
      I += F[i, j] / (mu0 * r[i]) * grid.h * grid.h;
    }
  }
  if (I == 0) sigma = 1;
  else sigma = param.I_p / (param.psi_bdry * I);
}

double get_Ip(Field<double> &F, double sigma) {
  // takes in unnormalized F
  double I = 0;
  const Grid &grid = F.grid;
  const std::vector<double> &r = grid.r;

  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      if (grid.boundary[i][j].domain != Domain::INT)
        continue;
      I += F[i, j] / (mu0 * r[i]) * grid.h * grid.h;
    }
  }
  return I;
}  

std::pair<int , int> normalize_old(Field<double> &psi, const Parameters &param) {
  const Grid &grid = psi.grid;
  // normalize psi between 0 at the minimum and 1 at the boundary
  double psi_min = INFINITY;
  int i_min = 0;
  int j_min = 0;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
	  if (grid.boundary[i][j].domain == Domain::EXT) continue;
      if (psi[i, j] < psi_min) {
        psi_min = psi[i, j];
        i_min = i;
        j_min = j;
      }
    }
  }
  int i = i_min;
  int j = j_min;
  if (grid.boundary[i][j].domain != Domain::INT) {
	throw std::logic_error("Minimum of psi is not interior");
  }
    
    // throw std::range_error("Minimum of psi is not interior");
  if (grid.boundary[i + 1][j + 1].domain == Domain::EXT ||
      grid.boundary[i + 1][j - 1].domain == Domain::EXT ||
      grid.boundary[i - 1][j + 1].domain == Domain::EXT ||
      grid.boundary[i - 1][j - 1].domain == Domain::EXT)
    std::cerr << "WARNING: The minimum is almost on the boundary. Do not trust "
                 "the result.\n"
              << std::endl;

  const std::vector<double> &r = grid.r;
  const std::vector<double> &z = grid.z;
  const double &h = grid.h;
  double r_min = grid.r[i];
  double z_min = grid.z[j]; // they are not minimum yet
  double r_delta;
  double z_delta;

  // finding the minimum point iteratively
  //
  // Succesive second-order univariate interpolation is used to
  // approximate the true minimum.
  // the solutions of the interpolation for r_min and z_min are rational
  // function of quartic polynomials, and analytic solution is not simple.
  //
  // r_min and z_min can be represented by each other, so this is iterating
  // between the two.  
  do {
    double r_min_p = r_min;
    double z_min_p = z_min;
    r_min =
        (r[i] + r[i - 1]) / 2 - h *
                                    (psi.interpolate_z(i, j, z_min) -
                                     psi.interpolate_z(i - 1, j, z_min)) /
                                    (psi.interpolate_z(i + 1, j, z_min) -
                                     2 * psi.interpolate_z(i, j, z_min) +
                                     psi.interpolate_z(i - 1, j, z_min));
    z_min =
        (z[j] + z[j - 1]) / 2 - h *
                                    (psi.interpolate_r(i, j, r_min) -
                                     psi.interpolate_r(i, j - 1, r_min)) /
                                    (psi.interpolate_r(i, j + 1, r_min) -
                                     2 * psi.interpolate_r(i, j, r_min) +
                                     psi.interpolate_r(i, j - 1, r_min));
    r_delta = r_min - r_min_p;
    z_delta = z_min - z_min_p;
  } while (std::abs(r_delta / r_min) > param.e_min ||
           std::abs(z_delta / z_min) > param.e_min);
  if (param.verbose) std::cout << "(r_min, z_min): (" << r_min << ", " << z_min << ")\n";

  double psi_min_p = psi_min;
  psi_min = psi.interpolate_z(i - 1, j, z_min) +
            (r_min - r[i - 1]) *
                (psi.interpolate_z(i, j, z_min) -
                 psi.interpolate_z(i - 1, j, z_min)) /
                h +
            (r_min - r[i - 1]) * (r_min - r[i]) *
                (psi.interpolate_z(i + 1, j, z_min) -
                 2 * psi.interpolate_z(i, j, z_min) +
                 psi.interpolate_z(i - 1, j, z_min)) /
                (2 * h * h);
  if (psi_min_p < psi_min)
    std::cerr << "WARNING: Minimum is not accurate. Might cause problem.\n";

  if (param.verbose)
	std::cout << "psi(r_min, z_min) = " << psi_min << "\n";

  // normalize ψ
  for (double &x : psi.value | std::views::join) {
    x = (x - psi_min) / (param.psi_l - psi_min) * param.psi_l;
  }
  return {i, j};
}

std::pair<int , int> normalize(Field<double> &psi, const Parameters &param, double &sigma) {
  const Grid &grid = psi.grid;
  // normalize psi between 0 at the minimum and 1 at the boundary
  double psi_min = INFINITY;
  int i_min = 0;
  int j_min = 0;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
	  if (grid.boundary[i][j].domain == Domain::EXT) continue;
      if (psi[i, j] < psi_min) {
        psi_min = psi[i, j];
        i_min = i;
        j_min = j;
      }
    }
  }
  int i = i_min;
  int j = j_min;
  if (grid.boundary[i][j].domain != Domain::INT) {
	throw std::logic_error("Minimum of psi is not interior");
  }
    
    // throw std::range_error("Minimum of psi is not interior");
  if (grid.boundary[i + 1][j + 1].domain == Domain::EXT ||
      grid.boundary[i + 1][j - 1].domain == Domain::EXT ||
      grid.boundary[i - 1][j + 1].domain == Domain::EXT ||
      grid.boundary[i - 1][j - 1].domain == Domain::EXT)
    std::cerr << "WARNING: The minimum is almost on the boundary. Do not trust "
                 "the result.\n"
              << std::endl;

  const std::vector<double> &r = grid.r;
  const std::vector<double> &z = grid.z;
  const double &h = grid.h;
  double r_min = grid.r[i];
  double z_min = grid.z[j]; // they are not minimum yet
  double r_delta;
  double z_delta;

  // finding the minimum point iteratively
  //
  // Succesive second-order univariate interpolation is used to
  // approximate the true minimum.
  // the solutions of the interpolation for r_min and z_min are rational
  // function of quartic polynomials, and analytic solution is not simple.
  //
  // r_min and z_min can be represented by each other, so this is iterating
  // between the two.  
  do {
    double r_min_p = r_min;
    double z_min_p = z_min;
    r_min =
        (r[i] + r[i - 1]) / 2 - h *
                                    (psi.interpolate_z(i, j, z_min) -
                                     psi.interpolate_z(i - 1, j, z_min)) /
                                    (psi.interpolate_z(i + 1, j, z_min) -
                                     2 * psi.interpolate_z(i, j, z_min) +
                                     psi.interpolate_z(i - 1, j, z_min));
    z_min =
        (z[j] + z[j - 1]) / 2 - h *
                                    (psi.interpolate_r(i, j, r_min) -
                                     psi.interpolate_r(i, j - 1, r_min)) /
                                    (psi.interpolate_r(i, j + 1, r_min) -
                                     2 * psi.interpolate_r(i, j, r_min) +
                                     psi.interpolate_r(i, j - 1, r_min));
    r_delta = r_min - r_min_p;
    z_delta = z_min - z_min_p;
  } while (std::abs(r_delta / r_min) > param.e_min ||
           std::abs(z_delta / z_min) > param.e_min);
  if (param.verbose) std::cout << "(r_min, z_min): (" << r_min << ", " << z_min << ")\n";

  double psi_min_p = psi_min;
  psi_min = psi.interpolate_z(i - 1, j, z_min) +
            (r_min - r[i - 1]) *
                (psi.interpolate_z(i, j, z_min) -
                 psi.interpolate_z(i - 1, j, z_min)) /
                h +
            (r_min - r[i - 1]) * (r_min - r[i]) *
                (psi.interpolate_z(i + 1, j, z_min) -
                 2 * psi.interpolate_z(i, j, z_min) +
                 psi.interpolate_z(i - 1, j, z_min)) /
                (2 * h * h);
  if (psi_min_p < psi_min)
    std::cerr << "WARNING: Minimum is not accurate. Might cause problem.\n";

  if (param.verbose)
	std::cout << "psi(r_min, z_min) = " << psi_min << "\n";

  // normalize ψ
  for (double &x : psi.value | std::views::join) {
    x = (x - psi_min) / (param.psi_l - psi_min) * param.psi_l;
  }
  sigma = std::pow( param.psi_l / (param.psi_l - psi_min), 2);
  return {i, j};
}

std::tuple<int, int, int, int, double, double> get_min_max(Field<double> &psi, const Parameters &param) {
  const Grid &grid = psi.grid;
  // normalize psi between 0 at the minimum and 1 at the boundary
  double psi_max = -INFINITY;
  double psi_min = INFINITY;
  int i_min = 0;
  int j_min = 0;
  int i_max = 0;
  int j_max = 0;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      if (psi[i, j] < psi_min) {
        psi_min = psi[i, j];
        i_min = i;
        j_min = j;
      }
      if (psi[i, j] > psi_max) {
		psi_max = psi[i, j];
        i_max = i;
        j_max = j;
	  }        
    }
  }
  if (grid.boundary[i_min][j_min].domain == Domain::EXT) {
	std::cerr << "WARNING: Minimum of psi is exterior.\n";
  }
  return {i_min, j_min, i_max, j_max, psi_min, psi_max};
}

std::tuple<int, int, int, int, double, double> get_range(Field<double> &psi,
                                              const Parameters &param) {
  const Grid &grid = psi.grid;
  auto [i_min, j_min, _, __, psi_min, psi_max] = get_min_max(psi, param);
  if (psi_max - param.psi_l > 1e-25)
    throw std::logic_error("The maximum is not psi_l");
  
  int i_max = i_min;
  bool i_max_set = false;
  for (int i = i_min; i < grid.N_r; i++) {
    if (psi[i, j_min] < psi_max)
      continue;
    else {
      i_max = i;
	  i_max_set = true;
      break;
    }
  }
  if (!i_max_set) throw std::logic_error("The maximum is not on the axis.");
  return {i_min, j_min, i_max, j_min, psi_min, psi_max};
}

Field<double> get_normalized(Field<double> &psi, const Parameters &param) {
  const Grid &grid = psi.grid;
  auto [i_min, j_min, _, __, psi_min, psi_max] = get_min_max(psi, param);
  if (psi_max - param.psi_l > 1e-25)
    throw std::logic_error("The maximum is not psi_l");
  Field<double> psi_N(grid, param, 1);
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      psi_N[i, j] = (psi[i, j] - psi_min) / (psi_max - psi_min);
    }
  }
  return psi_N;
}

void update_F_unnormalized(Field<double> &F, Field<double> &psi, InitialCondition &cond, const Parameters &param, const double &sigma) {
  const Grid &grid = psi.grid;
  const std::vector<double> &r = grid.r;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      F[i, j] = -(mu0 * r[i] * r[i] * cond.p_prime(psi[i, j]) +
                  cond.ff_prime(psi[i, j])) + 1e-100;
	  F[i, j] *= std::sqrt(std::abs(sigma));
    }
  }
}

void denormalize(Field<double> &psi, double sigma) {
  Grid grid = psi.grid;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
	  psi[i, j] = (psi[i, j] - 1) / std::sqrt(std::abs(sigma)) + 1;
	}
  }
}


void sor_old(Field<double> &psi, Field<double> &F,
         const double &sigma, const Parameters &param) {
  // compute psi_delta with less precise algorithm
  const Grid &grid = psi.grid;
  const Array<BoundaryInfo> &boundary = grid.boundary;
  const std::vector<double> &r = grid.r;
  double h = grid.h;

  bool converged;
  int n_it = 0;
  do {
    converged = true;
    for (int i = 0; i < grid.N_r; i++) {
      for (int j = 0; j < grid.N_z; j++) {
        if (boundary[i][j].domain == Domain::EXT)
          continue;
        double psi_delta;
        if (boundary[i][j].domain == Domain::INT) {
          psi_delta =
			((1 - h / (2*r[i])) * psi[i + 1, j] + (1 + h / (2*r[i])) * psi[i - 1, j] -
               4 * psi[i, j] + psi[i, j + 1] + psi[i, j - 1] -
               h * h * sigma * F[i, j]) /
              4;
        } else if (boundary[i][j].domain == Domain::BDRY) {
          double alpha1 = boundary[i][j].alpha1;
          if (alpha1 < 1)
            psi[i - 1, j] = param.psi_l;
          double alpha2 = boundary[i][j].alpha2;
          if (alpha2 < 1)
            psi[i + 1, j] = param.psi_l;
          double beta1 = boundary[i][j].beta1;
          if (beta1 < 1)
            psi[i, j - 1] = param.psi_l;
          double beta2 = boundary[i][j].beta2;
          if (beta2 < 1)
            psi[i, j + 1] = param.psi_l;
          double denom_a = alpha1 * alpha2 * (alpha1 + alpha2);
          double denom_b = beta1 * beta2 * (beta1 + beta2);
          double denom_t =
              2 * (alpha1 + alpha2) / denom_a + 2 * (beta1 + beta2) / denom_b -
              (alpha1 * alpha1 - alpha2 * alpha2) * h / (denom_a * r[i]);
          psi_delta =
              (((2 * alpha2 + alpha2 * alpha2 * h / r[i]) * psi[i - 1, j] +
                (2 * alpha1 - alpha1 * alpha1 * h / r[i]) * psi[i + 1, j] +
                (-2 * (alpha1 + alpha2) +
                 (alpha1 * alpha1 - alpha2 * alpha2) * h / r[i]) *
                    psi[i, j]) /
                   denom_a +
               (2 * beta2 * psi[i, j - 1] - 2 * (beta1 + beta2) * psi[i, j] +
                2 * beta1 * psi[i, j + 1]) /
                   denom_b -
               h * h * sigma * F[i, j]) /
              denom_t;
        } else
          throw std::invalid_argument("Domain is not known");
        psi[i, j] = psi[i, j] + param.omega * psi_delta;

        if (std::abs(psi_delta / psi[i, j]) > param.e_sor)
          converged = false;
      }
    }
  } while (!converged);
}


void sor(Field<double> &psi, Field<double> &F,
         const double &sigma, const Parameters &param) {
  // find solution to Grad-Shafranov equaiton using
  // successive over relation with fixed right-hand side evaluated
  // at the ψ from the previous step.  
  const Grid &grid = psi.grid;
  const Array<BoundaryInfo> &boundary = grid.boundary;
  const std::vector<double> &r = grid.r;
  double h = grid.h;

  bool converged;
  int n_it = 0;
  do {
    converged = true;
    for (int i = 0; i < grid.N_r; i++) {
      for (int j = 0; j < grid.N_z; j++) {
        if (boundary[i][j].domain == Domain::EXT)
          continue;
        double psi_delta;
        if (boundary[i][j].domain == Domain::INT) {
		  double a = 1 / (1 + h / (2 * r[i]));
		  double b = 1 / (1 - h / (2 * r[i]));
		  double c = a + b + 2;
		  psi_delta =
			  (a * psi[i + 1, j] + b * psi[i - 1, j] - c * psi[i, j] +
			   psi[i, j + 1] + psi[i, j - 1] - h * h * sigma * F[i, j]) /
			  c;
        } else if (boundary[i][j].domain == Domain::BDRY) {
          double alpha1 = boundary[i][j].alpha1;
          if (alpha1 < 1)
            psi[i - 1, j] = param.psi_l;
          double alpha2 = boundary[i][j].alpha2;
          if (alpha2 < 1)
            psi[i + 1, j] = param.psi_l;
          double beta1 = boundary[i][j].beta1;
          if (beta1 < 1)
            psi[i, j - 1] = param.psi_l;
          double beta2 = boundary[i][j].beta2;
          if (beta2 < 1)
            psi[i, j + 1] = param.psi_l;
          double denom_a = alpha1 * alpha2 * (alpha1 + alpha2);
          double denom_b = beta1 * beta2 * (beta1 + beta2);
          double denom_t =
              2 * (alpha1 + alpha2) / denom_a + 2 * (beta1 + beta2) / denom_b -
              (alpha1 * alpha1 - alpha2 * alpha2) * h / (denom_a * r[i]);
          psi_delta =
              (((2 * alpha2 + alpha2 * alpha2 * h / r[i]) * psi[i - 1, j] +
                (2 * alpha1 - alpha1 * alpha1 * h / r[i]) * psi[i + 1, j] +
                (-2 * (alpha1 + alpha2) +
                 (alpha1 * alpha1 - alpha2 * alpha2) * h / r[i]) *
                    psi[i, j]) /
                   denom_a +
               (2 * beta2 * psi[i, j - 1] - 2 * (beta1 + beta2) * psi[i, j] +
                2 * beta1 * psi[i, j + 1]) /
                   denom_b -
               h * h * sigma * F[i, j]) /
              denom_t;
        } else
          throw std::invalid_argument("Domain is not known");
        psi[i, j] = psi[i, j] + param.omega * psi_delta;

        if (std::abs(psi_delta / psi[i, j]) > param.e_sor)
          converged = false;
      }
    }
  } while (!converged);
}

void sor_unnormalized(Field<double> &psi, Field<double> &F,
         const Parameters &param) {
  // assume that psi and F are all unnormalized and just solves it.
  const Grid &grid = psi.grid;
  const Array<BoundaryInfo> &boundary = grid.boundary;
  const std::vector<double> &r = grid.r;
  double h = grid.h;

  bool converged;
  int n_it = 0;
  do {
    converged = true;
    for (int i = 0; i < grid.N_r; i++) {
      for (int j = 0; j < grid.N_z; j++) {
        if (boundary[i][j].domain == Domain::EXT)
          continue;
        double psi_delta;
        if (boundary[i][j].domain == Domain::INT) {
		  double a = 1 / (1 + h / (2 * r[i]));
		  double b = 1 / (1 - h / (2 * r[i]));
		  double c = a + b + 2;
		  psi_delta =
			  (a * psi[i + 1, j] + b * psi[i - 1, j] - c * psi[i, j] +
			   psi[i, j + 1] + psi[i, j - 1] - h * h * F[i, j]) /
			  c;
        } else if (boundary[i][j].domain == Domain::BDRY) {
          double alpha1 = boundary[i][j].alpha1;
          if (alpha1 < 1)
            psi[i - 1, j] = param.psi_l;
          double alpha2 = boundary[i][j].alpha2;
          if (alpha2 < 1)
            psi[i + 1, j] = param.psi_l;
          double beta1 = boundary[i][j].beta1;
          if (beta1 < 1)
            psi[i, j - 1] = param.psi_l;
          double beta2 = boundary[i][j].beta2;
          if (beta2 < 1)
            psi[i, j + 1] = param.psi_l;
          double denom_a = alpha1 * alpha2 * (alpha1 + alpha2);
          double denom_b = beta1 * beta2 * (beta1 + beta2);
          double denom_t =
              2 * (alpha1 + alpha2) / denom_a + 2 * (beta1 + beta2) / denom_b -
              (alpha1 * alpha1 - alpha2 * alpha2) * h / (denom_a * r[i]);
          psi_delta =
              (((2 * alpha2 + alpha2 * alpha2 * h / r[i]) * psi[i - 1, j] +
                (2 * alpha1 - alpha1 * alpha1 * h / r[i]) * psi[i + 1, j] +
                (-2 * (alpha1 + alpha2) +
                 (alpha1 * alpha1 - alpha2 * alpha2) * h / r[i]) *
                    psi[i, j]) /
                   denom_a +
               (2 * beta2 * psi[i, j - 1] - 2 * (beta1 + beta2) * psi[i, j] +
                2 * beta1 * psi[i, j + 1]) /
                   denom_b -
               h * h * F[i, j]) /
              denom_t;
        } else
          throw std::invalid_argument("Domain is not known");
        psi[i, j] = psi[i, j] + param.omega * psi_delta;

        if (std::abs(psi_delta / psi[i, j]) > param.e_sor)
          converged = false;
      }
    }
  } while (!converged);
}

double get_max_error(Field<double> &a, Field<double> &b) {
  double N_r = a.value.size();
  double N_z = a.value[0].size();
  double error = 0;
  const Grid &grid = a.grid;
  for (int i = 0; i < N_r; i++) {
    for (int j = 0; j < N_z; j++) {
	  if (grid.boundary[i][j].domain != Domain::INT)
		continue;
	  double error_tmp = std::abs((a[i, j] - b[i, j]) / a[i, j]);
	  if (error_tmp > error) error = error_tmp;
    }
  }
  return error;
}  

double get_l2_error(Field<double> &a, Field<double> &b) {
  double N_r = a.value.size();
  double N_z = a.value[0].size();
  double error_sq_sum = 0;
  const Grid &grid = a.grid;
  for (int i = 0; i < N_r; i++) {
    for (int j = 0; j < N_z; j++) {
	  if (grid.boundary[i][j].domain != Domain::INT)
		continue;
	  error_sq_sum += std::pow((a[i, j] - b[i, j]), 2) * grid.h * grid.h;
    }
  }
  return std::sqrt(error_sq_sum);
}  

void print_array(const Array<double> &a) {
  for (auto &b : a) {
    for (auto &c : b) {
      std::cout << c << "\t";
    }
    std::cout << "\n";
  }
}

std::string to_lower(std::string str) {
  auto str_v = str | std::views::transform([](auto c){ return static_cast<char>(std::tolower(c)); });
  std::string str_l(str_v.begin(), str_v.end());
  return str_l;  
}  

#define SET(SEC, VAR, TYPE, IVAR, DEFAULT) param.IVAR = config[#SEC][#VAR].value<TYPE>().value_or(DEFAULT);

int initialize(Parameters &param) {
  toml::table config = toml::parse_file("config.toml");
  param.select = config["select"].value<std::string>().value_or("default");
  SET(grid, type, std::string, gtype, "solovev")
  SET(grid, N, int, N, 11)
  SET(grid, R, double, R, 1.8)
  SET(grid, a, double, a, 0.5)
  SET(grid, b, double, b, 0.5)
  SET(grid, c0, double, c0, 0.1)
  SET(grid, k, double, k, 0.11)
  SET(grid, tolerance, double, tolerance, 0.1)
#ifdef USE_NETCDF
  SET(output, format, std::string, offormat, "netCDF")
  SET(output, name, std::string, ofname, "result.nc")
#else
  SET(output, format, std::string, offormat, "csv")
  SET(output, name, std::string, ofname, "result.csv")
#endif
  SET(output, print_psi, bool, print, false)
  SET(output, verbose, bool, verbose, false)
  SET(output, print_every, int, print_every, 50)
  SET(output, N_gr, int, N_gr, 17)
  SET(output, logh_low, double, logh_low, -2.4)
  SET(output, delta_logh, double, delta_logh, 0.1)
  SET(initial_condition, type, std::string, ictype, "solovev");
  SET(initial_condition, psi_start, double, psi_start, 0.4);
  SET(initial_condition, psi_end, double, psi_end, 0.7);
  SET(initial_condition, beta0, double, beta0, 0.5);
  SET(initial_condition, m, int, m, 2);
  SET(initial_condition, n, int, n, 2);
  SET(initial_condition, p0, double, p0, 1);
  SET(initial_condition, p_a, double, p_a, 1);
  SET(initial_condition, B0, double, B0, 1);
  SET(initial_condition, f_a, double, f_a, 1);
  SET(initial_condition, f_b, double, f_b, 1);
  SET(initial_condition, p, int, p, 2);
  SET(initial_condition, q, int, q, 2);
  SET(initial_condition, r, int, r, 2);
  SET(initial_condition, s, int, s, 2);
  SET(initial_condition, para, bool, para, false);
  SET(solver, e_fix, double, e_fix, 0.01);
  SET(solver, e_sor, double, e_sor, 0.01);
  SET(solver, e_min, double, e_min, 0.01);
  SET(solver, omega, double, omega, 1.9);
  SET(solver, I_p, double, I_p, 500e3);
  SET(solver, psi_l, double, psi_l, 1);
  SET(solver, psi_bdry, double, psi_bdry, 0.1);
  SET(solver, sigma0, double, sigma0, 1e8);
  return 0;
}

std::unique_ptr<InitialCondition>
get_initial_condition(const Parameters &param) {
  // *factory pattern*
  // this function returns a pointer to
  // a concrete class of `InitialCondition`.
  std::string ictype_l = to_lower(param.ictype);
  if (ictype_l == "polynomial") {
	std::cout << "Using Polynomial Profile" << std::endl;
    return std::make_unique<PolynomialCondition>(param);
  }
  if (ictype_l == "hmode") {
	std::cout << "Using H-Mode Profile" << std::endl;
    return std::make_unique<HModeCondition>(param);
  }    
  if (ictype_l == "solovev") {
	std::cout << "Using Solovev Profile" << std::endl;
    return std::make_unique<SolovevCondition>(param);
  }    
  if (ictype_l == "diamagnetic") {
	std::cout << "Using Diamagnetic Profile" << std::endl;
    return std::make_unique<DiamagneticCondition>(param);
  }    
  std::cout << "Initial condition type not known. Defaulting to Solovev Condition." << std::endl;
  return std::make_unique<SolovevCondition>(param);
}

std::pair<int, int> solve_default(const Parameters &param, const Grid &grid, Field<double> &psi, Field<double> &F, double &sigma, InitialCondition &cond) {
  // assumes the Ip and psi_bdry in order to estimate the sigma required for
  // normalization. They are wrong but it helps the convergence.  
  int n_fix = 0;
  Field<double> psi_p(grid, param, 1);

  double max_error;  
  double l2_error;
  std::pair<int, int> min_index{grid.N_r/2, grid.N_z/2};

  do {
	// The Grad-Shafranov can be rewritten:
	//     Δ*ψ = σ F(ψ)
	// where Δ* is the elliptic toroidal operator that is
	// simply equivalent to the left-hand side of Grad-Shafranov equation.
	//
	// Two nested iteration is used to solve this:
	// outer iteration is a fixed-point iteration where the right-hand side
	// is evaluated by ψ⁽ⁿ⁾ and the equation is solved to obtain ψ⁽ⁿ⁺¹⁾.
	//     Δ*ψ = σ F(ψ⁽ⁿ⁾)
	// This elliptic partial differential equation is then solved using
	// the inner iteration called succesive over relaxation (SOR) method.
	//
	// ψ stays normalized between 0 and 1, with the aid of
	// appropriate scaling by σ which is calculated such that the toroidal
	// current I_p is kept constant.
	// See /Accurate calculation of the gradients of the equilibrium poloidal
	// flux in tokamaks/ by M. Woo, G. Jo, B. H. Park, A. Y. Aydemir, J. H. Kim.
	n_fix++;
	if (param.verbose) {
		std::cout << "\nFixed-point Iteration: " << n_fix << "\n";
		std::cout << "Sigma: " << sigma << "\n";
	}

	update_F(F, psi, cond, param);
	update_sigma(sigma, F, param);
    psi_p = psi;
	min_index = normalize_old(psi, param);

    sor(psi, F, sigma, param);

	// max_error = get_max_error(psi, psi_p);
	l2_error = get_l2_error(psi, psi_p);
	if (param.verbose) {
		std::cout << "L2 Error: " << l2_error << "\n";
	} else {
	  std::cout << "# iter: " << n_fix << "\t" << "L2 error: " << l2_error << "\r";
	}          
  } while (l2_error > param.e_fix);
  std::cout << "\n";
  min_index = normalize(psi, param, sigma);
  return min_index;
}

std::pair<int, int> solve_default_unnormalized_output(const Parameters &param, const Grid &grid, Field<double> &psi, Field<double> &F, double &sigma, InitialCondition &cond) {
  // do the same thing as `solve_default`, but only outputs unnormalized psi.
  int n_fix = 0;
  Field<double> psi_p(grid, param, 1);

  double max_error;  
  double l2_error;
  std::pair<int, int> min_index{grid.N_r/2, grid.N_z/2};

  do {
	n_fix++;
	if (param.verbose) {
		std::cout << "\nFixed-point Iteration: " << n_fix << "\n";
		std::cout << "Sigma: " << sigma << "\n";
	}

	update_F(F, psi, cond, param);
	update_sigma(sigma, F, param);
    psi_p = psi;
	min_index = normalize_old(psi, param);

    sor(psi, F, sigma, param);

	l2_error = get_l2_error(psi, psi_p);
	if (param.verbose) {
		std::cout << "L2 Error: " << l2_error << "\n";
	} else {
	  std::cout << "# iter: " << n_fix << "\t" << "L2 error: " << l2_error << "\r" << std::flush;
	}          
  } while (l2_error > param.e_fix);
  std::cout << "\n";
  return min_index;
}

std::pair<int, int> solve_direct(const Parameters &param, const Grid &grid, Field<double> &psi, Field<double> &F, InitialCondition &cond) {
  // directly solve grad-shafranov without normalization.
  // CAUTION: This only works when the initial condition is given in terms of unnormalized psi.
  int n_fix = 0;
  Field<double> psi_p(grid, param, 1);

  double max_error;  
  double l2_error;
  std::pair<int, int> min_index{grid.N_r/2, grid.N_z/2};
  do {
	n_fix++;
	if (param.verbose) {
		std::cout << "\nFixed-point Iteration: " << n_fix << "\n";
	}
	update_F(F, psi, cond, param);

    psi_p = psi;

    sor_unnormalized(psi, F, param);

	l2_error = get_l2_error(psi, psi_p);
	if (param.verbose) {
		std::cout << "L2 Error: " << l2_error << "\n";
	} else {
	  std::cout << "# iter: " << n_fix << "\t" << "L2 error: " << l2_error << "\r";
	}          
  } while (l2_error > 1e-15);
  std::cout << "\n";
  update_F(F, psi, cond, param);
  return min_index;
}

void solve_unnormalized(const Parameters &param, const Grid &grid, Field<double> &psi, Field<double> &F, double &sigma, InitialCondition &cond) {
  // Take unnormalized psi and F as a function of normalized psi,
  // and calculate in the unnormalized realm.
  //
  // The value of F does not matter, since it is recomputed anyway.
  //
  // Outpus unnormalized psi and F.  

  int n_fix = 0;
  Field<double> psi_p(grid, param, 1);
  Field<double> psi_N(grid, param, 1);

  double max_error;  
  double l2_error;
  std::pair<int, int> min_index{grid.N_r/2, grid.N_z/2};

  do {
	n_fix++;
	if (param.verbose) {
		std::cout << "\nFixed-point Iteration: " << n_fix << "\n";
	}
    auto [i_min, j_min, _, __, psi_min, psi_max] = get_min_max(psi, param);
	if (param.verbose)
	std::cout << "psi_min: " << psi_min << "\npsi_max: " << psi_max << "\n";

	// normalize
    for (int i = 0; i < grid.N_r; i++) {
  	  for (int j = 0; j < grid.N_z; j++) {
		psi_N[i, j] = (psi[i, j] - psi_min) / (psi_max - psi_min);
	  }
	}

	update_F_unnormalized(F, psi_N, cond, param, 1/std::pow(psi_max - psi_min, 2));

	psi_p = psi;
    sor_unnormalized(psi, F, param);

	l2_error = get_l2_error(psi, psi_p);
	if (param.verbose) {
		std::cout << "L2 Error: " << l2_error << "\n";
	} else if (n_fix % param.print_every == 0) {
	  std::cout << "# iter: " << n_fix << "\t" << "L2 error: " << l2_error << "\r" << std::flush;
	}          
  } while (l2_error > param.e_fix);
  std::cout << "\n";
}

std::pair<int, int> solve_polynomial(const Parameters &param, const Grid &grid, Field<double> &psi, Field<double> &F, double &sigma, InitialCondition &cond, std::vector<double> &max_errors, std::vector<double> &l2_errors, netCDF::NcFile &file) {
  // The code originally used for testing polynomial initial condition, and measure its convergence rate.
  if (max_errors.size() > 0 || l2_errors.size() > 0) throw std::invalid_argument("The output arguments `max_errors` and `l2_errors` must be empty");
  // initializations
  update_F(F, psi, cond, param);
  update_sigma(sigma, F, param);

  int n_fix = 0;
  Field<double> psi_p(grid, param, 1);
  double max_error;
  double l2_error;
  std::pair<int, int> min_index{grid.N_r/2, grid.N_z/2};

  do {
	n_fix++;
	std::cout << "\nFixed-point Iteration: " << n_fix << "\n";

    psi_p = psi;
    std::cout << "Sigma: " << sigma << "\n";

    sor(psi, F, sigma, param);
	min_index = normalize(psi, param, sigma);

    update_F(F, psi, cond, param);
    update_sigma(sigma, F, param);
	max_error = get_max_error(psi, psi_p);
	l2_error = get_l2_error(psi, psi_p);
	max_errors.push_back(max_error);
	l2_errors.push_back(l2_error);        

    std::string i_str = std::to_string(n_fix);
    store_field(file, psi, "psi" + i_str, {"r" + i_str, "z" + i_str});
    store_vector(file, grid.r, "r" + i_str, grid.N_r, "r" + i_str);
    store_vector(file, grid.z, "z" + i_str, grid.N_z, "z" + i_str);
    store_field(file, F, "F" + i_str, {"r" + i_str, "z" + i_str});
    std::cout << "L2 Error: " << l2_error << "\n";        
  } while (l2_error > param.e_fix);
  return min_index;
}

void store_default(netCDF::NcFile &file, const Parameters &param, Field<double> &psi, Field<double> &F, double &sigma) {
  // output to a file 
  std::string offormat_l = to_lower(param.offormat);

  if (offormat_l == "netcdf") {
#ifdef USE_NETCDF
    store_netcdf(file, psi, F, sigma, param);
	std::cout << "Wrote to: " << param.ofname << std::endl;
#else
    throw std::invalid_argumet("Output to NetCDF is not supported. Change output.format to csv");
#endif
  } else if (offormat_l == "csv") {
	store_csv(psi, param.ofname);
	std::cout << "Wrote to: " << param.ofname << std::endl;
  } else {
	throw std::invalid_argument("No output is written. Current output.format is not supported.");
  }    
  if (param.print) print_array(psi.value);
}

void run_solovev(Parameters &param) {
  double sigma = 1;
  auto cond = get_initial_condition(param);
  Grid grid(param);

  std::vector<double> max_errors;
  std::vector<double> l2_errors;
  std::vector<double> time2run;

  int N_grid_analysis = param.N_gr;
  double logh_low = param.logh_low;
  double delta_logh = param.delta_logh;
  std::vector<double> h;
  std::vector<int> N;

  netCDF::NcFile file(param.ofname, netCDF::NcFile::replace);

  for (int i = 0; i < N_grid_analysis; i++) {
	h.push_back(std::pow(10, logh_low + i * delta_logh));
	N.push_back(static_cast<int>((grid.r[grid.N_r - 1] - grid.r[0]) / h[i]));
	param.N = N[i];
	Grid grid(param);
	Field<double> psi(grid, param, 0.11);
	Field<double> F(grid, param, 0);
	Field<double> psi_solovev(grid, param, 1);
	for (int i = 0; i < grid.N_r; i++) {
	  for (int j = 0; j < grid.N_z; j++) {
		psi_solovev[i, j] = grid.solovev(grid.r[i], grid.z[j]);
	  }
	}

    auto start = std::chrono::high_resolution_clock::now();
    update_F(F, psi, *cond, param);
    sor(psi, F, sigma, param);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    time2run.push_back(duration.count());

    std::string i_str = std::to_string(i);
    store_field(file, psi, "psi" + i_str, {"r" + i_str, "z" + i_str});
    store_vector(file, grid.r, "r" + i_str, grid.N_r, "r" + i_str);
    store_vector(file, grid.z, "z" + i_str, grid.N_z, "z" + i_str);
    store_field(file, F, "F" + i_str, {"r" + i_str, "z" + i_str});
    store_field(file, psi_solovev, "psi_t" + i_str, {"r" + i_str, "z" + i_str});

    double error_max = 0;
    double error_sq_cum = 0;
    for (int i = 0; i < grid.N_r; i++) {
      for (int j = 0; j < grid.N_z; j++) {
        if (grid.boundary[i][j].domain != Domain::INT)
          continue;
        double error_tmp =
            std::abs((psi[i, j] - psi_solovev[i, j]) / psi[i, j]);
        error_sq_cum += std::pow(psi[i, j] - psi_solovev[i, j], 2) * grid.h * grid.h;
        if (error_tmp > error_max)
          error_max = error_tmp;
      }
    }
    l2_errors.push_back(std::sqrt(error_sq_cum));
    max_errors.push_back(error_max);
    std::cout << "h: " << h[i] << "\n";
    std::cout << "N: " << N[i] << "\n";
    std::cout << "Max Error: " << max_errors[i] << "\n";
  }
  std::vector<double> N_f(N.size());
  for (int i = 0; i < N.size(); i++) {
	N_f[i] = static_cast<double>(N[i]);
  }

  store_vector(file, h, "h", N_grid_analysis, "h");
  store_vector(file, N_f, "N", N_grid_analysis, "h");
  store_vector(file, max_errors, "max_error", N_grid_analysis, "h");
  store_vector(file, l2_errors, "l2_error", N_grid_analysis, "h");
  store_vector(file, time2run, "time2run", N_grid_analysis, "h");
}

void run_polynomial(Parameters &param) {
  double sigma = 1;
  auto cond = get_initial_condition(param);
  Grid grid(param);
  netCDF::NcFile file(param.ofname, netCDF::NcFile::replace);

  Field<double> psi(grid, param, 0.11);
  Field<double> F(grid, param, 0);
  std::vector<double> max_errors;
  std::vector<double> l2_errors;

  solve_polynomial(param, grid, psi, F, sigma, *cond, max_errors, l2_errors, file);

  store_field(file, psi, "psi", {"r", "z"});
  store_field(file, F, "F", {"r", "z"});
  store_vector(file, grid.r, "r", grid.N_r, "r");
  store_vector(file, grid.z, "z", grid.N_z, "z");
  store_vector(file, max_errors, "max_error", max_errors.size(), "i");
  store_vector(file, l2_errors, "l2_error", l2_errors.size(), "i");
}

double get_f(double psi, std::vector<double> &psi_range,
             std::vector<double> &f) {
  int left = 0, right = psi_range.size() - 1;
  if (psi > psi_range[right] || psi < psi_range[left]) throw std::range_error("Get G: Out of range");
  while (left != right) {
	int middle = (left + right) / 2 + 1;
	if (psi_range[middle] > psi) {
	  right = middle - 1;
	} else {
	  left = middle;
	}
  }
  return f[left];
}

void calculate_output_old(netCDF::NcFile &file, Field<double> &J_r, Field<double> &J_phi, Field<double> &J_z, Field<double> &B_r, Field<double> &B_phi, Field<double> &B_z, Field<double> &psi, Field<double> &F, const InitialCondition *cond, const Parameters &param) {
    // absolute value is just for debugging. Sigma is required to be positive.
  const Grid &grid = J_r.grid;

  auto [i_min, j_min, _, __, psi_min, psi_max] = get_min_max(psi, param);
	if (psi_max - param.psi_l > 1e-25) throw std::logic_error("The maximum is not psi_l");
	double psi_bdry = psi_max - psi_min;
    double sigma = 1 / (psi_bdry * psi_bdry);
	double I = get_Ip(F, sigma);

	std::cout << "Sigma: " << sigma << "\n";
	std::cout << "I: " << I << "\n";
    int i_max = 0;
    for (int i = i_min; i < grid.N_r; i++) {
      if (psi[i, j_min] < psi_max) continue;
      else {
        i_max = i;
		break;
      }
    }
    int N_range = i_max - i_min + 1;
    std::vector<double> psi_range(N_range, 0);
    std::vector<double> psi_range_N(N_range, 0);
    std::vector<double> B_phi_range(N_range, 0);
    std::vector<double> p_prime(N_range, 0);
    std::vector<double> ff_prime(N_range, 0);
    std::vector<double> f_sq(N_range, 0);
    std::vector<double> f(N_range, 0);
    std::vector<double> f_prime(N_range, 0);
    std::vector<double> p(N_range, 0);

    for (int i = 0; i < N_range; i++) {
      int ri = i_min + i;
      int zi = j_min;
      psi_range[i] = psi[ri, zi];
	  psi_range_N[i] = (psi_range[i] - psi_min) / (psi_max - psi_min);
      p_prime[i] = (*cond).p_prime(psi_range_N[i]) * psi_bdry;
      ff_prime[i] = (*cond).ff_prime(psi_range_N[i]) * psi_bdry;
    }
	p[N_range - 1] = -p_prime[N_range - 1] *
						(psi_range[N_range - 1] - psi_range[N_range - 2]);
	f_sq[N_range - 1] = -2 * ff_prime[N_range - 1] *
						(psi_range[N_range - 1] - psi_range[N_range - 2]);
    for (int i = N_range - 2; i >= 0; i--) {
        p[i] = p[i + 1] - p_prime[i] * (psi_range[i + 1] - psi_range[i]);

        f_sq[i] =
            f_sq[i + 1] - 2 * ff_prime[i] * (psi_range[i + 1] - psi_range[i]);
    }
    double f_sq_min = *std::min_element(f_sq.begin(), f_sq.end());
    double p_min = *std::min_element(p.begin(), p.end());
    if (f_sq_min < 0) {
      for (int i = 0; i < N_range; i++) {
        f_sq[i] += - 2 * f_sq_min;
      }
    }
    if (p_min < 0) {
      for (int i = 0; i < N_range; i++) {
        p[i] += - 2 * p_min;
      }
    }
    for (int i = N_range - 1; i >= 0; i--) {
        f[i] = std::sqrt(f_sq[i]);
        f_prime[i] = ff_prime[i] / f[i];		
    }
    std::cout << "p and f computed\n";
	for (int i = 0; i < grid.N_r; i++) {
      for (int j = 0; j < grid.N_z; j++) {
		if (grid.boundary[i][j].domain != Domain::INT) continue;
        J_phi[i, j] = -F[i, j] / (grid.r[i] * mu0);
		B_phi[i, j] = get_f(psi[i, j], psi_range, f) / grid.r[i];
        double phi_r = (psi[i + 1, j] - psi[i - 1, j]) / (2 * grid.h);
        double phi_z = (psi[i, j + 1] - psi[i, j - 1]) / (2 * grid.h);
		B_r[i, j] = - phi_z / grid.r[i];
		B_z[i, j] = phi_r / grid.r[i];
		J_r[i, j] = - phi_z * get_f(psi[i, j], psi_range, f_prime) / (mu0 * grid.r[i]);
		J_z[i, j] = phi_r * get_f(psi[i, j], psi_range, f_prime) / (mu0 * grid.r[i]);
      }
    }
	std::cout << "B and J computed\n";
	store_vector(file, grid.r, "r", grid.N_r, "r");
	store_vector(file, grid.z, "z", grid.N_z, "z");
	store_field(file, psi, "psi", {"r", "z"});
	store_field(file, F, "F", {"r", "z"});
	store_field(file, J_r, "J_r", {"r", "z"});
	store_field(file, J_phi, "J_phi", {"r", "z"});
	store_field(file, J_z, "J_z", {"r", "z"});
	store_field(file, B_r, "B_r", {"r", "z"});
	store_field(file, B_phi, "B_phi", {"r", "z"});
	store_field(file, B_z, "B_z", {"r", "z"});

	store_vector(file, psi_range, "psi_range", psi_range.size(), "psi_v");        
    store_vector(file, p_prime, "p_prime", p_prime.size(), "psi_v");        
    store_vector(file, p, "p", p.size(), "psi_v");        
    store_vector(file, ff_prime, "ff_prime", ff_prime.size(), "psi_v");        
    store_vector(file, f_prime, "f_prime", f_prime.size(), "psi_v");        
    store_vector(file, f, "f", f.size(), "psi_v");        
}

void calculate_outputs(Field<Vector> &J, Field<Vector> &B, Field<double> &p, Field<double> &f, Field<double> &psi,
				  Field<double> &F, const InitialCondition &cond, const Parameters &param) {
  const Grid &grid = J.grid;
  Field<double> psi_N = get_normalized(psi, param);
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      if (grid.boundary[i][j].domain != Domain::INT)
        continue;
      J[i, j].phi = -F[i, j] / (grid.r[i] * mu0);
      B[i, j].phi = cond.f(psi_N[i, j]) / grid.r[i];
      double phi_r = (psi[i + 1, j] - psi[i - 1, j]) / (2 * grid.h);
      double phi_z = (psi[i, j + 1] - psi[i, j - 1]) / (2 * grid.h);
      B[i, j].r = -phi_z / grid.r[i];
      B[i, j].z = phi_r / grid.r[i];
      J[i, j].r = -phi_z * cond.f_prime(psi_N[i, j]) / (mu0 * grid.r[i]);
      J[i, j].z = phi_r * cond.f_prime(psi_N[i, j]) / (mu0 * grid.r[i]);
	  p[i ,j] = cond.p(psi_N[i, j]);
	  f[i ,j] = cond.f(psi_N[i, j]);
    }
  }
}

int main() {
  std::cout << std::fixed << std::setprecision(15);
  Parameters param;
  initialize(param);

  if (param.select == "solovev") {
	run_solovev(param);
	return EXIT_SUCCESS;
  } else if (param.select == "polynomial") {
	run_polynomial(param);
	return EXIT_SUCCESS;
  }

  double sigma = param.sigma0;
  auto cond = get_initial_condition(param);
  Grid grid(param);
  netCDF::NcFile file(param.ofname, netCDF::NcFile::replace);

  Field<double> psi(grid, param, 1);
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      if (grid.boundary[i][j].domain != Domain::INT)
        continue;
      psi[i, j] = (grid.r[i] - param.R) * (grid.r[i] - param.R) +
                  grid.z[j] * grid.z[j] + 0.9;
    }
  }
  Field<double> F(grid, param, 0);
  // outputs unnormalized psi
  std::pair<int, int> min_index =
      solve_default_unnormalized_output(param, grid, psi, F, sigma, *cond);
  // outputs unnormalized psi and F
  solve_unnormalized(param, grid, psi, F, sigma, *cond);

  Field<Vector> J(grid, param, {0, 0, 0});
  Field<Vector> B(grid, param, {0, 0, 0});
  Field<double> p(grid, param, 0);
  Field<double> f(grid, param, 0);
  calculate_outputs(J, B, p, f, psi, F, *cond, param);

  auto [i_min, j_min, i_max, _, psi_min, psi_max] = get_range(psi, param);  
  double psi_bdry = psi_max - psi_min;
  sigma = 1 / (psi_bdry * psi_bdry);
  double I = get_Ip(F, sigma);
  std::cout << "Sigma: " << sigma << "\n";
  std::cout << "I: " << I << "\n";

  int N_range = i_max - i_min + 1;
  std::vector<double> psi_range(N_range, 0);
  std::vector<double> psi_range_N(N_range, 0);
  std::vector<double> f_range(N_range, 0);
  std::vector<double> p_range(N_range, 0);
  for (int i = 0; i < N_range; i++) {
    int ri = i_min + i;
    int zi = j_min;
    psi_range[i] = psi[ri, zi];
    psi_range_N[i] = (psi_range[i] - psi_min) / (psi_max - psi_min);
    p_range[i] = (*cond).p(psi_range_N[i]);
    f_range[i] = (*cond).f(psi_range_N[i]);
  }
  store_vector(file, psi_range, "psi_range", psi_range.size(), "psi_range");
  store_vector(file, p_range, "p_range", p_range.size(), "psi_range");
  store_vector(file, f_range, "f_range", f_range.size(), "psi_range");

  store_vector(file, grid.r, "r", grid.N_r, "r");
  store_vector(file, grid.z, "z", grid.N_z, "z");
  store_field(file, psi, "psi", {"r", "z"});
  store_field(file, F, "F", {"r", "z"});
  store_field(file, p, "p", {"r", "z"});
  store_field(file, f, "f", {"r", "z"});
  store_vectorfield(file, J, "J", {"r", "z"});
  store_vectorfield(file, B, "B", {"r", "z"});
  std::cout << "\nWrote the result to: " << param.ofname << std::endl;
  return EXIT_SUCCESS;
}
