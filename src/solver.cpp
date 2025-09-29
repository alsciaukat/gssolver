#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <ranges>
#include <string>
#include <toml++/toml.hpp>

#include "base.hpp"
#include "dataio.hpp"

const double mu0 = 4 * std::numbers::pi * 1e-7;

void update_F(Field<double> &F, Field<double> &psi, InitialCondition &cond, const Parameters &param) {
  const Grid &grid = psi.grid;
  const std::vector<double> &r = grid.r;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      F[i, j] = -(mu0 * r[i] * r[i] * cond.p_prime(psi[i, j]) +
                  cond.gg_prime(psi[i, j]));
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

  sigma = param.I_p / (param.psi_bdry * I);
}

void normalize(Field<double> &psi, const Parameters &param) {
  const Grid &grid = psi.grid;
  // normalize psi between 0 at the minimum and 1 at the boundary
  double psi_min = INFINITY;
  int i_min;
  int j_min;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      if (psi[i, j] < psi_min) {
        psi_min = psi[i, j];
        i_min = i;
        j_min = j;
      }
    }
  }
  int i = i_min;
  int j = j_min;
  if (grid.boundary[i][j].domain != Domain::INT)
    throw std::range_error("Minimum of psi is not interior");
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
  std::cout << "(r_min, z_min): (" << r_min << ", " << z_min << ")\n";

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

  std::cout << "psi(r_min, z_min) = " << psi_min << "\n";

  // normalize ψ
  for (double &x : psi.value | std::views::join) {
    x = (x - psi_min) / (param.psi_l - psi_min) * param.psi_l;
  }
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
          psi_delta =
              ((1 - h / r[i]) * psi[i + 1, j] + (1 + h / r[i]) * psi[i - 1, j] -
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

bool is_converged(Field<double> &psi_p, Field<double> &psi, const Parameters &param) {
  // check the convergence of fixed-point iteration
  bool converged = true;
  double N_r = psi.value.size();
  double N_z = psi.value[0].size();
  for (int i = 0; i < N_r; i++) {
    for (int j = 0; j < N_z; j++) {
      if (std::abs((psi[i, j] - psi_p[i, j]) / psi[i, j]) > param.e_fix)
        converged = false;
    }
  }
  return converged;
}

void print_array(const Array<double> &a) {
  for (auto &b : a) {
    for (auto &c : b) {
      std::cout << c << "\t";
    }
    std::cout << "\n";
  }
}

#define SET(SEC, VAR, TYPE, IVAR, DEFAULT) param.IVAR = config[#SEC][#VAR].value<TYPE>().value_or(DEFAULT);

int initialize(Parameters &param) {
  toml::table config = toml::parse_file("config.toml");
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
  SET(initial_condition, type, std::string, ictype, "solovev");
  SET(solver, e_fix, double, e_fix, 0.01);
  SET(solver, e_sor, double, e_sor, 0.01);
  SET(solver, e_min, double, e_min, 0.01);
  SET(solver, omega, double, omega, 1.9);
  SET(solver, I_p, double, I_p, 500e3);
  SET(solver, psi_l, double, psi_l, 1);
  SET(solver, psi_bdry, double, psi_bdry, 0.1);
  return 0;
}

std::unique_ptr<InitialCondition>
get_initial_condition(const Parameters &param) {
  // *factory pattern*
  // this function returns a pointer to
  // a concrete class of `InitialCondition`.
  return std::make_unique<SolovevCondition>(param);
}

int main() {
  // initializations
  std::cout << std::fixed << std::setprecision(2);
  Parameters param;
  initialize(param);
  auto cond = get_initial_condition(param);
  Grid grid(param);
  Field<double> psi(grid, param, 1);
  Field<double> psi_p(grid, param, 1);
  Field<double> F(grid, param, 0);
  double sigma;

  int n_fix = 0;
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
    std::cout << "\nFixed-point Iteration: " << n_fix << "\n";

    psi_p = psi;
    update_F(F, psi, *cond, param);
    update_sigma(sigma, F, param);
    std::cout.unsetf(std::ios::fixed);
    std::cout << "Sigma: " << sigma << "\n" << std::fixed;

    sor(psi, F, sigma, param);
    normalize(psi, param);
  } while (!is_converged(psi_p, psi, param));

  // output to a file 
  auto offormat_v = param.offormat | std::views::transform([](auto c){ return static_cast<char>(std::tolower(c)); });
  std::string offormat_l(offormat_v.begin(), offormat_v.end());

  if (offormat_l == "netcdf") {
#ifdef USE_NETCDF
    store_netcdf(psi.value, param.ofname, "psi");
#else
    throw std::invalid_argumet("Output to NetCDF is not supported. Change output.format to csv");
#endif
  } else if (offormat_l == "csv") {
	store_csv(psi.value, param.ofname);
  } else {
	throw std::invalid_argument("No output is written. Current output.format is not supported.");
  }    
  print_array(psi.value);

  return EXIT_SUCCESS;
}
