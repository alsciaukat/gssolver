#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <ranges>
#include <string>
#include <chrono>
#include <toml++/toml.hpp>

#include "base.hpp"
#include "dataio.hpp"

void update_F(Field<double> &F, Field<double> &psi, InitialCondition &cond, const Parameters &param) {
  const Grid &grid = psi.grid;
  const std::vector<double> &r = grid.r;
  for (int i = 0; i < grid.N_r; i++) {
    for (int j = 0; j < grid.N_z; j++) {
      F[i, j] = -(mu0 * r[i] * r[i] * cond.p_prime(psi[i, j]) +
                  cond.gg_prime(psi[i, j])) + 1e-100;
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

void normalize(Field<double> &psi, const Parameters &param) {
  const Grid &grid = psi.grid;
  // normalize psi between 0 at the minimum and 1 at the boundary
  double psi_min = INFINITY;
  int i_min;
  int j_min;
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

void sor_old(Field<double> &psi, Field<double> &F,
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
          // psi_delta =
		  // 	((1 - h / (2*r[i])) * psi[i + 1, j] + (1 + h / (2*r[i])) * psi[i - 1, j] -
          //      4 * psi[i, j] + psi[i, j + 1] + psi[i, j - 1] -
          //      h * h * sigma * F[i, j]) /
          //     4;
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
	  error_sq_sum += std::abs((a[i, j] - b[i, j]) / a[i, j]);
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
  SET(output, print, bool, print, false)
  SET(output, N_gr, int, N_gr, 17)
  SET(output, logh_low, double, logh_low, -2.4)
  SET(output, delta_logh, double, delta_logh, 0.1)
  SET(initial_condition, type, std::string, ictype, "solovev");
  SET(initial_condition, beta0, double, beta0, 0.5);
  SET(initial_condition, m, double, m, 2);
  SET(initial_condition, n, double, n, 1);
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
  if (param.ictype == "polynomial") return std::make_unique<PolynomialCondition>(param);
  return std::make_unique<SolovevCondition>(param);
}

void solve_default(const Parameters &param, const Grid &grid, Field<double> &psi, Field<double> &F, double &sigma, InitialCondition &cond, std::vector<double> &max_errors, std::vector<double> &l2_errors) {
  if (max_errors.size() > 0 || l2_errors.size() > 0) throw std::invalid_argument("The output arguments `max_errors` and `l2_errors` must be empty");
  // initializations
  update_F(F, psi, cond, param);
  update_sigma(sigma, F, param);

  int n_fix = 0;
  Field<double> psi_p(grid, param, 1);
  double max_error;
  double l2_error;

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
    std::cout << "Sigma: " << sigma << "\n";

    sor(psi, F, sigma, param);
	if (param.ictype != "solovev") normalize(psi, param);

    update_F(F, psi, cond, param);
    update_sigma(sigma, F, param);
	max_error = get_max_error(psi, psi_p);
	l2_error = get_l2_error(psi, psi_p);
	max_errors.push_back(max_error);
	l2_errors.push_back(l2_error);

    std::cout << "Error: " << max_error << "\n";        
  } while (max_error > param.e_fix);
}

void store_default(netCDF::NcFile &file, const Parameters &param, Field<double> &psi, Field<double> &F, double &sigma) {
  // output to a file 
  auto offormat_v = param.offormat | std::views::transform([](auto c){ return static_cast<char>(std::tolower(c)); });
  std::string offormat_l(offormat_v.begin(), offormat_v.end());

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
        error_sq_cum += error_tmp * error_tmp;
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

  solve_default(param, grid, psi, F, sigma, *cond, max_errors, l2_errors);

  store_field(file, psi, "psi", {"r", "z"});
  store_field(file, F, "F", {"r", "z"});
  store_vector(file, grid.r, "r", grid.N_r, "r");
  store_vector(file, grid.z, "z", grid.N_z, "z");
  store_vector(file, max_errors, "max_error", max_errors.size(), "i");
  store_vector(file, l2_errors, "l2_error", l2_errors.size(), "i");
}

int main() {
  std::cout << std::fixed << std::setprecision(6);
  Parameters param;
  initialize(param);
  if (param.ictype == "solovev") run_solovev(param); 
  else if (param.ictype == "polynomial") run_polynomial(param); 
  return EXIT_SUCCESS;
}
