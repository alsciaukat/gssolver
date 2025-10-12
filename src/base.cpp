#include "base.hpp"
#include <cmath>
#include <stdexcept>
#include <numbers>

const double mu0 = 4 * std::numbers::pi * 1e-7;

Vector::Vector(double r, double phi, double z) : r(r), phi(phi), z(z) {}

Vector Vector::operator+(const Vector &v) const {
  return Vector(r + v.r, phi + v.phi, z + v.z);
}

Vector Vector::operator-(const Vector &v) const {
  return Vector(r - v.r, phi - v.phi, z - v.z);
}  

Vector Vector::operator*(const double a) const {
  return Vector(a * r, a * phi, a * z);
}

Vector operator*(double a, const Vector &v) {
  return Vector(a * v.r, a * v.phi, a * v.z);
}  

Vector Vector::operator/(const double a) const {
  return Vector(r / a, phi / a, z / a);
}

Grid::Grid(const Parameters &param)
	: R(param.R), a(param.a), b(param.b), c0(param.c0), k(param.k),
	  tolerance(param.tolerance) {
  /*
    Solov'ev's solution is used to set the boundary.
    The p(psi) and g(psi) might not match the
    Solov'ev's, but they are used anyway assuming that it is not going to be far
    from the reality.

    The boundary is given by the flux surface given by:
      ψ(r, z) = k
    where
      ψ(r, z) := 0.5 (b + c₀) R² z² + c₀ R ζ(r) z²+ 0.5 (a - c₀) R² ζ(r)²
	  ζ(r) := (r² - R²)/(2R)
    with default parameters based on the dimension of KSTAR:
      R = 1.8, a = 0.5, b = 0.5, c₀ = 0.1, k = 0.11

	See /Vaccum solution for Solov'ev's equilibraum configuration in tokamaks/
	by Tao Xu and Richard Fitzpatrick.
  */

  if (c0 > a || -b > c0)
    throw std::invalid_argument("Invalid parameters: a, b, c0");
  if (param.N < 2)
    throw std::invalid_argument("Size too small: N");

  // analytically find the bounding box of the boundary.
  r_max = std::sqrt(R * R + std::sqrt(8 * k / (a - c0)));
  r_min = std::sqrt(R * R - std::sqrt(8 * k / (a - c0)));
  double p = c0 * c0 / (2 * (a - c0));
  double q = -(b + c0) * R * R / 2;
  z_max = std::sqrt((-q - std::sqrt(q * q - 4 * p * k)) / (2 * p));
  z_min = -z_max;
  if (std::isnan(r_min) || std::isnan(r_max) || std::isnan(z_max))
	throw std::invalid_argument(
		"The parameters are unphysical. Try setting k lower.");

  // automatically calculate the number of points in z direction
  // to have the same spacing as in r.
  h = (r_max - r_min + 2 * tolerance) / (param.N - 1);
  N_r = param.N;
  N_z = static_cast<int>(ceil((z_max - z_min + 2 * tolerance) / h)) + 1;

  boundary.resize(N_r);
  for (auto &c : boundary)
    c.resize(N_z);
  r.resize(N_r);
  z.resize(N_z);

  for (int i = 0; i < N_r; i++)
    r[i] = r_min - tolerance + i * h;
  for (int j = 0; j < N_z; j++)
    z[j] = z_min - tolerance + j * h;

  // analytically find the inhomogeneity at the boundary.
  // `alpha` and `beta` represent the relative distance
  // from the last interior point to the boundary
  // compared to the regular spacing `h`.
  // 
  // If boundary point has less r coordinate then the last adjacent interior point (LAIP),
  // `alpha1` is set at the LAIP in the `BoundaryInfo`.
  // Similarly, `alpha2` if greater r coordinate, `beta_1` if less r coordinate,
  // `beta2` if greater r coordinate.
  for (int j = 0; j < N_z; j++) {
    if (z[j] < z_min || z[j] > z_max)
      continue;
    double r_l = std::sqrt(R * R - 2 * c0 * z[j] * z[j] / (a - c0) -
                      2 *
                          std::sqrt(c0 * c0 * std::pow(z[j], 4) + 2 * (a - c0) * k -
                               (a - c0) * (b + c0) * R * R * z[j] * z[j]) /
                          (a - c0));
    int i_l = static_cast<int>(ceil((r_l - r[0]) / h));
    double alpha_1l = (r[i_l] - r_l) / h;
    boundary[i_l][j] = {Domain::BDRY, alpha_1l, 1, 1, 1};

    double r_r = std::sqrt(R * R - 2 * c0 * z[j] * z[j] / (a - c0) +
                      2 *
                          std::sqrt(c0 * c0 * std::pow(z[j], 4) + 2 * (a - c0) * k -
                               (a - c0) * (b + c0) * R * R * z[j] * z[j]) /
                          (a - c0));
    int i_r = static_cast<int>(floor((r_r - r[0]) / h));
    double alpha_2r = (r_r - r[i_r]) / h;
    if (boundary[i_r][j].domain == Domain::BDRY) {
      boundary[i_r][j].alpha2 = alpha_2r;
    } else {
      boundary[i_r][j] = {Domain::BDRY, 1, 1, alpha_2r, 1};
    }
    for (int i = i_l + 1; i < i_r; i++)
      boundary[i][j] = {Domain::INT, 1, 1, 1, 1};
  }

  for (int i = 0; i < N_r; i++) {
    if (r[i] < r_min || r[i] > r_max)
      continue;
    double z_b =
        -std::sqrt((2 * k - (a - c0) * std::pow(r[i] * r[i] - R * R, 2) / 4) /
              ((b + c0) * R * R + c0 * (r[i] * r[i] - R * R)));
    int j_b = static_cast<int>(ceil((z_b - z[0]) / h));
    double beta_1b = (z[j_b] - z_b) / h;
    boundary[i][j_b].domain = Domain::BDRY;
    boundary[i][j_b].beta1 = beta_1b;

    double z_t = -z_b;
    int j_t = static_cast<int>(floor((z_t - z[0]) / h));
    double beta_2t = (z_t - z[j_t]) / h;
    boundary[i][j_t].domain = Domain::BDRY;
    boundary[i][j_t].beta2 = beta_2t;
  }
}

double Grid::solovev(double r, double z) const {
  return 0.5 * (b + c0) * R * R * z * z + 0.5 * c0 * (r * r - R * R) * z * z +
         (a - c0) * std::pow(r * r - R * R, 2) / 8;
}

template class Field<double>;


// interpolation is done using the second order Newton polynomials
template <typename T>
T Field<T>::interpolate_z(int i, int j, double zz) const {
  const std::vector<double> &z = grid.z;
  const std::vector<double> &r = grid.r;
  const double h = grid.h;
  return (zz - z[j - 1]) * (zz - z[j]) *
	(value[i][j + 1] - 2 * value[i][j] + value[i][j - 1]) / (2 * h * h) +
	(zz - z[j - 1]) * (value[i][j] - value[i][j - 1]) / h + value[i][j - 1];
}

template <typename T>
T Field<T>::interpolate_r(int i, int j, double rr) const {
  const std::vector<double> &z = grid.z;
  const std::vector<double> &r = grid.r;
  const double h = grid.h;
  return (rr - r[i - 1]) * (rr - r[i]) *
	(value[i + 1][j] - 2 * value[i][j] + value[i - 1][j]) / (2 * h * h) +
	(rr - r[i - 1]) * (value[i][j] - value[i - 1][j]) / h + value[i - 1][j];
}

template <typename T>
Field<T>::Field(const Grid &grid, const Parameters &param, T initial_value) : value(grid.N_r, std::vector<T>(grid.N_z, initial_value)), grid(grid) {};

template <typename T>
T &Field<T>::operator[](std::size_t i, std::size_t j) {
  return value[i][j];
}

template <typename T>
Field<T> &Field<T>::operator=(const Field<T> &other) {
  if (this != &other) {
	value = other.value;
  }    
  return *this;
}

// abstract class for initial conditions
InitialCondition::InitialCondition(const Parameters &param) : param(param) {};


double SolovevCondition::p_prime(double psi, double r) {
  return - param.a / mu0;
}

double SolovevCondition::gg_prime(double psi, double r) {
  return - param.b * param.R * param.R;
}

double PolynomialCondition::p_prime(double psi, double r) {
  return param.beta0 * std::pow(1 - std::pow(psi, param.m), param.n) / param.R;
}

double PolynomialCondition::gg_prime(double psi, double r) {
  return (1 - param.beta0) * mu0 * param.R * std::pow(1 - std::pow(psi, param.m), param.n);
}  

// H-mode p_prime, gg_prime
// You can see plot of this profile and adjust parameters:
// https://www.desmos.com/calculator/65uenvvx0u

const double d = 0.02;
const double betaH = 0.6;
const double p_prime0 = 0.18;
const double psi0 = 0.95;
const double h = 0.67;
const double D = 0.036;
const double j_bs0 = 0.4;

double HModeCondition::p_prime(double psi, double r) {
  return - ( betaH * std::pow(1 - std::pow(psi, param.m), param.n) + p_prime0) * 0.5*h * (std::tanh( (psi0 - psi) / d ) + 1) / r;
}  

double HModeCondition::gg_prime(double psi, double r) {
  return mu0 * r * ( std::pow(1 - std::pow(psi, param.m), param.n) - ( betaH * std::pow(1 - std::pow(psi, param.m), param.n) + p_prime0) * 0.5*h * (std::tanh( (psi0 - psi) / d ) + 1) );
}  
