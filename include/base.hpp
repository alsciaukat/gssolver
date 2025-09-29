#ifndef BASE_H
#define BASE_H

#include <cstddef>
#include <string>
#include <vector>
#include <numbers>

template <typename T>
using Array = std::vector<std::vector<T>>;

enum class Domain {
  EXT,  // exterior
  BDRY, // boundary
  INT,  // interior
};

struct BoundaryInfo {
  Domain domain = Domain::EXT;
  double alpha1 = 1;
  double beta1 = 1;
  double alpha2 = 1;
  double beta2 = 1;
};

struct Parameters {
  int N;
  double R, a, b, c0, k, tolerance;
  double e_fix = 0.01, e_sor = 0.01, e_min = 0.01, omega = 1.9; // relaxation coefficient. 2 diverges.
  double A = 1;
  double mu0 = 4 * std::numbers::pi * 1e-7;
  double sigma = 0.5;
  double Ip = 500e3, psi_bdry = 0.1;
  double psi_l = 1;
  std::string ofname;
  std::string offormat;
  std::string ictype;
};

class Grid {
public:
  double R, a, b, c0, k;
  std::vector<double> r;
  std::vector<double> z;
  double r_max, r_min, z_max, z_min, h;
  int N_r, N_z;
  double tolerance;
  Array<BoundaryInfo> boundary;

  Grid(const Parameters &param);

  double solovev(double r, double z) const;
};

template<typename T>
class Field {
public:
  Array<T> value;
  const Grid &grid;

  Field(const Grid &grid, const Parameters &param, T initial_value);

  T &operator[](std::size_t i, std::size_t j);
  Field<T> &operator=(const Field<T> &other);
  double interpolate_z(int i, int j, double zz) const;
  double interpolate_r(int i, int j, double rr) const;
};

class InitialCondition {
public:
  const Parameters &param;
  InitialCondition(const Parameters &param);
  virtual double p_prime(double psi) = 0;
  virtual double gg_prime(double psi) = 0;
  virtual ~InitialCondition() = default;
};

class SolovevCondition : public InitialCondition {
public:
  using InitialCondition::InitialCondition;
  double p_prime(double psi) override;  
  double gg_prime(double psi) override;  
};

class QuadraticCondition : public InitialCondition {
  using InitialCondition::InitialCondition;
  double p_prime(double psi) override;  
  double gg_prime(double psi) override;  
};


#endif
