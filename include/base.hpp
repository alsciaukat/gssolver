#ifndef BASE_H
#define BASE_H

#include <cstddef>
#include <string>
#include <vector>

extern const double mu0;

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
  std::string select;
  std::string gtype;
  int N;
  double R, a, b, c0, k, tolerance;
  double e_fix, e_sor, e_min, omega;
  double I_p, psi_bdry;
  double psi_l;
  std::string ofname;
  std::string offormat;
  int N_gr;
  double logh_low;
  double delta_logh;
  std::string ictype;
  bool print;
  double beta0, m, n;
  double A;
};

class Vector {
public:
  double r;
  double phi;
  double z;

  Vector(double r, double phi, double z);

  Vector operator+(const Vector &v) const;
  Vector operator-(const Vector &v) const;
  Vector operator*(const double a) const;
  Vector operator/(const double a) const;
  friend Vector operator*(double a, const Vector &v);
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

  Field(const Grid &grid, const Parameters &param, T initial_value = 0);

  T &operator[](std::size_t i, std::size_t j);
  Field<T> &operator=(const Field<T> &other);
  T interpolate_z(int i, int j, double zz) const;
  T interpolate_r(int i, int j, double rr) const;
};

template class Field<Vector>;

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

class PolynomialCondition : public InitialCondition {
  using InitialCondition::InitialCondition;
  double p_prime(double psi) override;  
  double gg_prime(double psi) override;  
};


#endif
