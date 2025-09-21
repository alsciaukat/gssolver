#include <vector>

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

class Grid {
public:
  double R, a, b, c0, k;
  std::vector<double> r;
  std::vector<double> z;
  double r_max, r_min, z_max, z_min, h;
  int N_r, N_z;
  double tolerance;
  Array<BoundaryInfo> boundary;

  Grid(size_t N, double R = 3, double a = 0.6, double b = 0.2, double c0 = 0.1,
       double k = 3, double tolerance = 0.1);

  double interpolate_z(const Array<double> &psi, int i, int j, double zz) const;
  double interpolate_r(const Array<double> &psi, int i, int j, double rr) const;
  double solovev(double r, double z) const;
};
