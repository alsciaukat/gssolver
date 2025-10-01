#include "dataio.hpp"
#include "base.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

using namespace netCDF;

int store_netcdf(Field<double> &psi, Field<double> &F, const double &sigma,
				 const InitialCondition &cond, const Parameters &param) {
  const Grid &grid = psi.grid;
  const std::vector<double> &r = grid.r;
  const std::vector<double> &z = grid.z;
  size_t N_r = r.size();
  size_t N_z = z.size();
  std::vector<double> J_phi(N_r * N_z);
  std::vector<double> B_r(N_r * (N_z - 1));
  std::vector<double> B_z((N_r - 1) * N_z);

  for (int i = 0; i < N_r; i++) {
	for (int j = 0; j < N_z; j++) {
	  // Do I not need a minus sign?
	  J_phi[N_z*i + j] = sigma * param.psi_bdry * F[i, j] / (mu0 * r[i]);
	  if (j < N_z - 1)
		B_r[(N_z - 1) * i + j] = - sigma * param.psi_bdry * (psi[i, j + 1] - psi[i, j]) / (grid.h * r[i]);
	  if (i < N_r - 1)
		B_z[N_z * i + j] = sigma * param.psi_bdry * (psi[i + 1, j] - psi[i, j]) / (grid.h * r[i]);
	}
  }    

  // flatten 2D array into 1D
  std::vector<double> psi_f;
  std::vector<double> F_f;
  for (double x : psi.value | std::views::join)
	psi_f.push_back(x);
  for (double x : F.value | std::views::join)
	F_f.push_back(x);
  NcFile file(param.ofname, NcFile::replace);
  NcDim rDim = file.addDim("r", N_r);
  NcDim rpDim = file.addDim("rp", N_r-1);
  NcDim zDim = file.addDim("z", N_z);
  NcDim zpDim = file.addDim("zp", N_z-1);
  std::vector<NcDim> dims = {rDim, zDim};
  std::vector<NcDim> rpdims = {rpDim, zDim};
  std::vector<NcDim> zpdims = {rDim, zpDim};

  NcVar r_v = file.addVar("r", ncDouble, rDim);
  NcVar z_v = file.addVar("z", ncDouble, zDim);
  NcVar psi_v = file.addVar("psi", ncDouble, dims);
  NcVar F_v = file.addVar("F", ncDouble, dims);
  NcVar J_phi_v = file.addVar("J_phi", ncDouble, dims);
  NcVar B_r_v = file.addVar("B_r", ncDouble, zpdims);
  NcVar B_z_v = file.addVar("B_z", ncDouble, rpdims);
  r_v.putVar(psi.grid.r.data());
  z_v.putVar(psi.grid.z.data());
  psi_v.putVar(psi_f.data());
  F_v.putVar(F_f.data());
  J_phi_v.putVar(J_phi.data());
  B_r_v.putVar(B_r.data());
  B_z_v.putVar(B_z.data());
  return EXIT_SUCCESS;
}

int store_csv(const Field<double> &field, const std::string &fname) {
  std::ofstream file(fname);
  if (!file) return EXIT_FAILURE;
  for (auto &v : field.value) {
	size_t n = v.size();
	for (int i = 0; i < n - 1; i++) {
	  file << v[i] << ", ";
	}
	file << v[n-1] << "\n";
  }
  file.close();
  return EXIT_SUCCESS;
}  
