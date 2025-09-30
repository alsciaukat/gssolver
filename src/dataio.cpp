#include "dataio.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

using namespace netCDF;

int store_netcdf(const Field<double> &psi, const Field<double> &F, const std::string &fname) {
  std::vector<double> psi_f;
  std::vector<double> F_f;
  for (double x : psi.value | std::views::join)
	psi_f.push_back(x);
  for (double x : F.value | std::views::join)
	F_f.push_back(x);
  size_t N_r = psi.value.size();
  size_t N_z = psi.value[0].size();
  NcFile file(fname, NcFile::replace);
  NcDim rDim = file.addDim("r", N_r);
  NcDim zDim = file.addDim("z", N_z);
  std::vector<NcDim> dims = {rDim, zDim};

  NcVar r = file.addVar("r", ncDouble, rDim);
  NcVar z = file.addVar("z", ncDouble, zDim);
  NcVar psi_v = file.addVar("psi", ncDouble, dims);
  NcVar F_v = file.addVar("F", ncDouble, dims);
  r.putVar(psi.grid.r.data());
  z.putVar(psi.grid.z.data());
  psi_v.putVar(psi_f.data());
  F_v.putVar(F_f.data());
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
