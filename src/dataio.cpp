#include "dataio.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

using namespace netCDF;

int store_netcdf(const Field<double> &field, const std::string &fname,
				 const std::string &vname) {
  std::vector<double> flat_arr;
  for (double x : field.value | std::views::join)
	flat_arr.push_back(x);
  size_t N_r = field.value.size();
  size_t N_z = field.value[0].size();
  NcFile file(fname, NcFile::replace);
  NcDim rDim = file.addDim("r", N_r);
  NcDim zDim = file.addDim("z", N_z);
  std::vector<NcDim> dims = {rDim, zDim};

  NcVar r = file.addVar("r", ncDouble, rDim);
  NcVar z = file.addVar("z", ncDouble, zDim);
  NcVar var = file.addVar(vname, ncDouble, dims);
  r.putVar(field.grid.r.data());
  z.putVar(field.grid.z.data());
  var.putVar(flat_arr.data());
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
