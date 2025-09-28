#include "dataio.hpp"
#include <cstdlib>
#include <ranges>
#include <vector>

using namespace netCDF;

int store_netcdf(const Array<double>& arr, const std::string& fname, const std::string& vname) {
  std::vector<double> flat_arr;
  for (double x : arr | std::views::join) flat_arr.push_back(x);
  size_t N_r = arr.size();
  size_t N_z = arr[0].size();
  NcFile file(fname, NcFile::replace);
  NcDim rDim = file.addDim("r", N_r);
  NcDim zDim = file.addDim("z", N_z);
  std::vector<NcDim> dims = {rDim, zDim};
  NcVar var = file.addVar(vname, ncDouble, dims);
  var.putVar(flat_arr.data());
  return EXIT_SUCCESS;
}

