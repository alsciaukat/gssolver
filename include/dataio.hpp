#ifndef DATAIO_H
#define DATAIO_H

#include <string>
#include "base.hpp"

#ifdef USE_NETCDF

#include <netcdf>
int store_netcdf(const Array<double>& arr, const std::string& fname, const std::string& vname);

#endif

int store_csv(const Array<double>& arr, const std::string& fname);

#endif
