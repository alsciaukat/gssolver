#ifndef DATAIO_H
#define DATAIO_H

#include <string>
#include "base.hpp"

#ifdef USE_NETCDF

#include <netcdf>
int store_netcdf(const Field<double>& psi, const Field<double> &F, const std::string& fname);

#endif

int store_csv(const Field<double>& field, const std::string& fname);

#endif
