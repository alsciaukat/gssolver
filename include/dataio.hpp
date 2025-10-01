#ifndef DATAIO_H
#define DATAIO_H

#include <string>
#include "base.hpp"

#ifdef USE_NETCDF

#include <netcdf>
int store_netcdf(Field<double>& psi, Field<double> &F, const double &sigma, const InitialCondition &cond, const Parameters &param);

#endif

int store_csv(const Field<double>& field, const std::string& fname);

#endif
