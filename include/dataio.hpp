#ifndef DATAIO_H
#define DATAIO_H

#include <cstddef>
#include <string>
#include "base.hpp"

#ifdef USE_NETCDF

#include <netcdf>
int store_netcdf(netCDF::NcFile &file, Field<double>& psi, Field<double> &F, const double &sigma, const Parameters &param);

int store_field(netCDF::NcFile &file, Field<double> &field,
				const std::string vname,
				const std::vector<std::string> axnames);

int store_vectorfield(netCDF::NcFile &file, Field<Vector> &field,
					  const std::string vname, const std::vector<std::string> axnames);

int store_vector(netCDF::NcFile &file, const std::vector<double> &vector,
				 const std::string vname, const size_t size,
				 const std::string axname);

#endif

int store_csv(const Field<double>& field, const std::string& fname);

#endif
