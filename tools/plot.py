import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

file = nc.Dataset("result.nc", "r")

psi = file.variables["psi"][:]
r = file.variables["r"][:]
z = file.variables["z"][:]

R, Z = np.meshgrid(r, z, indexing="ij")

plt.figure(figsize=(6, 5))
contours = plt.contour(R, Z, psi, levels=20, cmap='viridis')

plt.xlabel("r")
plt.ylabel("z")
plt.title("Contours of ψ(r,z)")
plt.colorbar(contours, label="ψ")
plt.axis("equal")
plt.show()
