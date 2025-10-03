import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

mu0 = 4 * np.pi * 1e-7

data_solovev = nc.Dataset("solovev.nc", "r")
data_poly = nc.Dataset("polynomial.nc", "r")
data_test = nc.Dataset("test.nc", "r")

# 0.5 * (b + c0) * R * R * z * z + 0.5 * c0 * (r * r - R * R) * z * z
# + (a - c0) * std::pow(r * r - R * R, 2) / 8;

RR = 1.8
a = 0.5
b = 0.5
c0 = 0.1


def solovev(r, z):
    return (b + c0) * RR**2 * z**2 / 2 + c0 * (r**2 - RR**2) * z**2 / 2 + \
     (a - c0) * (r**2 - RR**2)**2 / 8


psi_solovev = data_solovev.variables["psi"][:]
F_solovev = data_solovev.variables["F"][:]
J_phi_solovev = data_solovev.variables["J_phi"][:]
B_r_solovev = data_solovev.variables["B_r"][:]
B_z_solovev = data_solovev.variables["B_z"][:]
r_solovev = data_solovev.variables["r"][:]
z_solovev = data_solovev.variables["z"][:]

psi_poly = data_poly.variables["psi"][:]
J_phi_poly = data_poly.variables["J_phi"][:]
B_r_poly = data_poly.variables["B_r"][:]
B_z_poly = data_poly.variables["B_z"][:]
F_poly = data_poly.variables["F"][:]
r_poly = data_poly.variables["r"][:]
z_poly = data_poly.variables["z"][:]

R_solovev, Z_solovev = np.meshgrid(r_solovev, z_solovev, indexing="ij")
Rrp_solovev, Zrp_solovev = np.meshgrid(r_solovev[:-1], z_solovev, indexing="ij")
Rzp_solovev, Zzp_solovev = np.meshgrid(r_solovev, z_solovev[:-1], indexing="ij")
psi_theory = solovev(R_solovev, Z_solovev)

R_poly, Z_poly = np.meshgrid(r_poly, z_poly, indexing="ij")
Rrp_poly, Zrp_poly = np.meshgrid(r_poly[:-1], z_poly, indexing="ij")
Rzp_poly, Zzp_poly = np.meshgrid(r_poly, z_poly[:-1], indexing="ij")

psi_theory_poly = solovev(R_poly, Z_poly)


def compare(numerical, theory, R, Z):
    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    overlay = axes.contourf(R, Z, theory, levels=np.linspace(0, 0.11, 11),
                            cmap='viridis', alpha=0.5)
    contours = axes.contour(R, Z, numerical, levels=np.linspace(0, 1, 11),
                            cmap='plasma')
    axes.set(xlabel="r", ylabel="z")
    axes.axis("equal")

    plt.clabel(contours, inline=True)
    plt.colorbar(overlay)
    plt.show()


def show_psi(R, Z, field):
    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    contours = axes.contour(R, Z, field, levels=np.linspace(0, 1, 11),
                            cmap='plasma')
    axes.set(xlabel="r", ylabel="z", title="Contours of ψ(r,z)")
    axes.axis("equal")
    plt.colorbar(contours, label="ψ")
    plt.show()


def show_field(R, Z, field, title: str):
    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    contours = axes.contourf(R, Z, field, levels=40,
                             cmap='plasma')
    axes.set(xlabel="r", ylabel="z", title=title)
    axes.axis("equal")
    plt.colorbar(contours)
    plt.show()


if __name__ == "__main__":
    # compare(psi_solovev, psi_theory, R_solovev, Z_solovev)
    compare(psi_poly, psi_theory_poly, R_poly, Z_poly)
    # show_field(R_solovev, Z_solovev, J_phi_solovev, "J_ϕ")
    # show_field(Rrp_solovev, Zrp_solovev, B_z_solovev, "B_z")
    # show_psi(R_poly, Z_poly, psi_poly, "psi")
    # show_field(Rrp_poly, Zrp_poly, B_z_poly, "B_z")
    # show_field(R_poly, Z_poly, J_phi_poly, "J_ϕ")
