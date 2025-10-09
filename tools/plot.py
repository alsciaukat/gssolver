#!/usr/bin/python

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys


def compare_solovev(numerical, theory, R, Z):
    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    overlay = axes.contourf(R, Z, theory, levels=np.linspace(0, 0.11, 11),
                            cmap='viridis', alpha=0.5)
    contours = axes.contour(R, Z, numerical, levels=np.linspace(0, 0.11, 11),
                            cmap='plasma')
    axes.set(xlabel="r", ylabel="z")
    axes.axis("equal")

    plt.clabel(contours, inline=True)
    plt.colorbar(overlay)
    plt.show()


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


def plot_solovev():
    file = nc.Dataset("solovev.nc", "r")
    h = file.variables["h"][:]
    N = file.variables["N"][:]
    max_error = file.variables["max_error"][:]
    l2_error = file.variables["l2_error"][:]
    time2run = file.variables["time2run"][:]

    n = int(sys.argv[1])
    psi = file.variables[f"psi{n}"][:]
    psi_t = file.variables[f"psi_t{n}"][:]
    r = file.variables[f"r{n}"][:]
    z = file.variables[f"z{n}"][:]
    R, Z = np.meshgrid(r, z, indexing="ij")
    compare(psi, psi_t, R, Z)

    # logh = np.log10(h)
    # loge = np.log10(l2_error)

    # coeffs = np.polyfit(logh, loge, 1)
    # b = coeffs[0]
    # a = 10 ** coeffs[1]
    # e_fit = a * h ** 3

    # plt.loglog(h, l2_error, "o")
    # plt.loglog(h, e_fit, "-")

    # print(f"Fitted Exponent: {b}")
    # plt.plot(N, time2run)
    # plt.show()


def plot_iter_error():
    file = nc.Dataset("test.nc", "r")
    max_errors = file.variables["max_error"][:]
    l2_errors = file.variables["l2_error"][:]
    plt.plot(l2_errors)
    plt.yscale('log')
    plt.show()


def debug():
    file = nc.Dataset("test.nc", "r")
    file2 = nc.Dataset("solovev.nc", "r")
    n = int(sys.argv[1])
    psi = file.variables[f"psi{n}"][:]
    psi_t = file2.variables["psi_t1"][:]
    r = file.variables[f"r{n}"][:]
    z = file.variables[f"z{n}"][:]
    r2 = file2.variables["r1"][:]
    z2 = file2.variables["z1"][:]
    R, Z = np.meshgrid(r, z, indexing="ij")
    R2, Z2 = np.meshgrid(r2, z2, indexing="ij")
    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    overlay = axes.contourf(R2, Z2, psi_t, levels=np.linspace(0, 0.11, 11),
                            cmap='viridis', alpha=0.5)
    contours = axes.contour(R, Z, psi, levels=np.linspace(0, 1, 11),
                            cmap='plasma')
    axes.set(xlabel="r", ylabel="z")
    axes.axis("equal")

    plt.clabel(contours, inline=True)
    plt.colorbar(overlay)
    plt.show()


if __name__ == "__main__":
    # compare_solovev(psi_solovev, psi_theory, R_solovev, Z_solovev)
    # compare(psi_poly, psi_theory_poly, R_poly, Z_poly)
    # show_field(R_solovev, Z_solovev, J_phi_solovev, "J_ϕ")
    # show_field(Rrp_solovev, Zrp_solovev, B_z_solovev, "B_z")
    # show_psi(R_poly, Z_poly, psi_poly, "psi")
    # show_field(Rrp_poly, Zrp_poly, B_z_poly, "B_z")
    # show_field(R_poly, Z_poly, J_phi_poly, "J_ϕ")
    # compare_solovev(psi, psi_t, R, Z)
    plot_iter_error()
    # plot_solovev()
    # debug()
