#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import tomllib
from argparse import ArgumentParser

fig_w = 10
fig_h = 7


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


def show_field(R, Z, field, title: str):
    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    contours = axes.contourf(R, Z, field, levels=40,
                             cmap='plasma')
    axes.set(xlabel="r", ylabel="z", title=title)
    axes.axis("equal")
    plt.colorbar(contours)
    plt.show()


def plot_solovev():
    file = Dataset("solovev.nc", "r")
    h = file.variables["h"][:]
    N = file.variables["N"][:]
    max_error = file.variables["max_error"][:]
    l2_error = file.variables["l2_error"][:]
    time2run = file.variables["time2run"][:]

    logh = np.log10(h)
    loge = np.log10(l2_error)

    coeffs = np.polyfit(logh, loge, 1)
    b = coeffs[0]
    a = 10 ** coeffs[1]
    e_fit = a * h ** 3

    plt.loglog(h, l2_error, "o", label="Raw")
    plt.loglog(h, e_fit, "-", label="Fitted")
    plt.xlabel("h (Grid spacing)")
    plt.ylabel("L2 Error")
    plt.legend()

    print(f"Fitted Exponent: {b}")
    # plt.plot(N, time2run)
    plt.show()


def plot_iter_error():
    file = Dataset("test.nc", "r")
    max_errors = file.variables["max_error"][:]
    l2_errors = file.variables["l2_error"][:]
    plt.plot(l2_errors)
    plt.xlabel("# Iteration")
    plt.ylabel("L2 Error (from previous)")
    plt.yscale('log')
    plt.show()


def overlay_solovev(file: Dataset):
    file2 = Dataset("solovev.nc", "r")
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


def extract_test(fname: str):
    file = Dataset(fname, "r")
    prefix = "assets/" + fname + "_"
    r = file.variables["r"][:]
    z = file.variables["z"][:]
    R, Z = np.meshgrid(r, z, indexing="ij")
    stride = 7
    R_s, Z_s = np.meshgrid(r[::stride], z[::stride], indexing="ij")

    for name, colorlabel, title in [("J", r"$J$ [$\rm A/m^2$]", r"Poloidal $J$"), ("B", r"$B$ [T]", r"Poloidal $B$")]:
        v = file.variables[name][:]
        v_mag = np.sqrt(v[:, :, 0]**2 + v[:, :, 2]**2)
        v_mag_s = np.sqrt(v[::stride, ::stride, 0]**2
                          + v[::stride, ::stride, 2]**2)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        quivers = ax.quiver(R_s, Z_s,
                            v[::stride, ::stride, 0] /
                            v_mag_s, v[::stride, ::stride, 2] / v_mag_s,
                            v_mag_s, pivot='mid', width=0.006)
        fig.colorbar(quivers, ax=ax, label=colorlabel)
        ax.set(xlabel=r"$r$ [m]", ylabel=r"$z$ [m]", title=title, aspect="equal")
        fig.savefig(prefix + name + "_quiver.png")
        del fig, ax

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        streams = ax.streamplot(R.T, Z.T, v[:, :, 0].T, v[:, :, 2].T,
                                color=v_mag.T, density=0.9, arrowsize=3,
                                broken_streamlines=False, linewidth=2)
        fig.colorbar(plt.cm.ScalarMappable(
                norm=plt.Normalize(v_mag.min(), v_mag.max()),
                cmap=streams.lines.get_cmap()),
                     ax=ax, label=colorlabel)
        ax.set(xlabel=r"$r$ [m]", ylabel=r"$z$ [m]", title=title, aspect="equal")
        fig.savefig(prefix + name + ".png")
        del fig, ax


def extract_plots(fname: str, ictype: str):
    file = Dataset(fname, "r")
    prefix = "assets/" + fname + "_"
    r = file.variables["r"][:]
    z = file.variables["z"][:]
    psi_range = file.variables["psi_range"][:]
    R, Z = np.meshgrid(r, z, indexing="ij")
    stride = 7
    R_s, Z_s = np.meshgrid(r[::stride], z[::stride], indexing="ij")

    for name, colorlabel, title in [("J", r"$J$ [$\rm A/m^2$]", r"Poloidal $J$"), ("B", r"$B$ [T]", r"Poloidal $B$")]:
        v = file.variables[name][:]
        v_mag = np.sqrt(v[:, :, 0]**2 + v[:, :, 2]**2)
        v_mag_s = np.sqrt(v[::stride, ::stride, 0]**2
                          + v[::stride, ::stride, 2]**2)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        quivers = ax.quiver(R_s, Z_s,
                            v[::stride, ::stride, 0] /
                            v_mag_s, v[::stride, ::stride, 2] / v_mag_s,
                            v_mag_s, pivot='mid', width=0.006)
        fig.colorbar(quivers, ax=ax, label=colorlabel)
        ax.set(xlabel=r"$r$ [m]", ylabel=r"$z$ [m]", title=title, aspect="equal")
        fig.savefig(prefix + name + "_quiver.png")
        del fig, ax

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        streams = ax.streamplot(R.T, Z.T, v[:, :, 0].T, v[:, :, 2].T,
                                color=v_mag.T, density=0.9, arrowsize=3,
                                broken_streamlines=False, linewidth=2)
        fig.colorbar(plt.cm.ScalarMappable(
                norm=plt.Normalize(v_mag.min(), v_mag.max()),
                cmap=streams.lines.get_cmap()),
                     ax=ax, label=colorlabel)
        ax.set(xlabel=r"$r$ [m]", ylabel=r"$z$ [m]", title=title, aspect="equal")
        fig.savefig(prefix + name + ".png")
        del fig, ax

    zi = np.argmin(np.abs(z))
    for name, label, title in [("J", r"$J_\phi$ [$\rm A/m^2$]", "Toroidal $J$"), ("B", r"$B_\phi$ [T]", "Toroidal $B$")]:
        v = file.variables[name][:]
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.plot(r, v[:, zi, 1])
        ax.set(xlim=(1.35, 2.15), xlabel="$r$ [m]", ylabel=label, title=title)
        fig.savefig(prefix + name + "phi_axis.png")
        del fig, ax

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        contours = ax.contourf(R, Z, v[:, :, 1], levels=100,
                               cmap='viridis')
        ax.set(xlabel="$r$ [m]", ylabel="$z$ [m]", title=title)
        fig.colorbar(contours, label=label)
        fig.savefig(prefix + name + "phi.png")
        del fig, ax, contours

    for name, colorlabel, title in [("psi", r"$\psi$ [Tm]", r"$\psi$")]:
        v = file.variables[name][:]
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        contours = ax.contourf(R, Z, v, levels=100, cmap='viridis')
        ax.set(xlabel="$r$ [m]", ylabel="$z$ [m]", title=title)
        fig.colorbar(contours, label=colorlabel)
        fig.savefig(prefix + name + ".png")
        del fig, ax, contours
    for name, label, title in [("p", r"$p$ [Pa]", r"$p$"), ("f", r"$F$ [Tm]", r"$F$")]:
        v = file.variables[name][:]
        v_range = file.variables[name + "_range"][:]
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.plot(r, v[:, zi])
        ax.set(xlim=(1.35, 2.15), xlabel="$r$ [m]", ylabel=label, title=title)
        fig.savefig(prefix + name + "_axis.png")
        del fig, ax

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        contours = ax.contourf(R, Z, v, levels=100,
                               cmap='viridis')
        ax.set(xlabel="$r$ [m]", ylabel="$z$ [m]", title=title)
        fig.colorbar(contours, label=label)
        fig.savefig(prefix + name + ".png")
        del fig, ax, contours

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.plot(psi_range, v_range)
        ax.set(xlabel=r"$\psi$ [$\rm Tm^2$]", ylabel=label, title=title)
        fig.savefig(prefix + name + "_range.png")
        del fig, ax
    if ictype.lower() == "hmode":
        p_H = file.variables["p"][:]
        p_ped = file.variables["p_ped"][:]
        p_L = p_H - p_ped
        v_range = file.variables[name + "_range"][:]
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.plot(r, p_H[:, zi], label="H-mode")
        ax.plot(r, p_ped[:, zi], label="Pedestal")
        ax.plot(r, p_L[:, zi], label="L-mode")
        ax.set(xlim=(1.35, 2.15), xlabel="$r$ [m]", ylabel=r"$p$ [Pa]", title="High-confinement mode profile")
        ax.legend()
        fig.savefig(prefix + "p_ped_axis.png")
        del fig, ax
    if ictype.lower() == "diamagnetic":
        B_phi = file.variables["f"][:] / r[:, None]
        B_phi_delta = file.variables["f_delta"][:] / r[:, None]
        B_phi_bg = B_phi - B_phi_delta
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.plot(r, B_phi[:, zi], label="Total")
        ax.plot(r, B_phi_bg[:, zi], label="Background")
        ax.set(xlim=(1.35, 2.15), xlabel="$r$ [m]", ylabel=r"$B_\phi$ [T]", title="Diamagnetic")
        ax.legend()
        fig.savefig(prefix + "B_delta_axis.png")
        del fig, ax
    file.close()


if __name__ == "__main__":
    parser = ArgumentParser(prog="plot.py",
                            description="Convert between file formats.")
    parser.add_argument("-f")
    parser.add_argument("-t")
    args = parser.parse_args()
    while not os.path.exists(".git"):
        if os.path.abspath(os.curdir) == "/":
            raise Exception("Not in a valid project. Create Git repo.")
        os.chdir("..")
    with open("config.toml", "rb") as f:
        config = tomllib.load(f)
        if not args.f:
            fname = config["output"]["name"]
        else:
            fname = args.f
        if not args.t:
            ictype = config["initial_condition"]["type"]
        else:
            ictype = args.t

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 18,
        "figure.dpi": 200
        })

    extract_plots(fname, ictype)
