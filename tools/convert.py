#!/usr/bin/python

from netCDF4 import Dataset
import csv
import os
import tomllib
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser(prog="convert.py",
                            description="Convert between file formats.")
    parser.add_argument("variable")
    parser.add_argument("-f")
    parser.add_argument("-o")
    args = parser.parse_args()
    if (args.f is None):
        while not os.path.exists(".git"):
            if os.path.abspath(os.curdir) == "/":
                raise Exception("Not in a valid project. Create Git repo.")
            os.chdir("..")
        with open("config.toml", "rb") as f:
            config = tomllib.load(f)
            fname = config["output"]["name"]
    else:
        fname = args.f

    if (args.o is None):
        oname = "result.csv"
    else:
        oname = args.o
    vname = args.variable

    file = Dataset(fname, "r")
    v_nc = file.variables[vname]
    v = v_nc[:]
    f = open(oname, "w")
    writer = csv.writer(f, delimiter="\t")

    if len(v_nc.dimensions) == 1:
        axname, = v_nc.dimensions
        ax = file.variables[axname][:]
        writer.writerow([axname, vname])
        for axi, vi in zip(ax, v):
            writer.writerow([axi, vi])
    elif len(v_nc.dimensions) == 2:
        ax1name, ax2name = v_nc.dimensions
        ax1 = file.variables[ax1name][:]
        ax2 = file.variables[ax2name][:]
        writer.writerow([ax1name, ax2name, vname])
        for i in range(len(file.dimensions[ax1name])):
            for j in range(len(file.dimensions[ax2name])):
                writer.writerow([ax1[i], ax2[j], v[i, j]])
    elif len(v_nc.dimensions) == 3:
        ax1name, ax2name, _ = v_nc.dimensions
        ax1 = file.variables[ax1name][:]
        ax2 = file.variables[ax2name][:]
        writer.writerow([ax1name, ax2name, vname])
        for i in range(len(file.dimensions[ax1name])):
            for j in range(len(file.dimensions[ax2name])):
                writer.writerow([ax1[i], ax2[j], v[i, j, 0], v[i, j, 1], v[i, j, 2]])
