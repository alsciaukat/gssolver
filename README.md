# GS Solver

Solve for $\psi(r, z)$ in Grad-Shafranov equation:
```math
r\frac{\partial}{\partial r}\left(\frac{1}{r}\frac{\partial \psi}{\partial r}\right) + \frac{\partial^2\psi}{\partial z^2} = -\mu_0 r^2 \frac{dp}{d\psi} - \frac{1}{2} \frac{d g^2}{d\psi}
```
using finite difference method.


## Build

In order to compile the code yourself,
you need the following tools to be on your system:
- [CMake](https://cmake.org/) 3.12 or above
- C++ compiler. `g++` (from [GCC](https://gcc.gnu.org/)) or `clang++` (from [Clang](https://clang.llvm.org/)) would work.
- [pkgconf](https://github.com/pkgconf/pkgconf)

You also need these additional C++ libraries:
- NetCDF4
- [toml++](https://github.com/marzer/tomlplusplus)

If you are running Windows or Mac, consult the website for proper instruction for your system.

For Mac users, you might want to use [homebrew](https://brew.sh/).
The package names are `cmake`, `gcc`, `pkgconf`, `netcdf-cxx`, and `tomlplusplus`.

If you are on Linux, I strongly recommend using the package manager from your distribution.
For Debian, Ubuntu and Linux Mint, it is `apt`. For RedHat and Fedora, it is `dnf` (or the legacy `rpm`).


The libraries are linked using [pkgconf](https://github.com/pkgconf/pkgconf).
In order to check if `pkgconf` is installed, try running
```sh
$ pkgconf --version
> 2.5.1 # This indicate that it is installed correctly.
```
in your shell.

Additionally, `pkgconf` should be able to find `netcdf-cxx4` and `toml++` module.
Check by running the following.
```sh
$ pkgconf --list-all | grep netcdf
> netcdf-cxx4  # you should see this.
> netcdf

$ pkgconf --list-all | grep toml
> tomlplusplus
```

If everything is successfully installed,
execute the following to compile the executable `solver` in build directory `./build/`.
```sh
$ cmake -S . -B build
$ cmake --build build
```


## Debug

Every executable is compiled with debugging information by default.
You can use whatever debugger you prefer.

If you have `gdb` on your system,
`debug` target can be used to execute `gdb` on the executable `solver`.
```sh
$ cmake --build build --target debug
```


# Python Tools

Python scripts that visualize and export the data is provided.
They are stored in the `./tools/` directory.


## Setup

The dependencies are specified in the `pyproject.toml` file.

Run
```sh
$ pip install .
```
at the project root (`gssolver/`) to setup the environment.

The required packages will be automatically installed by the above command, along
with their dependencies.
