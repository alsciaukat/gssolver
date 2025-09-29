# GS Solver

Solve for $\psi(r, z)$ in Grad-Shafranov equation:
```math
r\frac{\partial}{\partial r}\left(\frac{1}{r}\frac{\partial \psi}{\partial r}\right) + \frac{\partial^2\psi}{\partial z^2} = -\mu_0 r^2 \frac{dp}{d\psi} - \frac{1}{2} \frac{d g^2}{d\psi}
```
using finite difference method.


## Build

CMake 3.12 or above is required.

> On Windows machine, NetCDF is not supported out of the box.
> The resulting binary would have limited functionalities.

For Linux or Mac, NetCDF is linked using [pkgconf](https://github.com/pkgconf/pkgconf).
In order to check if `pkgconf` is installed, try running
```sh
$ pkgconf --version
> 2.5.1 # This indicate that it is installed correctly.
```
in your shell.

Furthermore, `pkgconf` should be able to find NetCDF4 module.
Check by running the following.
```sh
$ pkgconf --list-all | grep netcdf
> netcdf-cxx4 ... # you should see this.
> netcdf      ...
```

Finally, execute the following to compile the executable `solver`
in build directory `./build/`.
```sh
$ cmake -S . -B build
$ cmake --build build
```


## Debug

For your convenience, `debug` target is provided.
```sh
cmake --build build --target debug
```
This will execute `gdb` on the compiled executable `solver`.

Every executable is compiled with debugging information by default.
Feel free to use whatever debugger you prefer.
