# GS Solver

Solve for $\psi(r, z)$ in Grad-Shafranov equation:
```math
r\frac{\partial}{\partial r}\left(\frac{1}{r}\frac{\partial \psi}{\partial r}\right) + \frac{\partial^2\psi}{\partial z^2} = -\mu_0 r^2 \frac{dp}{d\psi} - \frac{1}{2} \frac{d F^2}{d\psi}
```
using finite difference method.


## Build

This project uses CMake 4.

Execute the following to compile the executable `solver`
in build directory `./build/`.
```sh
cmake -S . -B build
cmake --build build
```


## Debug

For your convenience, `debug` target
is provided.
```sh
cmake --build build --target debug
```
This will execute `gdb` on the compiled executable `solver`.

Every executable is compiled with debugging information by default.
You can feel free to use what you prefer.
