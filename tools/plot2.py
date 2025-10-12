import netCDF4 as nc
import numpy as np
import math
import matplotlib.pyplot as plt

ds = nc.Dataset('hmode.nc', "r")
ds_sol = nc.Dataset('output_solovev.nc', "r")
ds_poly = nc.Dataset('output_polynomial.nc', "r")
print(ds)

"""
# Let's calculate & plot B_z, B_r, B_theta
# B_r, B_z are easy to calculrate, but B_theta is not.
# B_r = - 1/r * dpsi/dz
# B_z = 1/r * dpsi/dr
# B_theta = 1/r * g(psi)
"""

r = ds.variables['r'][:]
z = ds.variables['z'][:]  

# test
psi = ds.variables['psi'][:]

J_r = ds.variables["J_r"][:]
J_z = ds.variables["J_z"][:]
J_phi = ds.variables['J_phi'][:]

B_r = ds.variables["B_r"][:]
B_z = ds.variables["B_z"][:]
B_phi = ds.variables['B_phi'][:]

# solov'ev
psi_sol = ds_sol.variables['psi'][:]

J_r_sol = ds_sol.variables["J_r"][:]
J_z_sol = ds_sol.variables["J_z"][:]
J_phi_sol = ds_sol.variables['J_phi'][:]

B_r_sol = ds_sol.variables["B_r"][:]
B_z_sol = ds_sol.variables["B_z"][:]
B_phi_sol = ds_sol.variables['B_phi'][:]

# polynomial
psi_poly = ds_poly.variables['psi'][:]

J_r_poly = ds_poly.variables["J_r"][:]
J_z_poly = ds_poly.variables["J_z"][:]
J_phi_poly = ds_poly.variables['J_phi'][:]

B_r_poly = ds_poly.variables["B_r"][:]
B_z_poly = ds_poly.variables["B_z"][:]
B_phi_poly = ds_poly.variables['B_phi'][:]

try:
    R = ds.getncattr("R")
    a = ds.getncattr("a")
    b = ds.getncattr("b")
    n = int( ds.getncattr("n") )
    m = int( ds.getncattr("m") )
    type = ds.getncattr("type")
    beta0 = ds.getncattr("beta0")
except:
    print('use built-in parameters')
    R = 1.8
    a = 0.5
    b = 0.5
    n = 1
    m = 2
    type = 'polynomial'
    beta0 = 0.5

mu0 = 4 * np.pi * 1e-7
print('shape of psi:', psi.shape)

"""
for solov'ev solution, gg_prime = - b * R**2 (constant)
So we can determine g:
g(psi) = sqrt(c - 2 * b * R**2 * psi)
where c is integral constant

for polynomial solution, gg_prime = A * (1-psi**m)**n
where A = (1-beta0)*mu0*R
general solution for g is
sqrt[ 2A*sum_(k=0)^n{ binom(n, k) * (-1)^k / (1+mr) * psi^(mr+1) } + C]
where C is arbitrary constant

TODO: How to determine constant c?
"""

def calc_g(psi, type, b, R, c):
    if type == "solovev":
        return np.sqrt(c - 2*b*R**2 * psi)
    elif type == "polynomial":
        g = np.zeros_like(psi)
        for k in range(n):
            g += math.comb(n, k) * (-1)**k / (1 + m*k) * psi**(m*k+1)
        A = 2*(1-beta0) * mu0 * R
        return np.sqrt( A * g + c)
    
def plot2D(f, r=None, z=None, title=None):
    fig, ax = plt.subplots()
    if r is None or z is None:
        plot = ax.imshow(f.T)
    else:
        r, z = np.meshgrid(r, z, indexing='ij')
        plot = ax.pcolormesh(r, z, f[:-1, :-1])
    ax.set_title(title)
    fig.colorbar(plot, ax=ax, fraction=0.04, pad=0.02)
    plt.show()

def plots2D(fs, r=None, z=None, titles=None):
    n = len(fs)
    fig, axes = plt.subplots(ncols=n, figsize=(4*n, 3))

    R, Z = np.meshgrid(r, z, indexing='ij')
    for ax, f, title in zip(axes, fs, titles):
        if r is None or z is None:
            plot = ax.imshow(f.T)
        else:
            plot = ax.pcolormesh(R, Z, f[:-1, :-1])
        ax.set_title(title)
        ax.set_aspect('equal'); ax.set_xlabel('r'); ax.set_ylabel('z')
        
    fig.colorbar(plot, ax=axes, fraction=0.04, pad=0.02)
    plt.show()

def calc_B(psi, r, z, g=None):
    dpdz = np.gradient(psi, z, axis=1)
    dpdr = np.gradient(psi, r, axis=0)

    R, Z = np.meshgrid(r, z, indexing='ij')
    B_r = - 1/R * dpdz
    B_z = + 1/R * dpdr
    if g is not None:
        B_theta = 1/R * g
        return B_r, B_z, B_theta
    else:
        return B_r, B_z, np.zeros_lie(B_r)

def plot_vector(B_r, B_z, B_theta=None, r=None, z=None, step=5, mode='stream', title=None ):
    if r is None or z is None: # if r, z is not provided, make it
        r = np.arange(B_z.shape[0])
        z = np.arange(B_r.shape[1])

    # too many point is not good!
    r, z = r[::step], z[::step]
    B_r, B_z = B_r[::step, ::step], B_z[::step, ::step]

    B_r_T, B_z_T = B_r.T, B_z.T

    
    if B_theta is None: # draw at rz plane
        fig, ax = plt.subplots() 
        if mode == 'stream':
            mag = np.hypot(B_r_T, B_z_T)
            sp = ax.streamplot(r, z, B_r_T, B_z_T, color=mag, density=1, broken_streamlines=False, linewidth=1)
            fig.colorbar(plt.cm.ScalarMappable(
                norm=plt.Normalize(mag.min(), mag.max()), cmap=sp.lines.get_cmap()),
                ax=ax, label='|B|')
        else:  # quiver
            X, Y = np.meshgrid(r, z)
            q = ax.quiver(X, Y, B_r_T, B_z_T, np.hypot(B_r_T, B_z_T), pivot='mid')
            fig.colorbar(q, ax=ax, label='|B|')
    else:
        ax = plt.figure().add_subplot(projection='3d')

        X, Y, Z = np.meshgrid(r, z, 0)
        B_theta_T = B_theta[::step, ::step].T
        B_r_T, B_z_T, B_theta_T = np.expand_dims(B_r_T, axis=-1), np.expand_dims(B_z_T, axis=-1), np.expand_dims(B_theta_T, axis=-1)
        print(B_r.shape)
        print(X.shape)
        ax.quiver(X, Y, Z, B_r_T, B_z_T, B_theta_T, length=0.05, normalize=True)

    ax.set_aspect('equal'); ax.set_xlabel('r'); ax.set_ylabel('z')
    if title: ax.set_title(title)
    plt.show()
    return fig, ax

if __name__ == "__main__":
    g = calc_g(psi, type=type, b=b, R=R, c=0)
    B_r, B_z, B_theta = calc_B(psi=psi, r=r, z=z, g=g)
    plot2D(psi, r, z, title='psi')
    plots2D([J_phi, J_phi_sol, J_phi_poly], r=r, z=z, titles=['test', 'solovev', 'polynomial'])
    # plots2D([B_r, B_z, B_theta], r=r, z=z, titles=['B_r', 'B_z', 'B_theta'])
    # plot_vector(B_r, B_z, r=r, z=z, mode='qiuver')
    # plot_vector(B_r, B_z, B_theta=B_theta, r=r, z=z, mode='qiuver')
