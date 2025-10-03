import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

ds = nc.Dataset('result.nc', "r")

# Let's calculate B_z, B_r
# B_theta is a little treaky.

# B_r = - 1/r * dpsi/dz
# B_z = 1/r * dpsi/dr
# B_theta = 1/r * F(psi)

psi = ds.variables['psi'][:] # axis=0 is z. axis=1 is r
r = ds.variables['r'][:]
z = ds.variables['z'][:]    

print(psi.shape)

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
    fig.colorbar(plot, ax=axes, fraction=0.04, pad=0.02)
    plt.show()

def calc_B(psi, r, z):
    dpdz = np.gradient(psi, z, axis=1)
    dpdr = np.gradient(psi, r, axis=0)

    R, Z = np.meshgrid(r, z, indexing='ij')
    B_r = - 1/R * dpdz
    B_z = + 1/R * dpdr
    return B_r, B_z

def plot_B(B_r, B_z, B_theta=None, r=None, z=None, step=10, mode='stream', title=None ):
    if r is None or z is None:
        r = np.arange(B_z.shape[0])
        z = np.arange(B_r.shape[1])
    
    # too many point is not good!
    r, z = r[::step], z[::step]
    B_r, B_z = B_r[::step, ::step], B_z[::step, ::step]

    B_r_T, B_z_T = B_r.T, B_z.T

    fig, ax = plt.subplots()
    if B_theta is None: # draw at rz plane 
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
        pass # 3D plot is not supported for now

    ax.set_aspect('equal'); ax.set_xlabel('r'); ax.set_ylabel('z')
    if title: ax.set_title(title)
    plt.show()
    return fig, ax

if __name__ == "__main__":
    B_r, B_z = calc_B(psi, r, z)
    # plot2D(psi, r, z, title='psi')
    # plots2D([B_r, B_z], r=r, z=z, titles=['B_r', 'B_z'])
    print(r.shape)
    plot_B(B_r, B_z, r=r, z=z, mode='qiuver')
