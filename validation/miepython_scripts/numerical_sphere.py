import numpy as np

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
MU0 = 4e-7 * np.pi
EPS0 = 8.854187817e-12

# -----------------------------------------------------------------------------
# Complex permittivity (e^{+jωt})
# -----------------------------------------------------------------------------
def complex_permittivity(eps_r, sigma, omega):
    return eps_r * EPS0 - 1j * sigma / omega


def wavenumber(eps_c, omega):
    k = omega * np.sqrt(MU0 * eps_c)
    if np.imag(k) > 0:
        k = -k
    return k


# -----------------------------------------------------------------------------
# Dipole field (background field)
# -----------------------------------------------------------------------------
def dipole_field(r_obs, r_src, eps_bg, k, p):
    r = r_obs - r_src
    R = np.linalg.norm(r)
    r_hat = r / R
    kr = k * R

    exp_term = np.exp(-1j * kr) / (4 * np.pi)

    term1 = k**2 * np.cross(np.cross(r_hat, p), r_hat) / R
    term2 = (1/R**3 + 1j * k / R**2) * (3*r_hat*np.dot(r_hat, p) - p)

    return (1/eps_bg) * exp_term * (term1 + term2)


# -----------------------------------------------------------------------------
# Dyadic Green's function (homogeneous version)
# Replace THIS with Sommerfeld for half-space
# -----------------------------------------------------------------------------
def green_dyadic(r, rp, k):
    R_vec = r - rp
    R = np.linalg.norm(R_vec)
    r_hat = R_vec / R

    exp_term = np.exp(-1j * k * R) / (4 * np.pi * R)

    I = np.eye(3)

    term1 = (1 + 1j/(k*R) - 1/(k*R)**2) * I
    term2 = (-1 - 3j/(k*R) + 3/(k*R)**2) * np.outer(r_hat, r_hat)

    return exp_term * (term1 + term2)


# -----------------------------------------------------------------------------
# Discretize sphere into volume elements
# -----------------------------------------------------------------------------
def discretize_sphere(center, radius, N):
    """
    Simple cubic voxel discretization inside sphere
    """
    points = []
    xs = np.linspace(-radius, radius, N)

    for x in xs:
        for y in xs:
            for z in xs:
                r = np.array([x, y, z])
                if np.linalg.norm(r) <= radius:
                    points.append(center + r)

    points = np.array(points)
    dv = (2*radius/N)**3
    return points, dv


# -----------------------------------------------------------------------------
# Solve volume integral equation
# -----------------------------------------------------------------------------
def solve_vie(points, dv, r_s, p, eps_bg_c, eps_sph_c, k, omega):
    """
    Solve for E inside sphere
    """
    N = len(points)
    chi = 1j * omega * (eps_sph_c - eps_bg_c)

    A = np.zeros((3*N, 3*N), dtype=complex)
    b = np.zeros(3*N, dtype=complex)

    # Incident field
    for i, ri in enumerate(points):
        Ei = dipole_field(ri, r_s, eps_bg_c, k, p)
        b[3*i:3*i+3] = Ei

    # Build system
    for i, ri in enumerate(points):
        for j, rj in enumerate(points):

            if i == j:
                A[3*i:3*i+3, 3*j:3*j+3] = np.eye(3)
            else:
                G = green_dyadic(ri, rj, k)
                A[3*i:3*i+3, 3*j:3*j+3] -= chi * G * dv

    # Solve
    E_vec = np.linalg.solve(A, b)

    return E_vec.reshape((N, 3))


# -----------------------------------------------------------------------------
# Compute anomaly field at receiver
# -----------------------------------------------------------------------------
def compute_anomaly_field(r_obs, points, E_inside, dv, eps_bg_c, eps_sph_c, k, omega):
    chi = 1j * omega * (eps_sph_c - eps_bg_c)

    E = np.zeros(3, dtype=complex)

    for ri, Ei in zip(points, E_inside):
        G = green_dyadic(r_obs, ri, k)
        E += G @ (chi * Ei) * dv

    return E


# -----------------------------------------------------------------------------
# Main driver
# -----------------------------------------------------------------------------
def compute_full_field(
    freqs,
    r_s,
    r_c,
    r_receivers,
    p,
    eps_bg,
    sigma_bg,
    eps_sph,
    sigma_sph,
    radius,
    N_vox=6
):
    results = []

    for f in freqs:
        omega = 2*np.pi*f

        eps_bg_c = complex_permittivity(eps_bg, sigma_bg, omega)
        eps_sph_c = complex_permittivity(eps_sph, sigma_sph, omega)

        k = wavenumber(eps_bg_c, omega)

        # Discretize sphere
        points, dv = discretize_sphere(r_c, radius, N_vox)

        # Solve internal field
        E_inside = solve_vie(points, dv, r_s, p, eps_bg_c, eps_sph_c, k, omega)

        freq_vals = []

        for r_r in r_receivers:
            E = compute_anomaly_field(
                r_r, points, E_inside, dv,
                eps_bg_c, eps_sph_c, k, omega
            )
            freq_vals.append(E)

        results.append(freq_vals)

    return np.array(results)


from scipy.special import j0

def green_dyadic_halfspace(r, rp, k1, k2, omega):
    """
    Dyadic Green's function for half-space:
    region 1 (z > 0): air
    region 2 (z < 0): ground

    r  = observation point
    rp = source point

    k1 = wavenumber in upper half-space (air)
    k2 = wavenumber in lower half-space (ground)

    Assumes both r and rp are in lower half-space (typical GPR case).
    """

    MU0 = 4e-7 * np.pi

    x, y, z = r
    xp, yp, zp = rp

    rho = np.sqrt((x - xp)**2 + (y - yp)**2)
    z_sum = z + zp   # reflection path

    # --- Direct term (homogeneous ground) ---
    R_vec = r - rp
    R = np.linalg.norm(R_vec)
    r_hat = R_vec / R

    exp_term = np.exp(-1j * k2 * R) / (4*np.pi*R)

    I = np.eye(3)

    term1 = (1 + 1j/(k2*R) - 1/(k2*R)**2) * I
    term2 = (-1 - 3j/(k2*R) + 3/(k2*R)**2) * np.outer(r_hat, r_hat)

    G_direct = exp_term * (term1 + term2)

    # --- Sommerfeld integral (reflected + lateral) ---
    # numerical quadrature over k_rho

    def integrand(k_rho):
        kz2 = np.sqrt(k2**2 - k_rho**2 + 0j)
        kz1 = np.sqrt(k1**2 - k_rho**2 + 0j)

        # choose branch with Im(kz) >= 0
        if np.imag(kz2) < 0: kz2 = -kz2
        if np.imag(kz1) < 0: kz1 = -kz1

        # Fresnel coefficients
        R_TE = (kz2 - kz1) / (kz2 + kz1)
        R_TM = (k1**2 * kz2 - k2**2 * kz1) / (k1**2 * kz2 + k2**2 * kz1)

        J = j0(k_rho * rho)
        phase = np.exp(-1j * kz2 * z_sum)

        return k_rho * J * phase * (R_TE + R_TM) / (2 * kz2)

    # integrate numerically
    kmax = 10 * abs(k2)
    Nk = 200

    k_rhos = np.linspace(0, kmax, Nk)
    dk = k_rhos[1] - k_rhos[0]

    integral = 0.0 + 0j
    for kr in k_rhos:
        integral += integrand(kr) * dk

    G_ref_scalar = integral / (2*np.pi)

    # approximate dyadic form (dominant transverse behavior)
    G_ref = G_ref_scalar * I

    return G_direct + G_ref


# -----------------------------
# Usage
# -----------------------------

# frequency
freq = 100e6  # 100 MHz
freqs = [freq]

# Geometry (meters)
r_s = np.array([0.0, 0.0, 0.025])   # Dipole source 2.5 cm above the surface
r_c = np.array([0.0, 0.0, -0.7])    # @TODO: Need to fix to match elfe3D_GPR later.

# Medium properties
eps_bg = 4.0
sigma_bg = 1e-4

# Sphere properties (anomaly)
eps_sph = 20.0
sigma_sph = 1e-4

radius = (3e8 / freq) / 16
p = np.array([1.0, 0.0, 0.0]) # Calibrate to elfe3D_GPR later.

# Receivers for broadside: along y from 0.1 to 1 m, 48 points
distances_broadside = np.linspace(0.1, 1.0, 48)
r_receivers_broadside = [np.array([0.0, d, -0.0001]) for d in distances_broadside]

# Receivers for endfire: along y from 0.1 to 1 m at x=0.5, 48 points
distances_endfire = np.linspace(0.1, 1.0, 48)
r_receivers_endfire = [np.array([d, 0.0, -0.0001]) for d in distances_endfire]

# Compute for broadside
E_broadside = compute_full_field(
    freqs,
    r_s,
    r_c,
    r_receivers_broadside,
    p,
    eps_bg,
    sigma_bg,
    eps_sph,
    sigma_sph,
    radius,
)

# Compute for endfire
E_endfire = compute_full_field(
    freqs,
    r_s,
    r_c,
    r_receivers_endfire,
    p,
    eps_bg,
    sigma_bg,
    eps_sph,
    sigma_sph,
    radius,
)

# E_broadside and E_endfire shape: (1, 48, 3) since one freq, 48 receivers, 3 components
import csv
from matplotlib import pyplot as plt
# Save to CSV
# For broadside
with open('broadside_electric_field_num.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Distance (m)', 'Ex_real', 'Ex_imag', 'Ey_real', 'Ey_imag', 'Ez_real', 'Ez_imag'])
    for i, d in enumerate(distances_broadside):
        row = [d, E_broadside[0, i, 0].real, E_broadside[0, i, 0].imag,
               E_broadside[0, i, 1].real, E_broadside[0, i, 1].imag,
               E_broadside[0, i, 2].real, E_broadside[0, i, 2].imag]
        writer.writerow(row)

# For endfire
with open('endfire_electric_field_num.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Distance (m)', 'Ex_real', 'Ex_imag', 'Ey_real', 'Ey_imag', 'Ez_real', 'Ez_imag'])
    for i, d in enumerate(distances_endfire):
        row = [d, E_endfire[0, i, 0].real, E_endfire[0, i, 0].imag,
               E_endfire[0, i, 1].real, E_endfire[0, i, 1].imag,
               E_endfire[0, i, 2].real, E_endfire[0, i, 2].imag]
        writer.writerow(row)

print("Data saved to broadside_electric_field_num.csv and endfire_electric_field_num.csv")

# Create three subplots, one for each measurement (inline and endfire)

# Common styling (match class defaults)
label_fs = 18/2
tick_fs = 18/2
legend_fs = 15/2
title_fs = 18/2
lw = 2.5

component_names = ['|E_x|', '|E_y|', '|E_z|']

# Create 3x2 subplot grid (3 components x 2 configurations)
fig, axes = plt.subplots(3, 2, figsize=(14, 12))

component_names = ['|E_x|', '|E_y|', '|E_z|']

for component in range(3):
    # Broadside configuration
    ax = axes[component, 0]
    ax.plot(
        distances_broadside,
        np.abs(E_broadside[0, :, component]),
        linewidth=lw,
        color='C0'
    )
    ax.set_xlabel("Distance (m)", fontsize=label_fs)
    ax.set_ylabel("Electric Field Magnitude", fontsize=label_fs)
    ax.set_title(f"Broadside: {component_names[component]}", fontsize=title_fs, fontweight="bold")
    ax.tick_params(labelsize=tick_fs)
    ax.grid(True, linestyle="--", linewidth=0.5)

    # Endfire configuration
    ax = axes[component, 1]
    ax.plot(
        distances_endfire,
        np.abs(E_endfire[0, :, component]),
        linewidth=lw,
        color='C1'
    )
    ax.set_xlabel("Distance (m)", fontsize=label_fs)
    ax.set_ylabel("Electric Field Magnitude", fontsize=label_fs)
    ax.set_title(f"Endfire: {component_names[component]}", fontsize=title_fs, fontweight="bold")
    ax.tick_params(labelsize=tick_fs)
    ax.grid(True, linestyle="--", linewidth=0.5)

plt.tight_layout()
plt.savefig("numerical_green_radial_survey.png", dpi=300)
plt.show()
