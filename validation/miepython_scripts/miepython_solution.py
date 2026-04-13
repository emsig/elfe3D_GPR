import numpy as np
import matplotlib.pyplot as plt
import miepython as mie
import csv

# Physical constants
MU0 = 4e-7 * np.pi
EPS0 = 8.854187817e-12
C0 = 1 / np.sqrt(MU0 * EPS0)

def complex_permittivity(eps_r: float, sigma: float, omega: float) -> np.complex128:
    """
    Complex permittivity in a lossy medium.
    """
    return eps_r * EPS0 - 1j * sigma / omega


def wavenumber(eps_complex: np.complex128, omega: float) -> np.complex128:
    """
    Complex wavenumber in a lossy medium.
    """
    return omega * np.sqrt(MU0 * eps_complex)


def dipole_field(r_obs: np.ndarray, r_src: np.ndarray, eps_bg: np.complex128, k: np.complex128, p: np.ndarray) -> np.ndarray:
    """
    Electric field from a Hertzian dipole of finite length.

    r_obs: observation point
    r_src: source location
    eps_bg: background permittivity
    k: complex wavenumber
    p: dipole moment vector
    """

    r = r_obs - r_src
    R = np.linalg.norm(r)

    r_hat = r / R
    kr = k * R

    # E-field terms
    exp_term = np.exp(1j * kr) / (4 * np.pi)
    term1 = k**2 * np.cross(np.cross(r_hat, p), r_hat) / R
    term2 = (1/R**3 - 1j * k / R**2) * (3 * r_hat * np.dot(r_hat, p) - p)

    # Total field
    E = (1/eps_bg) * exp_term * (term1 + term2)

    return E


def spherical_basis(r_hat: np.ndarray) -> np.ndarray:
    """"
    Given a unit vector r_hat, compute the corresponding spherical basis vectors theta_hat and phi_hat.
    """

    z = np.array([0.0, 0.0, 1.0])

    if np.allclose(r_hat, z):
        theta_hat = np.array([1.0, 0.0, 0.0])
    else:
        theta_hat = (z - np.dot(z, r_hat) * r_hat)
        theta_hat /= np.linalg.norm(theta_hat)

    phi_hat = np.cross(r_hat, theta_hat)

    return theta_hat, phi_hat


def mie_scatter_vector(E_inc: np.ndarray, r_hat_scat: np.ndarray, m: np.complex128, k: np.complex128, radius: float, theta: float) -> np.ndarray:
    """
    Compute the scattered electric field vector using Mie scattering amplitudes.
    """

    x = np.real(k) * radius      # Mie size parameter
    mu = np.cos(theta)  # Mie scattering angle parameter

    S1, S2 = mie.S1_S2(m, x, mu)

    theta_hat, phi_hat = spherical_basis(r_hat_scat)

    # Project incident field
    E_theta = np.dot(E_inc, theta_hat)
    E_phi   = np.dot(E_inc, phi_hat)

    # Apply Mie scattering using the S1 and S2 amplitudes
    E_scat = S2 * E_theta * theta_hat + S1 * E_phi * phi_hat

    return E_scat


def propagate(r_from: np.ndarray, r_to: np.ndarray, k: np.complex128) -> np.complex128:
    """
    Green's function for wave propagation in a homogeneous medium.
    """

    r = r_to - r_from
    R = np.linalg.norm(r)
    return np.exp(1j * k * R) / (4 * np.pi * R)


def reflect_point(r: np.ndarray) -> np.ndarray:
    """
    Reflect a point across the z=0 plane.
    """
    return np.array([r[0], r[1], -r[2]])


def fresnel_reflection_coefficient(eps_bg: np.complex128, omega: float, theta_i: float) -> np.complex128:
    """
    Angle-dependent Fresnel reflection (averaged TE + TM contributions @TODO: Check for generalization.)
    """
    eps_air = EPS0

    k1 = omega * np.sqrt(MU0 * eps_bg)
    k2 = omega * np.sqrt(MU0 * eps_air)

    sin_t = (k1/k2)*np.sin(theta_i)
    cos_t = np.sqrt(1 - sin_t**2)
    cos_i = np.cos(theta_i)

    Gamma_TE = (k1*cos_i - k2*cos_t) / (k1*cos_i + k2*cos_t)
    Gamma_TM = (eps_bg*k2*cos_i - eps_air*k1*cos_t) / (eps_bg*k2*cos_i + eps_air*k1*cos_t)

    return 0.5*(Gamma_TE + Gamma_TM)


# ---------------------------------
# Main Mie + Green's function model
# ---------------------------------

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
    include_direct=True,
    include_reflect_tx=False,
    include_reflect_rx=False
):
    """
    Compute the total electric field at receiver locations due to a dipole s
    ource and Mie scattering from a sphere embedded in a lossy medium, 
    including optional paths that don't just involve direct 
    scattering from the sphere (e.g. reflections from the surface).

    No full Sommerfeld integral yet (unsure if needed).
    Polarization mixing simplified.
    """

    results = []

    # Loop over frequencies
    for f in freqs:
        omega = 2*np.pi*f

        eps_bg_c = complex_permittivity(eps_bg, sigma_bg, omega)
        eps_sph_c = complex_permittivity(eps_sph, sigma_sph, omega)

        k = wavenumber(eps_bg_c, omega)
        m = np.sqrt(eps_sph_c / eps_bg_c)   # Mie parameter: relative refractive index

        freq_vals = []

        # Loop over receivers
        for r_r in r_receivers:

            E_total = np.zeros(3, dtype=complex)

            paths = []

            # Direct path from source to sphere to receiver
            if include_direct:
                paths.append((r_s, r_c, r_r, 1.0))

            # TX reflection from ground
            if include_reflect_tx:
                paths.append((reflect_point(r_s), r_c, r_r, "tx"))

            # RX reflection from ground
            if include_reflect_rx:
                paths.append((r_s, r_c, reflect_point(r_r), "rx"))

            for r_s_eff, r_c_eff, r_r_eff, tag in paths:

                # Incident field at sphere
                E_inc = dipole_field(r_c_eff, r_s_eff, eps_bg_c, k, p)

                # Direction to receiver
                r_vec = r_r_eff - r_c_eff
                r_hat = r_vec / np.linalg.norm(r_vec)

                # Scattering angle
                v_inc = r_c_eff - r_s_eff
                v_scat = r_vec

                theta = np.arccos(
                    np.dot(v_inc/np.linalg.norm(v_inc),
                           v_scat/np.linalg.norm(v_scat))
                )

                # Mie scattering
                E_scat = mie_scatter_vector(E_inc, r_hat, m, k, radius, theta)

                # Propagation
                G = propagate(r_c_eff, r_r_eff, k)

                # Fresnel scaling if needed
                if tag == "tx" or tag == "rx":
                    theta_i = theta  # approximate
                    Gamma = fresnel_reflection_coefficient(eps_bg_c, omega, theta_i)
                    E_scat *= Gamma

                # Final contribution
                E_total += 1j * omega * MU0 * G * E_scat

            freq_vals.append(E_total)

        results.append(freq_vals)

    return np.array(results)


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
    include_direct=True,
    include_reflect_tx=False,
    include_reflect_rx=False
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
    include_direct=True,
    include_reflect_tx=False,
    include_reflect_rx=False
)

# E_broadside and E_endfire shape: (1, 48, 3) since one freq, 48 receivers, 3 components

# Save to CSV
# For broadside
with open('broadside_electric_field.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Distance (m)', 'Ex_real', 'Ex_imag', 'Ey_real', 'Ey_imag', 'Ez_real', 'Ez_imag'])
    for i, d in enumerate(distances_broadside):
        row = [d, E_broadside[0, i, 0].real, E_broadside[0, i, 0].imag,
               E_broadside[0, i, 1].real, E_broadside[0, i, 1].imag,
               E_broadside[0, i, 2].real, E_broadside[0, i, 2].imag]
        writer.writerow(row)

# For endfire
with open('endfire_electric_field.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Distance (m)', 'Ex_real', 'Ex_imag', 'Ey_real', 'Ey_imag', 'Ez_real', 'Ez_imag'])
    for i, d in enumerate(distances_endfire):
        row = [d, E_endfire[0, i, 0].real, E_endfire[0, i, 0].imag,
               E_endfire[0, i, 1].real, E_endfire[0, i, 1].imag,
               E_endfire[0, i, 2].real, E_endfire[0, i, 2].imag]
        writer.writerow(row)

print("Data saved to broadside_electric_field.csv and endfire_electric_field.csv")

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
plt.savefig("mie_green_radial_survey.png", dpi=300)
plt.show()
