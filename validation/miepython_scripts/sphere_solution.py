import csv
import math

import numpy as np
from scipy.special import lpmv, spherical_jn, hankel1

# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------
MU0 = 4e-7 * np.pi
EPS0 = 8.854187817e-12
C0 = 1.0 / np.sqrt(MU0 * EPS0)


# -----------------------------------------------------------------------------
# Convention: phasors ~ exp(+j omega t)
# Passive media should use exp(-j k r) with Im(k) <= 0
# -----------------------------------------------------------------------------
def complex_permittivity(eps_r: float, sigma: float, omega: float) -> np.complex128:
    """
    Complex permittivity for exp(+j omega t):
        eps_tilde = eps0 * eps_r - j sigma / omega
    """
    return np.complex128(eps_r * EPS0 - 1j * sigma / omega)


def robust_wavenumber(eps_complex: np.complex128, omega: float) -> np.complex128:
    """
    Robust complex wavenumber for exp(+j omega t), using
        k = omega * sqrt(mu0 * eps_tilde)

    The branch is chosen so that Im(k) <= 0, giving exp(-j k r) decay.
    """
    k = omega * np.sqrt(MU0 * eps_complex)
    if (np.imag(k) > 0) or (abs(np.imag(k)) < 1e-14 and np.real(k) < 0):
        k = -k
    return np.complex128(k)


# -----------------------------------------------------------------------------
# Geometry helpers
# -----------------------------------------------------------------------------
def normalize(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Cannot normalize zero vector.")
    return v / n


def reflect_point(r: np.ndarray) -> np.ndarray:
    """Reflect a point across the z = 0 plane."""
    return np.array([r[0], r[1], -r[2]], dtype=float)


def cartesian_to_spherical(v: np.ndarray):
    """
    Return (r, theta, phi) for vector v in Cartesian coordinates.
    theta: polar angle from +z
    phi: azimuth from +x toward +y
    """
    r = np.linalg.norm(v)
    if r == 0:
        raise ValueError("Spherical coordinates undefined at the origin.")
    x, y, z = v
    theta = np.arccos(np.clip(z / r, -1.0, 1.0))
    phi = np.arctan2(y, x)
    return r, theta, phi


def spherical_unit_vectors(theta: float, phi: float):
    """
    Standard spherical basis vectors. This handles both +z and -z without
    special casing.
    """
    st = np.sin(theta)
    ct = np.cos(theta)
    sp = np.sin(phi)
    cp = np.cos(phi)

    r_hat = np.array([st * cp, st * sp, ct], dtype=float)
    theta_hat = np.array([ct * cp, ct * sp, -st], dtype=float)
    phi_hat = np.array([-sp, cp, 0.0], dtype=float)
    return r_hat, theta_hat, phi_hat


def spherical_basis(r_hat: np.ndarray):
    """
    Given a unit vector r_hat, compute spherical basis vectors theta_hat and phi_hat.
    Handles both +z and -z.
    """
    r_hat = normalize(r_hat)
    _, theta, phi = cartesian_to_spherical(r_hat)
    _, theta_hat, phi_hat = spherical_unit_vectors(theta, phi)
    return theta_hat, phi_hat


# -----------------------------------------------------------------------------
# Exact dipole field in a homogeneous lossy medium, exp(+j omega t)
# -----------------------------------------------------------------------------
def dipole_field(
    r_obs: np.ndarray,
    r_src: np.ndarray,
    eps_bg: np.complex128,
    k: np.complex128,
    p: np.ndarray,
) -> np.ndarray:
    """
    Electric field from a Hertzian electric dipole in a homogeneous medium.

    Convention: exp(+j omega t), outgoing wave ~ exp(-j k R).
    This is the exact full near/intermediate/far field expression in a
    homogeneous medium.

    p is the electric dipole moment vector (or an equivalent calibration vector).
    """
    r = r_obs - r_src
    R = np.linalg.norm(r)
    if R == 0:
        raise ValueError("Dipole field singular at the source location.")

    r_hat = r / R
    kr = k * R

    # Outgoing factor for exp(+j omega t)
    exp_term = np.exp(-1j * kr) / (4.0 * np.pi)

    # Standard dyadic Green's function form adapted to exp(+j omega t)
    term1 = k**2 * np.cross(np.cross(r_hat, p), r_hat) / R
    term2 = (1.0 / R**3 + 1j * k / R**2) * (3.0 * r_hat * np.dot(r_hat, p) - p)

    E = (1.0 / eps_bg) * exp_term * (term1 + term2)
    return E


def propagate(r_from: np.ndarray, r_to: np.ndarray, k: np.complex128) -> np.complex128:
    """
    Scalar outgoing Green's function factor. Kept for compatibility,
    but not used in the rigorous sphere-scattering path.
    """
    R = np.linalg.norm(r_to - r_from)
    return np.exp(-1j * k * R) / (4.0 * np.pi * R)


# -----------------------------------------------------------------------------
# Spherical harmonics and vector spherical wave functions (VSWFs)
# -----------------------------------------------------------------------------
def spherical_harmonic_Y(n: int, m: int, theta: float, phi: float) -> complex:
    """
    Normalized complex spherical harmonic Y_n^m(theta, phi).
    Uses scipy's sph_harm convention.
    """
    if m < 0:
        mp = -m
        return ((-1) ** mp) * np.conj(spherical_harmonic_Y(n, mp, theta, phi))

    x = np.clip(np.cos(theta), -1.0 + 1e-15, 1.0 - 1e-15)
    norm = np.sqrt(
        (2 * n + 1) / (4.0 * np.pi) * math.factorial(n - m) / math.factorial(n + m)
    )
    Pnm = lpmv(m, n, x)
    return norm * Pnm * np.exp(1j * m * phi)


def dY_dtheta(n: int, m: int, theta: float, phi: float) -> complex:
    """
    Derivative dY_n^m / dtheta.
    Uses associated Legendre identities; exact for our numerical purposes.
    """
    if m < 0:
        mp = -m
        return ((-1) ** mp) * np.conj(dY_dtheta(n, mp, theta, phi))

    x = np.clip(np.cos(theta), -1.0 + 1e-15, 1.0 - 1e-15)
    if n == 0:
        return 0.0 + 0.0j

    Pnm = lpmv(m, n, x)
    Pn1m = lpmv(m, n - 1, x) if n - 1 >= 0 else 0.0
    denom = x * x - 1.0
    dPdx = (n * x * Pnm - (n + m) * Pn1m) / denom

    norm = np.sqrt(
        (2 * n + 1) / (4.0 * np.pi) * math.factorial(n - m) / math.factorial(n + m)
    )
    return norm * np.exp(1j * m * phi) * (-np.sin(theta) * dPdx)


def spherical_jn_derivative(n: int, z: complex) -> complex:
    """Derivative of spherical Bessel j_n(z) via recurrence."""
    if n == 0:
        return -spherical_jn(1, z)
    return spherical_jn(n - 1, z) - (n + 1) / z * spherical_jn(n, z)


def spherical_hankel1(n: int, z: complex) -> complex:
    """Spherical Hankel function of the first kind h_n^(1)(z)."""
    return np.sqrt(np.pi / (2.0 * z)) * hankel1(n + 0.5, z)


def spherical_hankel1_derivative(n: int, z: complex) -> complex:
    """Derivative of spherical Hankel h_n^(1)(z) via recurrence."""
    if n == 0:
        return -spherical_hankel1(1, z)
    return spherical_hankel1(n - 1, z) - (n + 1) / z * spherical_hankel1(n, z)


def vswf_MN(
    n: int,
    m: int,
    k: np.complex128,
    r: float,
    theta: float,
    phi: float,
    kind: str = "regular",
):
    """
    Return the vector spherical wave functions M_n^m and N_n^m in Cartesian form.

    kind:
        "regular" -> spherical Bessel j_n
        "outgoing" -> spherical Hankel h_n^(1)

    This formulation is what makes the dipole-to-sphere coupling rigorous.
    """
    r_hat, theta_hat, phi_hat = spherical_unit_vectors(theta, phi)
    kr = k * r

    if kind == "regular":
        z = spherical_jn(n, kr)
        dz = spherical_jn_derivative(n, kr)
    elif kind == "outgoing":
        z = spherical_hankel1(n, kr)
        dz = spherical_hankel1_derivative(n, kr)
    else:
        raise ValueError("kind must be 'regular' or 'outgoing'.")

    Y = spherical_harmonic_Y(n, m, theta, phi)
    dY = dY_dtheta(n, m, theta, phi)

    sin_theta = max(np.sin(theta), 1e-14)
    nfac = np.sqrt(n * (n + 1.0))

    # Tangential vector spherical harmonics
    X = ((1j * m / sin_theta) * Y * theta_hat - dY * phi_hat) / nfac
    Psi = (dY * theta_hat + (1j * m / sin_theta) * Y * phi_hat) / nfac

    # M is purely tangential; N has radial + tangential parts
    M = z * X
    N = (n * (n + 1.0) * z / kr) * Y * r_hat + ((z + kr * dz) / kr) * Psi
    return M, N


# -----------------------------------------------------------------------------
# Sphere Mie coefficients (exact for a sphere in a homogeneous background)
# -----------------------------------------------------------------------------
def mie_coefficients_sphere(
    k_bg: np.complex128,
    k_sph: np.complex128,
    radius: float,
    nmax: int,
):
    """
    Standard sphere Mie coefficients a_n, b_n using complex k and complex x.
    These are exact for an isotropic sphere in a homogeneous background.

    Note: no taking Re(k) here; that would throw away attenuation.
    """
    x = k_bg * radius
    mx = k_sph * radius
    m_rel = k_sph / k_bg

    an = np.zeros(nmax + 1, dtype=complex)
    bn = np.zeros(nmax + 1, dtype=complex)

    for n in range(1, nmax + 1):
        psi_x = x * spherical_jn(n, x)
        psi_mx = mx * spherical_jn(n, mx)

        dpsi_x = spherical_jn(n, x) + x * spherical_jn_derivative(n, x)
        dpsi_mx = spherical_jn(n, mx) + mx * spherical_jn_derivative(n, mx)

        xi_x = x * spherical_hankel1(n, x)
        dxi_x = spherical_hankel1(n, x) + x * spherical_hankel1_derivative(n, x)

        an[n] = (
            m_rel * psi_mx * dpsi_x - psi_x * dpsi_mx
        ) / (
            m_rel * psi_mx * dxi_x - xi_x * dpsi_mx
        )

        bn[n] = (
            psi_mx * dpsi_x - m_rel * psi_x * dpsi_mx
        ) / (
            psi_mx * dxi_x - m_rel * xi_x * dpsi_mx
        )

    return an, bn


# -----------------------------------------------------------------------------
# Plane-interface Fresnel helper (polarization-resolved; still ray-based)
# -----------------------------------------------------------------------------
def fresnel_coefficients(
    eps_inc: np.complex128,
    eps_trn: np.complex128,
    omega: float,
    theta_i: float,
):
    """
    TE/TM Fresnel reflection coefficients using impedance form.

    This is polarization-resolved (not averaged).
    Still a ray-based correction; exact half-space physics would require
    Sommerfeld integrals / layered Green's functions.
    """
    k1 = omega * np.sqrt(MU0 * eps_inc)
    k2 = omega * np.sqrt(MU0 * eps_trn)

    sin_t = (k1 / k2) * np.sin(theta_i)
    cos_t = np.lib.scimath.sqrt(1.0 - sin_t**2)
    cos_i = np.cos(theta_i)

    eta1 = np.sqrt(MU0 / eps_inc)
    eta2 = np.sqrt(MU0 / eps_trn)

    Gamma_TE = (eta2 * cos_i - eta1 * cos_t) / (eta2 * cos_i + eta1 * cos_t)
    Gamma_TM = (eta2 * cos_t - eta1 * cos_i) / (eta2 * cos_t + eta1 * cos_i)

    return Gamma_TE, Gamma_TM


def apply_planar_reflection_to_vector(
    E_vec: np.ndarray,
    k_hat_inc: np.ndarray,
    eps_inc: np.complex128,
    eps_trn: np.complex128,
    omega: float,
    normal: np.ndarray = np.array([0.0, 0.0, 1.0]),
) -> np.ndarray:
    """
    Decompose a vector field into TE/TM with respect to the plane of incidence,
    apply Fresnel coefficients, and reconstruct.

    This is a polarization-correct ray approximation, not an exact half-space
    solution.
    """
    n_hat = normalize(np.asarray(normal, dtype=float))
    k_hat_inc = normalize(np.asarray(k_hat_inc, dtype=float))

    s_hat = np.cross(n_hat, k_hat_inc)
    s_norm = np.linalg.norm(s_hat)
    if s_norm < 1e-12:
        # Near normal incidence: choose an arbitrary transverse direction
        tmp = np.array([1.0, 0.0, 0.0]) if abs(np.dot(n_hat, [1.0, 0.0, 0.0])) < 0.9 else np.array([0.0, 1.0, 0.0])
        s_hat = normalize(np.cross(n_hat, tmp))
    else:
        s_hat = s_hat / s_norm

    p_hat = normalize(np.cross(s_hat, k_hat_inc))

    # Angle of incidence from the normal
    cos_theta_i = np.clip(abs(np.dot(-k_hat_inc, n_hat)), 0.0, 1.0)
    theta_i = np.arccos(cos_theta_i)

    Gamma_TE, Gamma_TM = fresnel_coefficients(eps_inc, eps_trn, omega, theta_i)

    E_s = np.dot(E_vec, s_hat)
    E_p = np.dot(E_vec, p_hat)

    # Reconstruct with reflected TE/TM amplitudes
    return Gamma_TE * E_s * s_hat + Gamma_TM * E_p * p_hat


# -----------------------------------------------------------------------------
# Numerical VSWF projection of an arbitrary source field onto a sphere
# -----------------------------------------------------------------------------
def incident_field_on_surface(
    point: np.ndarray,
    source_pos: np.ndarray,
    p: np.ndarray,
    eps_bg_c: np.complex128,
    k_bg: np.complex128,
    omega: float,
    reflected: bool = False,
) -> np.ndarray:
    """
    Exact dipole field at a point on the projection sphere.
    If reflected=True, apply the polarization-resolved planar reflection
    correction to that local ray direction (approximation for the interface).
    """
    E = dipole_field(point, source_pos, eps_bg_c, k_bg, p)

    if reflected:
        k_hat_inc = normalize(point - source_pos)
        E = apply_planar_reflection_to_vector(
            E_vec=E,
            k_hat_inc=k_hat_inc,
            eps_inc=eps_bg_c,
            eps_trn=EPS0,   # air
            omega=omega,
        )
    return E


def mode_list(nmax: int):
    modes = []
    for n in range(1, nmax + 1):
        for m in range(-n, n + 1):
            modes.append((n, m))
    return modes


def projection_radius(radius: float, source_dist: float) -> float:
    """
    Choose a projection sphere radius outside the scatterer and inside the
    source location, so the incident field region is source-free.
    """
    r_proj = 1.15 * radius
    if source_dist <= r_proj:
        r_proj = 0.80 * source_dist
    if r_proj <= 1.01 * radius:
        raise ValueError(
            "Source is too close to the sphere for a clean source-free projection surface."
        )
    return r_proj


def project_incident_field_to_vswf(
    source_pos: np.ndarray,
    sphere_center: np.ndarray,
    p: np.ndarray,
    eps_bg_c: np.complex128,
    k_bg: np.complex128,
    omega: float,
    radius: float,
    nmax: int,
    reflected: bool = False,
):
    """
    Numerically project the exact dipole field onto regular VSWFs around the sphere.

    We fit the field on a spherical surface r = r_proj > radius, using the
    full vector field components in local spherical basis.
    """
    source_dist = np.linalg.norm(source_pos - sphere_center)
    r_proj = projection_radius(radius, source_dist)

    # Quadrature grid on the sphere
    n_theta = max(2 * nmax + 4, 12)
    n_phi = max(4 * nmax + 8, 24)

    mu, w_mu = np.polynomial.legendre.leggauss(n_theta)  # mu = cos(theta)
    phis = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    w_phi = 2.0 * np.pi / n_phi

    samples = []
    for i, mui in enumerate(mu):
        theta = np.arccos(np.clip(mui, -1.0, 1.0))
        for phi in phis:
            r_hat, theta_hat, phi_hat = spherical_unit_vectors(theta, phi)
            point = sphere_center + r_proj * r_hat
            E = incident_field_on_surface(
                point=point,
                source_pos=source_pos,
                p=p,
                eps_bg_c=eps_bg_c,
                k_bg=k_bg,
                omega=omega,
                reflected=reflected,
            )
            w = w_mu[i] * w_phi
            samples.append((point, theta, phi, w, E, r_hat, theta_hat, phi_hat))

    modes = mode_list(nmax)
    n_unknowns = 2 * len(modes)
    n_rows = 3 * len(samples)

    A = np.zeros((n_rows, n_unknowns), dtype=complex)
    b = np.zeros(n_rows, dtype=complex)

    sqrt_w_cache = [np.sqrt(s[3]) for s in samples]

    for idx, (point, theta, phi, w, E, r_hat, theta_hat, phi_hat) in enumerate(samples):
        sw = sqrt_w_cache[idx]

        # Data vector: [Er, Etheta, Ephi]
        b[3 * idx + 0] = sw * np.dot(E, r_hat)
        b[3 * idx + 1] = sw * np.dot(E, theta_hat)
        b[3 * idx + 2] = sw * np.dot(E, phi_hat)

        for j, (n, m) in enumerate(modes):
            M_reg, N_reg = vswf_MN(
                n=n,
                m=m,
                k=k_bg,
                r=r_proj,
                theta=theta,
                phi=phi,
                kind="regular",
            )

            A[3 * idx + 0, 2 * j + 0] = sw * np.dot(M_reg, r_hat)
            A[3 * idx + 1, 2 * j + 0] = sw * np.dot(M_reg, theta_hat)
            A[3 * idx + 2, 2 * j + 0] = sw * np.dot(M_reg, phi_hat)

            A[3 * idx + 0, 2 * j + 1] = sw * np.dot(N_reg, r_hat)
            A[3 * idx + 1, 2 * j + 1] = sw * np.dot(N_reg, theta_hat)
            A[3 * idx + 2, 2 * j + 1] = sw * np.dot(N_reg, phi_hat)

    coeff_vec, *_ = np.linalg.lstsq(A, b, rcond=None)

    coeff_M = {}
    coeff_N = {}
    for j, (n, m) in enumerate(modes):
        coeff_M[(n, m)] = coeff_vec[2 * j + 0]
        coeff_N[(n, m)] = coeff_vec[2 * j + 1]

    return coeff_M, coeff_N, nmax


def evaluate_vswf_field(
    r_obs: np.ndarray,
    sphere_center: np.ndarray,
    k_bg: np.complex128,
    coeff_M: dict,
    coeff_N: dict,
    nmax: int,
    kind: str = "outgoing",
) -> np.ndarray:
    """
    Evaluate a VSWF expansion at an observation point.
    """
    r_vec = r_obs - sphere_center
    r, theta, phi = cartesian_to_spherical(r_vec)

    E = np.zeros(3, dtype=complex)
    for n in range(1, nmax + 1):
        for m in range(-n, n + 1):
            M, N = vswf_MN(n, m, k_bg, r, theta, phi, kind=kind)
            E += coeff_M[(n, m)] * M + coeff_N[(n, m)] * N
    return E


def evaluate_scattered_field_from_source(
    r_obs: np.ndarray,
    r_src: np.ndarray,
    sphere_center: np.ndarray,
    p: np.ndarray,
    eps_bg_c: np.complex128,
    k_bg: np.complex128,
    an: np.ndarray,
    bn: np.ndarray,
    radius: float,
    omega: float,
    reflected_incident: bool = False,
    cache=None,
):
    """
    Compute the exact full-wave scattered field from the sphere for a given source.
    This:
      1) projects the exact dipole field onto regular VSWFs,
      2) multiplies by sphere Mie coefficients a_n, b_n,
      3) evaluates the outgoing VSWFs at the receiver.
    """
    if cache is None:
        cache = {}

    source_key = (
        tuple(np.asarray(r_src, dtype=float).tolist()),
        bool(reflected_incident),
        float(np.real(k_bg)),
        float(np.imag(k_bg)),
        float(radius),
        float(omega),
    )

    if source_key not in cache:
        x = k_bg * radius
        xmag = abs(x)
        nmax = max(3, int(np.ceil(xmag + 4.0 * xmag ** (1.0 / 3.0) + 2.0)))
        coeff_M_inc, coeff_N_inc, nmax = project_incident_field_to_vswf(
            source_pos=r_src,
            sphere_center=sphere_center,
            p=p,
            eps_bg_c=eps_bg_c,
            k_bg=k_bg,
            omega=omega,
            radius=radius,
            nmax=nmax,
            reflected=reflected_incident,
        )
        cache[source_key] = (coeff_M_inc, coeff_N_inc, nmax)

    coeff_M_inc, coeff_N_inc, nmax = cache[source_key]

    E_scat = np.zeros(3, dtype=complex)
    for n in range(1, nmax + 1):
        a_n = an[n]
        b_n = bn[n]
        for m in range(-n, n + 1):
            M_out, N_out = vswf_MN(
                n=n,
                m=m,
                k=k_bg,
                r=np.linalg.norm(r_obs - sphere_center),
                theta=cartesian_to_spherical(r_obs - sphere_center)[1],
                phi=cartesian_to_spherical(r_obs - sphere_center)[2],
                kind="outgoing",
            )
            # Standard sphere mapping: M-type gets b_n, N-type gets a_n
            E_scat += b_n * coeff_M_inc[(n, m)] * M_out + a_n * coeff_N_inc[(n, m)] * N_out

    return E_scat


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
    include_direct=True,
    include_reflect_tx=False,
    include_reflect_rx=False,
    include_direct_line=False,   # optional, off by default to preserve your old usage
):
    """
    Compute the total electric field at receiver locations.

    include_direct:
        Keeps your original "source -> sphere -> receiver" unreflected contribution.

    include_direct_line:
        Optional direct source-to-receiver field from the exact dipole Green's function.
        Off by default so existing usage remains unchanged.

    Reflection branches:
        Polarization-resolved but still ray-based. Exact half-space coupling is a
        separate Sommerfeld/layered-Green-function problem.
    """
    results = []

    for f in freqs:
        omega = 2.0 * np.pi * f

        eps_bg_c = complex_permittivity(eps_bg, sigma_bg, omega)
        eps_sph_c = complex_permittivity(eps_sph, sigma_sph, omega)

        k_bg = robust_wavenumber(eps_bg_c, omega)
        k_sph = robust_wavenumber(eps_sph_c, omega)

        # Exact sphere Mie coefficients in lossy background
        x = k_bg * radius
        xmag = abs(x)
        nmax = max(3, int(np.ceil(xmag + 4.0 * xmag ** (1.0 / 3.0) + 2.0)))
        an, bn = mie_coefficients_sphere(k_bg, k_sph, radius, nmax)

        freq_vals = []
        cache = {}

        for r_r in r_receivers:
            E_total = np.zeros(3, dtype=complex)

            # Optional direct source-to-receiver field
            if include_direct_line:
                E_total += dipole_field(r_r, r_s, eps_bg_c, k_bg, p)

            # Unreflected source -> sphere -> receiver path
            if include_direct:
                E_total += evaluate_scattered_field_from_source(
                    r_obs=r_r,
                    r_src=r_s,
                    sphere_center=r_c,
                    p=p,
                    eps_bg_c=eps_bg_c,
                    k_bg=k_bg,
                    an=an,
                    bn=bn,
                    radius=radius,
                    omega=omega,
                    reflected_incident=False,
                    cache=cache,
                )

            # Reflected TX path (source mirrored across z=0, then reflection-aware incident field)
            if include_reflect_tx:
                r_s_eff = reflect_point(r_s)
                E_total += evaluate_scattered_field_from_source(
                    r_obs=r_r,
                    r_src=r_s_eff,
                    sphere_center=r_c,
                    p=p,
                    eps_bg_c=eps_bg_c,
                    k_bg=k_bg,
                    an=an,
                    bn=bn,
                    radius=radius,
                    omega=omega,
                    reflected_incident=True,
                    cache=cache,
                )

            # Reflected RX path (compute at mirrored receiver and reflect back vectorially)
            if include_reflect_rx:
                r_r_eff = reflect_point(r_r)
                E_img = evaluate_scattered_field_from_source(
                    r_obs=r_r_eff,
                    r_src=r_s,
                    sphere_center=r_c,
                    p=p,
                    eps_bg_c=eps_bg_c,
                    k_bg=k_bg,
                    an=an,
                    bn=bn,
                    radius=radius,
                    omega=omega,
                    reflected_incident=False,
                    cache=cache,
                )
                k_hat_out = normalize(r_r_eff - r_c)
                E_total += apply_planar_reflection_to_vector(
                    E_vec=E_img,
                    k_hat_inc=k_hat_out,
                    eps_inc=eps_bg_c,
                    eps_trn=EPS0,
                    omega=omega,
                )

            freq_vals.append(E_total)

        results.append(freq_vals)

    return np.array(results, dtype=complex)


# -----------------------------------------------------------------------------
# Usage example (same style / same outputs as before)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # frequency
    freq = 100e6  # 100 MHz
    freqs = [freq]

    # Geometry (meters)
    r_s = np.array([0.0, 0.0, 0.025])   # Dipole source 2.5 cm above the surface
    r_c = np.array([0.0, 0.0, -0.7])    # Sphere center

    # Medium properties
    eps_bg = 4.0
    sigma_bg = 1e-4

    # Sphere properties (anomaly)
    eps_sph = 20.0
    sigma_sph = 1e-4

    radius = (3e8 / freq) / 16.0
    p = np.array([1.0, 0.0, 0.0])  # Calibration vector

    # Receivers for broadside
    distances_broadside = np.linspace(0.1, 1.0, 48)
    r_receivers_broadside = [np.array([0.0, d, -0.0001]) for d in distances_broadside]

    # Receivers for endfire
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
        include_reflect_rx=False,
        include_direct_line=False,   # keep old usage unchanged
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
        include_reflect_rx=False,
        include_direct_line=False,
    )

    # E_broadside and E_endfire shape: (1, 48, 3)

    # Save to CSV
    with open("broadside_electric_field_sph.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Distance (m)", "Ex_real", "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag"])
        for i, d in enumerate(distances_broadside):
            row = [
                d,
                E_broadside[0, i, 0].real, E_broadside[0, i, 0].imag,
                E_broadside[0, i, 1].real, E_broadside[0, i, 1].imag,
                E_broadside[0, i, 2].real, E_broadside[0, i, 2].imag,
            ]
            writer.writerow(row)

    with open("endfire_electric_field_sph.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Distance (m)", "Ex_real", "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag"])
        for i, d in enumerate(distances_endfire):
            row = [
                d,
                E_endfire[0, i, 0].real, E_endfire[0, i, 0].imag,
                E_endfire[0, i, 1].real, E_endfire[0, i, 1].imag,
                E_endfire[0, i, 2].real, E_endfire[0, i, 2].imag,
            ]
            writer.writerow(row)

    print("Data saved to broadside_electric_field_sph.csv and endfire_electric_field_sph.csv")

    # Create three subplots, one for each measurement (inline and endfire)

    from matplotlib import pyplot as plt
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
    plt.savefig("sphere_green_radial_survey.png", dpi=300)
    plt.show()
