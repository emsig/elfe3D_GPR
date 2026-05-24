% gpr_sphere_anomaly.m
% =========================================================================
% GPR forward model: dielectric sphere anomaly in a layered earth.
%
% Computes the inline Ex response of a horizontal electric dipole (HED)
% over a sphere embedded in the first subsurface layer, using the
% wavenumber-domain (Hankel-transform) kernels hfhorx / hfhorxr.
%
% Pipeline:
%   1. Primary field at surface receivers  (background, no sphere)
%   2. Primary field at sphere centre      (downward continuation)
%   3. Induced dipole on sphere            (quasi-static polarisability)
%   4. Scattered field back to receivers   (buried-dipole upward kernel)
%   5. Total field = primary + scattered
%
% Conventions follow elfe3D / hfhorx / hfhorxr as closely as possible.
% All %% UNCERTAIN tags mark sections where the derivation is approximate
% or where the exact elfe3D normalisation is not confirmed.
%
% Parameters are set to match the Python SphereAnomaly definition:
%   SphereAnomaly(center=(0,0,-0.7), radius=c0/f/16,
%                 properties=(20, 1e-4, 1.0, 0.0))
% =========================================================================

clear; clc; close all;

% =========================================================================
%% 1.  PHYSICAL CONSTANTS
% =========================================================================
c0   = 3e8;               % speed of light [m/s]
mu0  = 4*pi*1e-7;         % free-space permeability [H/m]
eps0 = 1/(mu0*c0^2);      % free-space permittivity [F/m]

% =========================================================================
%% 2.  FREQUENCY & WAVELENGTH
% =========================================================================
f      = 1e8;             % centre frequency [Hz]  — matches ricker_central_f
omega  = 2*pi*f;
lambda_air = c0 / f;      % 3 m in air

% =========================================================================
%% 3.  SPHERE ANOMALY (matches Python SphereAnomaly definition)
% =========================================================================
sphere_center  = [0.0, 0.0, -0.7];   % (x, y, z) [m]
sphere_radius  = lambda_air / 16;    % 0.1875 m
sphere_eps_r   = 20.0;
sphere_sigma   = 1e-4;               % [S/m]
sphere_mu_r    = 1.0;
sphere_sigma_m = 0.0;                % magnetic conductivity [Ω/m], unused here

fprintf('Sphere radius : %.4f m  (lambda_air/16)\n', sphere_radius);

% =========================================================================
%% 4.  LAYER MODEL
%      Model: air (z > 0) | layer 1 (z < 0) [| layer 2 below depth d12]
%      The sphere is always in layer 1.
% =========================================================================
two_layer = false;         % set true to activate layer 2 + hfhorxr kernels

% Air
air_eps_r  = 1.0;
air_sigma  = 1e-16;        % effectively zero

% Layer 1  (contains the sphere)
layer1_eps_r   = 4.0;      % relative permittivity
layer1_sigma   = 1e-3;     % conductivity [S/m]
layer1_mu_r    = 1.0;

% Layer 2  (only used when two_layer = true)
layer2_eps_r   = 10.0;
layer2_sigma   = 1e-2;
layer2_mu_r    = 1.0;
d12            = 1.5;      % depth to layer1/layer2 interface [m]
                           % must be deeper than the sphere centre |z| = 0.7 m

% =========================================================================
%% 5.  SOURCE
%      Horizontal electric dipole, x-directed, unit moment.
%      Height above surface matches survey.py default: z = 0.025 m.
% =========================================================================
src_pos = [0.0, 0.0, 0.025];   % (x, y, z) [m]
src_a   = src_pos(3);           % height 'a' used in hfhorx / hfhorxr

% =========================================================================
%% 6.  RECEIVER ARRAY  (inline, y = 0, at surface z = 0)
% =========================================================================
rx_x = linspace(-3.0, 3.0, 121);   % [m]
rx_y = zeros(size(rx_x));
rx_z = zeros(size(rx_x));

% =========================================================================
%% 7.  DERIVED EM QUANTITIES
% =========================================================================
smu  = 1j * omega * mu0;           % j·ω·μ₀  [used as 'smu' in kernels]

% Complex admittivities: η = σ + jωε  [S/m]  — convention matches hfhorx
eta0 = air_sigma    + 1j*omega * air_eps_r    * eps0;
eta1 = layer1_sigma + 1j*omega * layer1_eps_r * eps0;
if two_layer
    eta2 = layer2_sigma + 1j*omega * layer2_eps_r * eps0;
else
    eta2 = eta1;   % placeholder — not used when two_layer = false
end

% Complex permittivity of layer 1 and sphere (used for polarisability)
%   εc = ε_r·ε₀ − j·σ/ω   [F/m]
eps1_c      = layer1_eps_r * eps0 - 1j * layer1_sigma  / omega;
eps_sphere  = sphere_eps_r * eps0 - 1j * sphere_sigma  / omega;

% Wavenumber in layer 1 (complex, for free-space Green's function fallback)
k1 = sqrt(-smu * eta1);            % k₁ = ω√(μ₁ε₁)
lambda1 = 2*pi / real(k1);

fprintf('Wavelength in layer 1 : %.3f m\n', lambda1);
fprintf('sphere_radius / lambda1 = %.3f  ', sphere_radius/lambda1);
if sphere_radius/lambda1 < 0.1
    fprintf('(quasi-static approx. valid)\n');
else
    fprintf('[WARN: quasi-static approx. may be inaccurate]\n');
end

% =========================================================================
%% 8.  SPHERE POLARISABILITY  (Clausius-Mossotti, quasi-static)
%
%   α = 4π ε₁ a³ · (εs − ε₁)/(εs + 2ε₁)
%
%   Valid when a << λ₁. For a Mie-theory correction use Mie coefficients.
%   Here we use the electric-dipole (lowest-order) term only.
%   Because sphere_sigma_m = 0 and sphere_mu_r = 1, the magnetic
%   polarisability is zero and only the electric dipole is induced.
% =========================================================================
alpha = 4*pi * eps1_c * sphere_radius^3 * ...
        (eps_sphere - eps1_c) ./ (eps_sphere + 2*eps1_c);

fprintf('|alpha| = %.4e  [C·m / (V/m)]\n', abs(alpha));

% =========================================================================
%% 9.  WAVENUMBER INTEGRATION GRID
%
% %% UNCERTAIN:  kappa_max and nk should be tuned so that
%   (a) the Bessel functions have decayed sufficiently by kappa_max,
%   (b) the shortest spatial scale (sphere depth, antenna height) is sampled.
%   A rule of thumb: kappa_max ~ 10/min(src_a, |z_sphere|).
%   Use log-spacing for better sampling near kappa=0.
% =========================================================================
kappa_max = 20 / min(src_a, abs(sphere_center(3)));   % adaptive upper limit
nk        = 4000;
kappa     = logspace(log10(1e-4), log10(kappa_max), nk);   % log-spaced [1/m]

% =========================================================================
%% 10.  PRIMARY FIELD AT SURFACE RECEIVERS
%
%   For an x-directed unit HED at height src_a, the surface electric field is
%   built from two Hankel-transform integrals (J₀ and J₂ kernels):
%
%     Ex = (1/4π) · [ I₀ − I₂·cos(2φ) ]
%     Ey = −(1/4π) · I₂·sin(2φ)
%
%   where φ is the azimuth from source to receiver and I₀, I₂ are
%   integrals of the kernels returned by hfhorx / hfhorxr.
%
% %% UNCERTAIN:  The (1/4π) prefactor follows Slob et al. conventions.
%   Confirm against your elfe3D source.txt normalisation (current moment,
%   segment length) if an absolute-amplitude comparison is needed.
% =========================================================================
nrx        = numel(rx_x);
Ex_primary = zeros(1, nrx);
Ey_primary = zeros(1, nrx);

for irx = 1:nrx
    dx      = rx_x(irx) - src_pos(1);
    dy      = rx_y(irx) - src_pos(2);
    r_horiz = sqrt(dx^2 + dy^2);

    if r_horiz < 1e-9
        Ex_primary(irx) = NaN;   % skip source location
        continue
    end

    phi = atan2(dy, dx);

    if two_layer
        k_J0 = hfhorxr(kappa, eta0, eta1, eta2, smu, r_horiz, src_a, d12, 1);
        k_J2 = hfhorxr(kappa, eta0, eta1, eta2, smu, r_horiz, src_a, d12, 2);
    else
        k_J0 = hfhorx(kappa, eta0, eta1, smu, r_horiz, src_a, 1);
        k_J2 = hfhorx(kappa, eta0, eta1, smu, r_horiz, src_a, 2);
    end

    I0 = trapz(kappa, k_J0);
    I2 = trapz(kappa, k_J2);

    Ex_primary(irx) = (1/(4*pi)) * (I0 - I2 * cos(2*phi));
    Ey_primary(irx) = -(1/(4*pi)) *  I2 * sin(2*phi);
end

% =========================================================================
%% 11.  PRIMARY FIELD AT SPHERE CENTRE  (downward continuation into layer 1)
%
%   The sphere is at depth z_sph = |sphere_center(3)| below the interface.
%   The field inside layer 1 requires downward-continued kernels that
%   include the air→layer1 transmission and propagation exp(-γ₁·z_sph).
%
% %% UNCERTAIN:  The transmission kernels below are derived from the
%   standard Sommerfeld formulation for a HED above a half-space.
%   They have NOT been verified against the elfe3D reference output.
%   See hfhorx_depth() local function for the derivation.
% =========================================================================
z_sph   = abs(sphere_center(3));            % depth [m]
dx_sph  = sphere_center(1) - src_pos(1);
dy_sph  = sphere_center(2) - src_pos(2);
r_sph   = sqrt(dx_sph^2 + dy_sph^2);       % horizontal distance src→sphere
phi_sph = atan2(dy_sph, dx_sph);

if r_sph < 1e-9
    % Sphere directly below source: only J₀ contributes (J₂ → 0)
    phi_sph = 0;
    warning('Sphere is on the source axis; J₂ integral vanishes by symmetry.');
end

if two_layer
    k_J0_sph = hfhorxr_depth(kappa, eta0, eta1, eta2, smu, ...
                              r_sph, src_a, d12, z_sph, 1);
    k_J2_sph = hfhorxr_depth(kappa, eta0, eta1, eta2, smu, ...
                              r_sph, src_a, d12, z_sph, 2);
else
    k_J0_sph = hfhorx_depth(kappa, eta0, eta1, smu, r_sph, src_a, z_sph, 1);
    k_J2_sph = hfhorx_depth(kappa, eta0, eta1, smu, r_sph, src_a, z_sph, 2);
end

I0_sph = trapz(kappa, k_J0_sph);
I2_sph = trapz(kappa, k_J2_sph);

Ex_sph = (1/(4*pi)) * (I0_sph - I2_sph * cos(2*phi_sph));
Ey_sph = -(1/(4*pi)) *  I2_sph * sin(2*phi_sph);
% Ez at sphere: %% UNCERTAIN — requires a separate vertical-component kernel
%   (involves J₁ Bessel function); set to 0 for now.
Ez_sph = 0;

E_at_sphere = [Ex_sph; Ey_sph; Ez_sph];   % [V/m per A·m]

fprintf('\nPrimary field at sphere centre:\n');
fprintf('  Ex = %.4e + %.4e i\n', real(Ex_sph), imag(Ex_sph));
fprintf('  Ey = %.4e + %.4e i\n', real(Ey_sph), imag(Ey_sph));

% =========================================================================
%% 12.  INDUCED DIPOLE MOMENT ON SPHERE
%
%   p = α · E_inc   [C·m]
%
%   The x-directed primary field dominates for an inline geometry.
%   y and z components are included for completeness.
% =========================================================================
p = alpha * E_at_sphere;   % induced dipole moment vector [C·m]

fprintf('\nInduced dipole moment:\n');
fprintf('  |px| = %.4e  Cm\n', abs(p(1)));
fprintf('  |py| = %.4e  Cm\n', abs(p(2)));
fprintf('  |pz| = %.4e  Cm  (from Ez_sph=0, unreliable)\n', abs(p(3)));

% =========================================================================
%% 13.  SCATTERED FIELD FROM INDUCED DIPOLE TO SURFACE RECEIVERS
%
%   The scattered field is computed as the field of an electric dipole p
%   buried at (sphere_center) and observed at the surface receivers.
%
%   Approach A — PREFERRED but %% UNCERTAIN:
%     Use buried-dipole upward-continuation kernels (hfhorx_buried).
%     These are the reciprocal of the downward kernels in Section 11.
%     By source-receiver reciprocity, the kernel for "dipole at depth z_sph,
%     field at surface" equals the downward kernel with src_a ↔ z_sph swapped.
%     Only the x-component of p is handled here (dominant for inline).
%
%   Approach B — fallback, free-space Green's function in layer 1:
%     Used only if buried_kernel = false. Ignores the air/earth interface
%     on the upward path. Valid only if z_sph >> skin_depth.
% =========================================================================
buried_kernel = true;     % toggle: true = reciprocity kernel, false = free-space

Ex_scattered = zeros(1, nrx);

if buried_kernel
    % --- Approach A: reciprocity (swap src_a <-> z_sph in downward kernel) ---
    %
    % %% UNCERTAIN:  The sign of the exponential in hfhorx_depth uses
    %   exp(-γ₀·a)·exp(-γ₁·z_sph). Swapping a and z_sph is exact for
    %   a half-space only when the medium above is homogeneous.
    %   For two_layer, this approach needs further derivation.

    for irx = 1:nrx
        if isnan(Ex_primary(irx)), continue; end

        dx      = rx_x(irx) - sphere_center(1);
        dy      = rx_y(irx) - sphere_center(2);
        r_rx    = sqrt(dx^2 + dy^2);
        phi_rx  = atan2(dy, dx);

        if r_rx < 1e-9, continue; end

        % Buried x-dipole field at surface: use upward kernel
        % By reciprocity: swap src_a (height of transmitter in air) with
        % depth of "virtual transmitter" = z_sph
        k_J0_sc = hfhorx_depth(kappa, eta0, eta1, smu, r_rx, z_sph, src_a, 1);
        k_J2_sc = hfhorx_depth(kappa, eta0, eta1, smu, r_rx, z_sph, src_a, 2);
        %  ^^ Note argument order swap: src_height=z_sph, eval_depth=src_a
        %     src_a acts as the "depth" in air (height above interface → 0
        %     for receiver exactly at surface; use src_a=0 for true surface
        %     receivers). %% UNCERTAIN: check sign convention.

        I0_sc = trapz(kappa, k_J0_sc);
        I2_sc = trapz(kappa, k_J2_sc);

        % Scattered field from px component of induced dipole only
        %   (py contribution requires cross-term; pz requires Ez kernel)
        Ex_scattered(irx) = (1/(4*pi)) * p(1) * ...
                            (I0_sc - I2_sc * cos(2*phi_rx));
    end

else
    % --- Approach B: free-space Green's function in layer 1 ---------------
    %   E = (1/4π) (k₁² G + ∂²G/∂x²) · px / ε₁
    %   where G = exp(-jk₁R) / R,  R = distance from sphere to receiver.
    %
    % %% UNCERTAIN:  Ignores interface reflection on upward path.

    for irx = 1:nrx
        if isnan(Ex_primary(irx)), continue; end

        dx = rx_x(irx) - sphere_center(1);
        dy = rx_y(irx) - sphere_center(2);
        dz = rx_z(irx) - sphere_center(3);   % dz > 0: receiver above sphere
        R  = sqrt(dx^2 + dy^2 + dz^2);
        if R < 1e-9, continue; end

        jkR = 1j * k1 * R;
        G   = exp(-jkR) / (4*pi*R);

        % Dyadic factor G_xx = k₁² G + ∂²G/∂x²
        % ∂²G/∂x² = G·[(3x²/R² − 1)·(1/R² + jk₁/R) − k₁²·x²/R²]
        fac = (3*dx^2/R^2 - 1) * (1/R^2 + jkR/R^2) - k1^2 * dx^2/R^2;
        Gxx = k1^2 * G + G * fac;

        Ex_scattered(irx) = Gxx * p(1) / eps1_c;
    end
end

% =========================================================================
%% 14.  TOTAL FIELD
% =========================================================================
Ex_total = Ex_primary + Ex_scattered;

% =========================================================================
%% 15.  PLOTS
% =========================================================================
valid = ~isnan(Ex_primary);

figure('Name', 'GPR Sphere Anomaly — Surface Ex', 'Position', [100 100 900 700]);

subplot(3,1,1);
plot(rx_x(valid), real(Ex_primary(valid)), 'b',  'LineWidth', 1.5); hold on;
plot(rx_x(valid), imag(Ex_primary(valid)), 'b--','LineWidth', 1.5);
xline(sphere_center(1), 'k:', 'Sphere'); grid on;
title(sprintf('Primary E_x  (f=%g MHz, layer1 \\epsilon_r=%.0f)', f/1e6, layer1_eps_r));
ylabel('E_x  [V/m per A·m]'); legend('Real','Imag');

subplot(3,1,2);
plot(rx_x(valid), real(Ex_scattered(valid)), 'r',  'LineWidth', 1.5); hold on;
plot(rx_x(valid), imag(Ex_scattered(valid)), 'r--','LineWidth', 1.5);
xline(sphere_center(1), 'k:', 'Sphere'); grid on;
title(sprintf('Scattered E_x  (sphere \\epsilon_r=%g, a=%.4f m)', ...
              sphere_eps_r, sphere_radius));
ylabel('E_x  [V/m per A·m]'); legend('Real','Imag');

subplot(3,1,3);
plot(rx_x(valid), real(Ex_total(valid)), 'k',  'LineWidth', 1.5); hold on;
plot(rx_x(valid), imag(Ex_total(valid)), 'k--','LineWidth', 1.5);
xline(sphere_center(1), 'k:', 'Sphere'); grid on;
title('Total E_x  (primary + scattered)');
xlabel('x_{receiver}  [m]'); ylabel('E_x  [V/m per A·m]'); legend('Real','Imag');

% Ratio plot — highlight anomaly contrast
figure('Name', 'Anomaly Contrast');
contrast = abs(Ex_total) ./ abs(Ex_primary) - 1;
plot(rx_x(valid), contrast(valid)*100, 'k', 'LineWidth', 1.5);
xline(sphere_center(1), 'r--', sprintf('z_{sphere}=%.2f m', sphere_center(3)));
xlabel('x_{receiver}  [m]');
ylabel('|E_{total}|/|E_{primary}| − 1  [%]');
title('Relative anomaly contrast'); grid on;

% =========================================================================
%% LOCAL FUNCTIONS  (kernels mirror hfhorx / hfhorxr from user-supplied code)
% =========================================================================

% -------------------------------------------------------------------------
function exyz = hfhorx(kappa, eta0, eta1, smu, r, a, vc)
% HFHORX  Half-space Hankel kernel for HED, field at surface (z=0 in air).
%   One subsurface layer.  Matches user-supplied hfhorx.m exactly.
%
%   NOTE: J2kc uses kappa.^3 in the original hfhorx vs kappa in hfhorxr.
%   This may reflect a different normalisation of the J2 integral.
grg0 = sqrt(smu*eta0 + kappa.^2);
grg1 = sqrt(smu*eta1 + kappa.^2);
J0kc = besselj(0, kappa*r) .* kappa;
J2kc = besselj(2, kappa*r) .* kappa.^3;   % %% UNCERTAIN: kappa^3 vs kappa^1
if vc == 2
    exyz = J2kc .* exp(-grg0*a) ./ (eta1.*grg0 + eta0.*grg1);
else
    exyz = J0kc .* exp(-grg0*a) .* ...
           (2*smu ./ (grg0+grg1) + kappa.^2 ./ (eta1.*grg0 + eta0.*grg1));
end
end

% -------------------------------------------------------------------------
function exyz = hfhorxr(kappa, eta0, eta1, eta2, smu, r, a, d, vc)
% HFHORXR  Two-layer Hankel kernel, field at surface.
%   Matches user-supplied hfhorxr.m exactly.
grg0 = sqrt(smu*eta0 + kappa.^2);
grg1 = sqrt(smu*eta1 + kappa.^2);
grg2 = sqrt(smu*eta2 + kappa.^2);

r0te = (grg0-grg1) ./ (grg0+grg1);
r1te = (grg1-grg2) ./ (grg1+grg2);
r0tm = (eta1.*grg0-eta0.*grg1) ./ (eta1.*grg0+eta0.*grg1);
r1tm = (eta2.*grg1-eta1.*grg2) ./ (eta2.*grg1+eta1.*grg2);

e2   = exp(-2*grg1*d);
glre = (1-r0te.^2) .* r1te .* e2 ./ (1+r0te.*r1te.*e2);
glrm = (1-r0tm.^2) .* r1tm .* e2 ./ (1+r0tm.*r1tm.*e2);

J0kc = besselj(0, kappa*r) .* kappa;
J2kc = besselj(2, kappa*r) .* kappa;    % kappa^1 — matches user's hfhorxr
if vc == 2
    exyz = J2kc .* exp(-grg0*a) .* (grg0.*glrm/eta0 + smu*glre./grg0);
else
    exyz = J0kc .* exp(-grg0*a) .* (grg0.*glrm/eta0 - smu*glre./grg0);
end
end

% -------------------------------------------------------------------------
function exyz = hfhorx_depth(kappa, eta0, eta1, smu, r, a, z1, vc)
% HFHORX_DEPTH  Downward-continued kernel for HED field at depth z1
%               inside layer 1 (half-space model).
%
% Derivation (Sommerfeld half-space):
%   Source in air at height a; observation point at depth z1 > 0 in layer 1.
%   The transmitted field uses transmission coefficients:
%     T_TE = 2 γ₀ / (γ₀ + γ₁)         (TE mode)
%     T_TM = 2 η₁ γ₀ / (η₁γ₀ + η₀γ₁)  (TM mode)
%   Propagation factor: exp(−γ₀ a) · exp(−γ₁ z₁)
%
%   J₀ kernel combines TM and TE contributions:
%     I₀ ~ J₀ · exp(−γ₀a−γ₁z₁) · [smu·T_TE/γ₁ + κ²·T_TM/(η₁·γ₁·2)]
%          %% UNCERTAIN: the η₁ and factor-of-2 in the TM denominator
%          are best guesses matching the surface kernel structure.
%          Cross-check against Chave & Cox (1982) or Ward & Hohmann (1988).
%
%   J₂ kernel (TM only):
%     I₂ ~ J₂ · exp(−γ₀a−γ₁z₁) · T_TM / (η₁·γ₁)
%          %% UNCERTAIN: normalisation factor.
%
%   For the reciprocal path (buried dipole → surface), call with
%   a  = depth of buried dipole  [m]
%   z1 = height of receiver above interface (0 for surface receiver in air)
%   %% UNCERTAIN: the exp(−γ₀·z1) factor then equals 1 when z1=0.

grg0  = sqrt(smu*eta0 + kappa.^2);
grg1  = sqrt(smu*eta1 + kappa.^2);
prop  = exp(-grg0*a) .* exp(-grg1*z1);   % source-to-interface-to-depth

T_TE  = 2*grg0 ./ (grg0 + grg1);
T_TM  = 2*eta1.*grg0 ./ (eta1.*grg0 + eta0.*grg1);

J0kc  = besselj(0, kappa*r) .* kappa;
J2kc  = besselj(2, kappa*r) .* kappa;    % %% UNCERTAIN: kappa vs kappa^3

if vc == 2
    exyz = J2kc .* prop .* T_TM ./ (eta1.*grg1);
else
    exyz = J0kc .* prop .* ...
           (smu .* T_TE ./ grg1  +  kappa.^2 .* T_TM ./ (2*eta1.*grg1));
end
end

% -------------------------------------------------------------------------
function exyz = hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, a, d, z1, vc)
% HFHORXR_DEPTH  Downward-continued kernel for two-layer model.
%                Observation depth z1 in layer 1 (must be < d).
%
% Uses the same transmission coefficients as hfhorxr but adds exp(−γ₁z₁).
% The reflection response glre / glrm is absent on the downward path
% (no upward-going wave in a purely transmitted field).
%
% %% UNCERTAIN:  This does not account for reflections from the layer1/layer2
%   boundary propagating back upward and then downward again. These
%   multiple reflections are neglected. Valid when z1 << d or when
%   layer2 contrast is small.

grg0  = sqrt(smu*eta0 + kappa.^2);
grg1  = sqrt(smu*eta1 + kappa.^2);
grg2  = sqrt(smu*eta2 + kappa.^2);  %#ok<NASGU>  % included for completeness

T_TE  = 2*grg0 ./ (grg0 + grg1);
T_TM  = 2*eta1.*grg0 ./ (eta1.*grg0 + eta0.*grg1);

prop  = exp(-grg0*a) .* exp(-grg1*z1);

J0kc  = besselj(0, kappa*r) .* kappa;
J2kc  = besselj(2, kappa*r) .* kappa;

if vc == 2
    exyz = J2kc .* prop .* T_TM ./ (eta1.*grg1);
else
    exyz = J0kc .* prop .* ...
           (smu .* T_TE ./ grg1  +  kappa.^2 .* T_TM ./ (2*eta1.*grg1));
end
end
