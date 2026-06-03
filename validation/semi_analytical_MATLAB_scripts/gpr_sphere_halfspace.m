% gpr_sphere_halfspace.m
% =========================================================================
% Ex response of an x-directed HED with a dielectric sphere anomaly
% embedded in a half-space, single frequency.
%
% Three output matrices written to separate CSVs:
%   output_primary   : background half-space field (= exx_radial_halfspace result)
%   output_scattered : scattered field from induced sphere dipole only
%   output_total     : primary + scattered
%
% Sphere is on the source vertical axis (x=y=0), depth z_sph.
% Polarisability: quasi-static Clausius-Mossotti (valid when a << lambda_1).
%
% Requires: hfhorx.m, hfhorx_depth.m
% =========================================================================

clear variables
clc

% =========================================================================
% Parameters
% =========================================================================
z      = 0.025;                    % Source height above surface [m]
r_all  = linspace(0.1, 1.0, 48); % Receiver horizontal offsets [m]
epsr1  = 4;                        % Layer 1 relative permittivity
sigma1 = 1e-4;                     % Layer 1 conductivity [S/m]
c0     = 299792458;                % Free-space velocity [m/s]
mu     = 4e-7 * pi;                % Magnetic permeability [H/m]
fc     = 1e8;                      % Centre frequency of Ricker wavelet [Hz]
freq   = fc;
phi    = pi * [0, 45, 90] / 180;   % Endfire, oblique, broadside [rad]

% Sphere anomaly parameters (matches Python SphereAnomaly definition)
sphere_epsr   = 20;                % Sphere relative permittivity
sphere_sigma  = 1e-4;              % Sphere conductivity [S/m]
sphere_mur    = 1.0;               % Sphere relative permeability (no magnetic contrast)
sphere_sigmam = 0.0;               % Sphere magnetic conductivity (unused)
sphere_z      = 0.7;              % Sphere centre z-coordinate [m] (negative = below surface)
sphere_radius = c0 / (freq * 16);  % lambda_air / 16 at fc [m]

z_sph = abs(sphere_z);             % Positive depth used in kernels [m]

fprintf('=== Sphere anomaly — half-space model ===\n');
fprintf('Sphere: a = %.4f m,  z = %.3f m,  epsr = %g,  sigma = %g S/m\n', ...
        sphere_radius, sphere_z, sphere_epsr, sphere_sigma);

% =========================================================================
% EM parameters
% =========================================================================
smu  = 2i * pi * freq * mu;
eta0 = 2i * pi * freq / (c0^2 * mu);
eta1 = eta0 * epsr1  + sigma1;
eta_sph = eta0 * sphere_epsr + sphere_sigma;

% =========================================================================
% Sphere polarisability — Clausius-Mossotti, quasi-static
%
%   alpha = 4*pi*eps1 * a^3 * (eps_sph - eps1) / (eps_sph + 2*eps1)
%
% Complex permittivity from the admittivity eta = sigma + j*omega*eps:
%   eps_c = eta / (j*omega) = epsr*eps0 - j*sigma/omega
%
% Valid when sphere_radius << wavelength_in_layer1.
% =========================================================================
lambda1 = c0 / (freq * sqrt(epsr1));   % approximate (lossless) lambda in layer 1
fprintf('sphere_radius / lambda_1 = %.3f', sphere_radius / lambda1);
if sphere_radius / lambda1 < 0.1
    fprintf('  (quasi-static OK)\n');
else
    fprintf('  [WARNING: quasi-static approximation may be inaccurate]\n');
end

eps1_c  = eta1    / (2i * pi * freq);  % complex permittivity layer 1 [F/m]
eps_sph = eta_sph / (2i * pi * freq);  % complex permittivity sphere   [F/m]
alpha   = 4*pi * eps1_c * sphere_radius^3 * ...
          (eps_sph - eps1_c) ./ (eps_sph + 2*eps1_c);

fprintf('|alpha|   = %.4e  [C m / (V/m)]\n', abs(alpha));

% =========================================================================
% Integration bounds
% =========================================================================
a_bd  = imag(sqrt(eta0 * smu));
b_bd  = imag(sqrt(eta1 * smu));
mxbnd = sqrt((40 / abs(z))^2 + abs(eta0 * smu));

% =========================================================================
% Primary field at sphere centre (r_sph = 0: sphere directly below source)
%
%   J2(kappa * 0) = 0 for all kappa, so vc2_sph = 0 exactly.
%   Only the J0 integral (vc1) contributes; phi_sph = 0.
%
%   exx_sph = (cos(0)*0 - vc1_sph) / (4*pi) = -vc1_sph / (4*pi)
%
% Kernel: hfhorx_depth with a = z (source height), z1 = z_sph.
% See hfhorx_depth.m for derivation.
% =========================================================================
r_sph = 0;

vc1_sph = integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r_sph, z, z_sph, 1), 0, a_bd) + ...
          integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r_sph, z, z_sph, 1), a_bd, b_bd) + ...
          integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r_sph, z, z_sph, 1), b_bd, mxbnd);
vc2_sph = 0;  % exact: J2(0) = 0

Ex_sph = (cos(2 * 0) .* vc2_sph - vc1_sph) / (4*pi);

% Induced dipole moment: p = alpha * E_inc
%   - Only Ex drives px for an inline geometry with sphere on source axis.
%   - Ey = 0 by symmetry (sphere on x-axis of source at y=0).
%   - Ez = 0 by symmetry (source directly overhead, no vertical dipole term).
%   - For an off-axis sphere, Ey and Ez would also induce py, pz components.
%     (%% UNCERTAIN: multi-component case not handled here.)
p_x = alpha * Ex_sph;

fprintf('Primary Ex at sphere:  %.4e + %.4e i  [V/m per A*m]\n', real(Ex_sph), imag(Ex_sph));
fprintf('Induced p_x:           %.4e + %.4e i  [C*m]\n',          real(p_x),    imag(p_x));

% =========================================================================
% Main loop — primary, scattered, total at each receiver offset
% =========================================================================
output_primary   = zeros(length(r_all), 13);
output_scattered = zeros(length(r_all), 13);
output_total     = zeros(length(r_all), 13);

for ir = 1:length(r_all)
    r = r_all(ir);

    % ---- Primary field (identical to exx_radial_halfspace.m) ----------------
    vc1 = integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 1), 0, a_bd) + ...
          integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 1), a_bd, b_bd) + ...
          integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 1), b_bd, mxbnd);
    vc2 = integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 2), 0, a_bd) + ...
          integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 2), a_bd, b_bd) + ...
          integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 2), b_bd, mxbnd);

    exx_p = (cos(2*phi) .* vc2 - vc1) / (4*pi);

    % ---- Scattered field from induced x-dipole at sphere --------------------
    % Kernel: hfhorx_depth with a = 0 (receivers AT the surface; a=0 follows
    % from source-receiver reciprocity — see hfhorx_depth.m header).
    %
    % phi-dependence via cos(2*phi) carries through for all receiver azimuths:
    %   phi=0   (endfire):   exx_sc = p_x*(+vc2_sc - vc1_sc)/(4*pi)
    %   phi=pi/4 (oblique):  exx_sc = p_x*(      0 - vc1_sc)/(4*pi)
    %   phi=pi/2 (broadside):exx_sc = p_x*(-vc2_sc - vc1_sc)/(4*pi)
    % This is valid because the sphere and source share the same x-y position,
    % so the azimuth from sphere to receiver equals the azimuth from source.
    vc1_sc = integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r, 0, z_sph, 1), 0, a_bd) + ...
             integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r, 0, z_sph, 1), a_bd, b_bd) + ...
             integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r, 0, z_sph, 1), b_bd, mxbnd);
    vc2_sc = integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r, 0, z_sph, 2), 0, a_bd) + ...
             integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r, 0, z_sph, 2), a_bd, b_bd) + ...
             integral(@(kappa) hfhorx_depth(kappa, eta0, eta1, smu, r, 0, z_sph, 2), b_bd, mxbnd);

    exx_sc = p_x .* (cos(2*phi) .* vc2_sc - vc1_sc) / (4*pi);

    % ---- Total field --------------------------------------------------------
    exx_t = exx_p + exx_sc;

    % ---- Pack output — column order matches exx_radial_halfspace.m ----------
    output_primary(ir, 1)  = r;
    output_primary(ir, 2)  = abs(exx_p(1));    output_primary(ir, 3)  = angle(exx_p(1));
    output_primary(ir, 4)  = real(exx_p(1));   output_primary(ir, 5)  = imag(exx_p(1));
    output_primary(ir, 6)  = abs(exx_p(3));    output_primary(ir, 7)  = angle(exx_p(3));
    output_primary(ir, 8)  = real(exx_p(3));   output_primary(ir, 9)  = imag(exx_p(3));
    output_primary(ir, 10) = abs(exx_p(2));    output_primary(ir, 11) = angle(exx_p(2));
    output_primary(ir, 12) = real(exx_p(2));   output_primary(ir, 13) = imag(exx_p(2));

    output_scattered(ir, 1)  = r;
    output_scattered(ir, 2)  = abs(exx_sc(1));  output_scattered(ir, 3)  = angle(exx_sc(1));
    output_scattered(ir, 4)  = real(exx_sc(1)); output_scattered(ir, 5)  = imag(exx_sc(1));
    output_scattered(ir, 6)  = abs(exx_sc(3));  output_scattered(ir, 7)  = angle(exx_sc(3));
    output_scattered(ir, 8)  = real(exx_sc(3)); output_scattered(ir, 9)  = imag(exx_sc(3));
    output_scattered(ir, 10) = abs(exx_sc(2));  output_scattered(ir, 11) = angle(exx_sc(2));
    output_scattered(ir, 12) = real(exx_sc(2)); output_scattered(ir, 13) = imag(exx_sc(2));

    output_total(ir, 1)  = r;
    output_total(ir, 2)  = abs(exx_t(1));   output_total(ir, 3)  = angle(exx_t(1));
    output_total(ir, 4)  = real(exx_t(1));  output_total(ir, 5)  = imag(exx_t(1));
    output_total(ir, 6)  = abs(exx_t(3));   output_total(ir, 7)  = angle(exx_t(3));
    output_total(ir, 8)  = real(exx_t(3));  output_total(ir, 9)  = imag(exx_t(3));
    output_total(ir, 10) = abs(exx_t(2));   output_total(ir, 11) = angle(exx_t(2));
    output_total(ir, 12) = real(exx_t(2));  output_total(ir, 13) = imag(exx_t(2));
end

% =========================================================================
% Write CSVs
% =========================================================================
header = {'Offset', ...
    'Endfire_Amplitude','Endfire_Phase','Endfire_Real','Endfire_Imag', ...
    'Broadside_Amplitude','Broadside_Phase','Broadside_Real','Broadside_Imag', ...
    'Oblique_Amplitude','Oblique_Phase','Oblique_Real','Oblique_Imag'};

csv_p = sprintf('Exx_sphere2_primary_%g_%.0fMHz.csv',   epsr1, freq/1e6);
csv_s = sprintf('Exx_sphere2_scattered_%g_%.0fMHz.csv', epsr1, freq/1e6);
csv_t = sprintf('Exx_sphere2_total_%g_%.0fMHz.csv',     epsr1, freq/1e6);

mats  = {output_primary, output_scattered, output_total};
fnames = {csv_p, csv_s, csv_t};
for ii = 1:3
    writecell(header, fnames{ii});
    writematrix(mats{ii}, fnames{ii}, 'Delimiter', ',', 'WriteMode', 'append');
end

% =========================================================================
% Figure 1 — Primary vs Total
% =========================================================================
figure;
subplot(2, 1, 1);
semilogy(r_all, output_primary(:, 6),  'b-',  'LineWidth', 1.5); hold on;
semilogy(r_all, output_total(:, 6),    'b--', 'LineWidth', 1.5);
semilogy(r_all, output_primary(:, 10), 'r-',  'LineWidth', 1.5);
semilogy(r_all, output_total(:, 10),   'r--', 'LineWidth', 1.5);
semilogy(r_all, output_primary(:, 2),  'g-',  'LineWidth', 1.5);
semilogy(r_all, output_total(:, 2),    'g--', 'LineWidth', 1.5);
hold off; grid on;
xlabel('Radial Offset (m)'); ylabel('|Exx| (V/m)');
title(sprintf('|Exx| primary vs total — sphere epsr=%g, z=%.2f m, a=%.4f m', ...
              sphere_epsr, sphere_z, sphere_radius));
legend('Broadside primary','Broadside total','45° primary','45° total', ...
       'Endfire primary','Endfire total','Location','Best');

subplot(2, 1, 2);
plot(r_all, output_primary(:, 7),  'b-',  'LineWidth', 1.5); hold on;
plot(r_all, output_total(:, 7),    'b--', 'LineWidth', 1.5);
plot(r_all, output_primary(:, 11), 'r-',  'LineWidth', 1.5);
plot(r_all, output_total(:, 11),   'r--', 'LineWidth', 1.5);
plot(r_all, output_primary(:, 3),  'g-',  'LineWidth', 1.5);
plot(r_all, output_total(:, 3),    'g--', 'LineWidth', 1.5);
hold off; grid on;
xlabel('Radial Offset (m)'); ylabel('Phase (radians)');
title('Phase primary vs total');
legend('Broadside primary','Broadside total','45° primary','45° total', ...
       'Endfire primary','Endfire total','Location','Best');

% =========================================================================
% Figure 2 — Scattered field and relative anomaly contrast
% =========================================================================
figure;
subplot(2, 1, 1);
semilogy(r_all, output_scattered(:, 6),  'b-',  'LineWidth', 1.5); hold on;
semilogy(r_all, output_scattered(:, 10), 'r--', 'LineWidth', 1.5);
semilogy(r_all, output_scattered(:, 2),  'g-.',  'LineWidth', 1.5);
hold off; grid on;
xlabel('Radial Offset (m)'); ylabel('|Exx| (V/m)');
title('Scattered field magnitude');
legend('Broadside','45°','Endfire','Location','Best');

subplot(2, 1, 2);
contrast_bs = (output_total(:,6)  - output_primary(:,6))  ./ output_primary(:,6)  * 100;
contrast_ef = (output_total(:,2)  - output_primary(:,2))  ./ output_primary(:,2)  * 100;
plot(r_all, contrast_bs, 'b-',  'LineWidth', 1.5); hold on;
plot(r_all, contrast_ef, 'g-.', 'LineWidth', 1.5);
hold off; grid on;
xlabel('Radial Offset (m)'); ylabel('(|E_{total}|-|E_{primary}|)/|E_{primary}| (%)');
title('Relative anomaly contrast');
legend('Broadside','Endfire','Location','Best');
