% gpr_sphere_two_layers.m
% =========================================================================
% Ex response of an x-directed HED with a dielectric sphere anomaly
% embedded in layer 1 of a two-layer earth, single frequency.
%
% Sphere is always in layer 1 (|sphere_z| < d).
%
% %% UNCERTAIN — depth kernels for two-layer model:
%   hfhorxr_depth uses hfhorxr * exp(-grg1*z1). This is approximate when
%   z_sph is comparable to d, because the layer-2 reflection generates an
%   upgoing wave in layer 1 that is NOT captured by simply multiplying the
%   surface reflection kernel by exp(-grg1*z1).
%   It is exact in the limit z_sph << d. Validated path: half-space only.
%
% Requires: hfhorx.m, hfhorxr.m, hfhorx_depth.m, hfhorxr_depth.m
% =========================================================================

clear variables
clc

% =========================================================================
% Parameters
% =========================================================================
z      = 0.025;                    % Source height above surface [m]
r_all  = linspace(0.1, 3.0, 100); % Receiver horizontal offsets [m]
epsr1  = 4;   epsr2  = 9;         % Layer 1, 2 relative permittivity
sigma1 = 1e-4; sigma2 = 1e-3;     % Layer 1, 2 conductivity [S/m]
d      = 1;                        % Thickness of layer 1 [m]
c0     = 299792458;
mu     = 4e-7 * pi;
fc     = 1e8;
freq   = fc;
phi    = pi * [0, 45, 90] / 180;  % Endfire, oblique, broadside [rad]

% Sphere anomaly parameters
sphere_epsr   = 20;
sphere_sigma  = 1e-4;
sphere_mur    = 1.0;
sphere_sigmam = 0.0;
sphere_z      = -0.7;             % Must satisfy |sphere_z| < d
sphere_radius = c0 / (freq * 16);

z_sph = abs(sphere_z);

if z_sph >= d
    error('Sphere depth (%.3f m) must be less than layer thickness d (%.3f m).', z_sph, d);
end

fprintf('=== Sphere anomaly — two-layer model ===\n');
fprintf('Layer 1: epsr=%g, sigma=%g  |  Layer 2: epsr=%g, sigma=%g  |  d=%.2f m\n', ...
        epsr1, sigma1, epsr2, sigma2, d);
fprintf('Sphere: a=%.4f m,  z=%.3f m,  epsr=%g,  z_sph/d=%.2f\n', ...
        sphere_radius, sphere_z, sphere_epsr, z_sph/d);
if z_sph / d > 0.5
    fprintf('[NOTE: z_sph/d = %.2f — depth kernel approximation less accurate]\n', z_sph/d);
end

% =========================================================================
% EM parameters
% =========================================================================
smu  = 2i * pi * freq * mu;
eta0 = 2i * pi * freq / (c0^2 * mu);
eta1 = eta0 * epsr1 + sigma1;
eta2 = eta0 * epsr2 + sigma2;
eta_sph = eta0 * sphere_epsr + sphere_sigma;

% =========================================================================
% Sphere polarisability
% =========================================================================
lambda1 = c0 / (freq * sqrt(epsr1));
fprintf('sphere_radius / lambda_1 = %.3f', sphere_radius / lambda1);
if sphere_radius / lambda1 < 0.1
    fprintf('  (quasi-static OK)\n');
else
    fprintf('  [WARNING: quasi-static approximation may be inaccurate]\n');
end

eps1_c  = eta1    / (2i * pi * freq);
eps_sph = eta_sph / (2i * pi * freq);
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
% Primary field at sphere centre (two-layer: half-space + layer-2 reflection)
%   r_sph = 0: sphere on source vertical axis → vc2 = 0 exactly.
% =========================================================================
r_sph = 0;

vc1_sph = integral(@(kappa) hfhorx_depth(kappa,  eta0, eta1,       smu, r_sph, z, z_sph, 1), 0, a_bd) + ...
          integral(@(kappa) hfhorx_depth(kappa,  eta0, eta1,       smu, r_sph, z, z_sph, 1), a_bd, b_bd) + ...
          integral(@(kappa) hfhorx_depth(kappa,  eta0, eta1,       smu, r_sph, z, z_sph, 1), b_bd, mxbnd);
vc1r_sph= integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r_sph, z, d, z_sph, 1), 0, a_bd) + ...
          integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r_sph, z, d, z_sph, 1), a_bd, b_bd) + ...
          integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r_sph, z, d, z_sph, 1), b_bd, mxbnd);
vc2_sph  = 0;  % J2(0) = 0 exactly
vc2r_sph = 0;

Ex_sph = (cos(2*0) .* (vc2_sph + vc2r_sph) - (vc1_sph + vc1r_sph)) / (4*pi);
p_x    = alpha * Ex_sph;

fprintf('Primary Ex at sphere:  %.4e + %.4e i  [V/m per A*m]\n', real(Ex_sph), imag(Ex_sph));
fprintf('Induced p_x:           %.4e + %.4e i  [C*m]\n',          real(p_x),    imag(p_x));

% =========================================================================
% Main loop
% =========================================================================
output_primary   = zeros(length(r_all), 13);
output_scattered = zeros(length(r_all), 13);
output_total     = zeros(length(r_all), 13);

for ir = 1:length(r_all)
    r = r_all(ir);

    % ---- Primary field (identical to exx_radial_two_layers.m) --------------
    vc1 = integral(@(kappa) hfhorx( kappa, eta0, eta1,       smu, r, z, 1), 0, a_bd) + ...
          integral(@(kappa) hfhorx( kappa, eta0, eta1,       smu, r, z, 1), a_bd, b_bd) + ...
          integral(@(kappa) hfhorx( kappa, eta0, eta1,       smu, r, z, 1), b_bd, mxbnd);
    vc1r= integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 1), 0, a_bd) + ...
          integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 1), a_bd, b_bd) + ...
          integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 1), b_bd, mxbnd);
    vc2 = integral(@(kappa) hfhorx( kappa, eta0, eta1,       smu, r, z, 2), 0, a_bd) + ...
          integral(@(kappa) hfhorx( kappa, eta0, eta1,       smu, r, z, 2), a_bd, b_bd) + ...
          integral(@(kappa) hfhorx( kappa, eta0, eta1,       smu, r, z, 2), b_bd, mxbnd);
    vc2r= integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 2), 0, a_bd) + ...
          integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 2), a_bd, b_bd) + ...
          integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 2), b_bd, mxbnd);

    exx_p = (cos(2*phi) .* (vc2 + vc2r) - (vc1 + vc1r)) / (4*pi);

    % ---- Scattered field (two-layer: half-space + reflection correction) ----
    vc1_sc = integral(@(kappa) hfhorx_depth( kappa, eta0, eta1,       smu, r, 0, z_sph, 1), 0, a_bd) + ...
             integral(@(kappa) hfhorx_depth( kappa, eta0, eta1,       smu, r, 0, z_sph, 1), a_bd, b_bd) + ...
             integral(@(kappa) hfhorx_depth( kappa, eta0, eta1,       smu, r, 0, z_sph, 1), b_bd, mxbnd);
    vc1r_sc= integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, 0, d, z_sph, 1), 0, a_bd) + ...
             integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, 0, d, z_sph, 1), a_bd, b_bd) + ...
             integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, 0, d, z_sph, 1), b_bd, mxbnd);
    vc2_sc = integral(@(kappa) hfhorx_depth( kappa, eta0, eta1,       smu, r, 0, z_sph, 2), 0, a_bd) + ...
             integral(@(kappa) hfhorx_depth( kappa, eta0, eta1,       smu, r, 0, z_sph, 2), a_bd, b_bd) + ...
             integral(@(kappa) hfhorx_depth( kappa, eta0, eta1,       smu, r, 0, z_sph, 2), b_bd, mxbnd);
    vc2r_sc= integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, 0, d, z_sph, 2), 0, a_bd) + ...
             integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, 0, d, z_sph, 2), a_bd, b_bd) + ...
             integral(@(kappa) hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, 0, d, z_sph, 2), b_bd, mxbnd);

    exx_sc = p_x .* (cos(2*phi) .* (vc2_sc + vc2r_sc) - (vc1_sc + vc1r_sc)) / (4*pi);
    exx_t  = exx_p + exx_sc;

    % ---- Pack output --------------------------------------------------------
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

csv_p = sprintf('Exx_sphere_primary_%g_%g_%.0fMHz.csv',   epsr1, epsr2, freq/1e6);
csv_s = sprintf('Exx_sphere_scattered_%g_%g_%.0fMHz.csv', epsr1, epsr2, freq/1e6);
csv_t = sprintf('Exx_sphere_total_%g_%g_%.0fMHz.csv',     epsr1, epsr2, freq/1e6);

mats  = {output_primary, output_scattered, output_total};
fnames = {csv_p, csv_s, csv_t};
for ii = 1:3
    writecell(header, fnames{ii});
    writematrix(mats{ii}, fnames{ii}, 'Delimiter', ',', 'WriteMode', 'append');
end

% =========================================================================
% Figures
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
title(sprintf('|Exx| primary vs total — epsr1=%g, epsr2=%g, sphere epsr=%g, z=%.2f m', ...
              epsr1, epsr2, sphere_epsr, sphere_z));
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

figure;
subplot(2, 1, 1);
semilogy(r_all, output_scattered(:, 6),  'b-',  'LineWidth', 1.5); hold on;
semilogy(r_all, output_scattered(:, 10), 'r--', 'LineWidth', 1.5);
semilogy(r_all, output_scattered(:, 2),  'g-.', 'LineWidth', 1.5);
hold off; grid on;
xlabel('Radial Offset (m)'); ylabel('|Exx| (V/m)');
title('Scattered field magnitude'); legend('Broadside','45°','Endfire','Location','Best');

subplot(2, 1, 2);
contrast_bs = (output_total(:,6)  - output_primary(:,6))  ./ output_primary(:,6)  * 100;
contrast_ef = (output_total(:,2)  - output_primary(:,2))  ./ output_primary(:,2)  * 100;
plot(r_all, contrast_bs, 'b-',  'LineWidth', 1.5); hold on;
plot(r_all, contrast_ef, 'g-.', 'LineWidth', 1.5); hold off; grid on;
xlabel('Radial Offset (m)');
ylabel('(|E_{total}|-|E_{primary}|)/|E_{primary}| (%)');
title('Relative anomaly contrast'); legend('Broadside','Endfire','Location','Best');
