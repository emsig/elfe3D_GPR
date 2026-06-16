% Simplified routine to compute the x-component of the electric field for
% an x-directed electric dipole source, assuming a single frequency.

clear variables
clc

% Parameters
z = 0.025; % Height of source above surface
r_all = linspace(0.1, 1, 48); % Radial offsets along surface
epsr1 = 4; epsr2 = 9; % Relative permittivity
sigma1 = 1e-4; sigma2 = 1e-3; % Conductivity
d = 1; % Thickness of first layer
c0 = 299792458; % Free space velocity
mu = 4e-7 * pi; % Magnetic permeability
fc = 1e8; % Center frequency of Ricker wavelet
freq = fc; % Single frequency
dt = 1 / (8 * freq); % Time step
phi = pi * [0, 45, 90] / 180; % Angles in radians

% Zeta and eta values
smu = 2i * pi * freq * mu;
eta0 = 2i * pi * freq / (c0^2 * mu);
eta1 = eta0 * epsr1 + sigma1;
eta2 = eta0 * epsr2 + sigma2;

% Preallocate results
output = zeros(length(r_all), 1 + 4 * length(phi)); % Columns: r, Real, Imag, Abs, Phase for each angle

% Loop over radial offsets
for ir = 1:length(r_all)
  r = r_all(ir); % Current radial offset

  % Bounds for integrals
  a = imag(sqrt(eta0 * smu));
  b = imag(sqrt(eta1 * smu));
  mxbnd = sqrt((40 / abs(z))^2 + abs(eta0 * smu));

  % Compute integrals for J0 and J2
  vc1 = integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 1), 0, a) + ...
      integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 1), a, b) + ...
      integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 1), b, mxbnd);
  vc1r = integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 1), 0, a) + ...
       integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 1), a, b) + ...
       integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 1), b, mxbnd);
  vc2 = integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 2), 0, a) + ...
      integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 2), a, b) + ...
      integral(@(kappa) hfhorx(kappa, eta0, eta1, smu, r, z, 2), b, mxbnd);
  vc2r = integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 2), 0, a) + ...
       integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 2), a, b) + ...
       integral(@(kappa) hfhorxr(kappa, eta0, eta1, eta2, smu, r, z, d, 2), b, mxbnd);

  % Compute electric field for the three angles
  exx = (cos(2 * phi).*(vc2 + vc2r) - (vc1 + vc1r)) / (4 * pi);

  % Save results
  output(ir, 1)  = r;                % Offset
  % Endfire (phi=0, index 1)
  output(ir, 2)  = abs(exx(1));      % Endfire_Amplitude
  output(ir, 3)  = angle(exx(1)); % Endfire_Phase (radians)
  output(ir, 4)  = real(exx(1));     % Endfire_Real
  output(ir, 5)  = imag(exx(1));     % Endfire_Imag
  % Broadside (phi=90, index 3)
  output(ir, 6)  = abs(exx(3));      % Broadside_Amplitude
  output(ir, 7)  = angle(exx(3)); % Broadside_Phase (radians)
  output(ir, 8)  = real(exx(3));     % Broadside_Real
  output(ir, 9)  = imag(exx(3));     % Broadside_Imag
  % Oblique (phi=45, index 2)
  output(ir,10)  = abs(exx(2));      % Oblique_Amplitude
  output(ir,11)  = angle(exx(2)); % Oblique_Phase (radians)
  output(ir,12)  = real(exx(2));     % Oblique_Real
  output(ir,13)  = imag(exx(2));     % Oblique_Imag
end

% Write results to a CSV file with new header and order
csv_filename = sprintf('Exx_single_freq_4_9_%.0fMHz_NR.csv', (freq/1e6));
header = {'Offset', ...
    'Endfire_Amplitude','Endfire_Phase','Endfire_Real','Endfire_Imag', ...
    'Broadside_Amplitude','Broadside_Phase','Broadside_Real','Broadside_Imag', ...
    'Oblique_Amplitude','Oblique_Phase','Oblique_Real','Oblique_Imag'};
writecell(header, csv_filename); % Write the header
writematrix(output, csv_filename, 'Delimiter', ',', 'WriteMode', 'append'); % Append the data

% Plot the absolute value and phase
figure;

% Absolute value plot
subplot(2, 1, 1);
semilogy(r_all, output(:, 6), 'b-', 'LineWidth', 1.5); % Broadside Amplitude
hold on;
semilogy(r_all, output(:, 10), 'r--', 'LineWidth', 1.5); % Oblique (45 deg) Amplitude
semilogy(r_all, output(:, 2), 'g-.', 'LineWidth', 1.5); % Endfire Amplitude
hold off;
grid on;
xlabel('Radial Offset (m)');
ylabel('|Exx| (V/m)');
title('Magnitude of Exx at Single Frequency');
legend('Broadside', '45°', 'Endfire', 'Location', 'Best');

% Phase plot
subplot(2, 1, 2);
plot(r_all, output(:, 7), 'b-', 'LineWidth', 1.5); % Broadside Phase (rad)
hold on;
plot(r_all, output(:, 11), 'r--', 'LineWidth', 1.5); % Oblique (45 deg) Phase (rad)
plot(r_all, output(:, 3), 'g-.', 'LineWidth', 1.5); % Endfire Phase (rad)
hold off;
grid on;
xlabel('Radial Offset (m)');
ylabel('Phase (radians)');
title('Phase of Exx at Single Frequency');
legend('Broadside', '45°', 'Endfire', 'Location', 'Best');