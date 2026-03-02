clear variables
clc
% plan view just below source level
z=0.01;
% compute the fields in the horizontal plane in a 3x3 m2 grid
xmin=0.2;
xmax=4.8;
nx=256;
xx=linspace(xmin,xmax,nx);
yy=linspace(xmin,xmax,nx);
zl = zeros(1, nx);

% Broadside view
x_bs = zl;
y_bs = yy;

% Endfire view
x_ef = xx;
y_ef = zl;

% Oblique view
x_ob = xx/sqrt(2);
y_ob = yy/sqrt(2);

% Common parameters
c0 = 299792458;
mu = 4e-7*pi;
freq = 1e8;
smu = 2i*pi*freq*mu;
eta0 = 2i*pi*freq/(c0^2*mu);
grg0 = sqrt(smu.*eta0);

% Broadside data set
r_bs = sqrt(x_bs.^2 + y_bs.^2);
R_bs = sqrt(r_bs.^2 + z^2);
exx_bs = ((3*(x_bs./R_bs).^2-1).*(1 + grg0*R_bs) ...
         + ((x_bs./R_bs).^2-1).*(grg0*R_bs).^2).*exp(-grg0*R_bs)./(4*pi*eta0*R_bs.^3);
eyx_bs = x_bs.*y_bs.*(3 + 3*grg0*R_bs + (grg0*R_bs).^2).*exp(-grg0*R_bs)./(4*pi*eta0*R_bs.^5);
ezx_bs = x_bs.*z.*(3 + 3*grg0*R_bs + (grg0*R_bs).^2).*exp(-grg0*R_bs)./(4*pi*eta0*R_bs.^5);

% Endfire data set
r_ef = sqrt(x_ef.^2 + y_ef.^2);
R_ef = sqrt(r_ef.^2 + z^2);
exx_ef = ((3*(x_ef./R_ef).^2-1).*(1 + grg0*R_ef) ...
         + ((x_ef./R_ef).^2-1).*(grg0*R_ef).^2).*exp(-grg0*R_ef)./(4*pi*eta0*R_ef.^3);
eyx_ef = x_ef.*y_ef.*(3 + 3*grg0*R_ef + (grg0*R_ef).^2).*exp(-grg0*R_ef)./(4*pi*eta0*R_ef.^5);
ezx_ef = x_ef.*z.*(3 + 3*grg0*R_ef + (grg0*R_ef).^2).*exp(-grg0*R_ef)./(4*pi*eta0*R_ef.^5);

% Oblique data set
r_ob = sqrt(x_ob.^2 + y_ob.^2);
R_ob = sqrt(r_ob.^2 + z^2);
exx_ob = ((3*(x_ob./R_ob).^2-1).*(1 + grg0*R_ob) ...
         + ((x_ob./R_ob).^2-1).*(grg0*R_ob).^2).*exp(-grg0*R_ob)./(4*pi*eta0*R_ob.^3);
eyx_ob = x_ob.*y_ob.*(3 + 3*grg0*R_ob + (grg0*R_ob).^2).*exp(-grg0*R_ob)./(4*pi*eta0*R_ob.^5);
ezx_ob = x_ob.*z.*(3 + 3*grg0*R_ob + (grg0*R_ob).^2).*exp(-grg0*R_ob)./(4*pi*eta0*R_ob.^5);

% Ref Solution
rr=xx;
x=rr;
R=sqrt(rr.^2+z^2);
exxe=((3*(x./R).^2-1).*(1 + grg0*R) ...
         + ((x./R).^2-1).*(grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^3);
% broadside
x=0;
exxb=((3*(x./R).^2-1).*(1 + grg0*R) ...
         + ((x./R).^2-1).*(grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^3);
% diagonal
x=rr/sqrt(2);
exxd=((3*(x./R).^2-1).*(1 + grg0*R) ...
         + ((x./R).^2-1).*(grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^3);

% Output
% --- Exx ---
output_exx = zeros(nx, 13);
output_exx(:,1)  = r_ef(:); % Offset (use Endfire offset)
output_exx(:,2)  = abs(exx_ef(:));   % Endfire_Amplitude
output_exx(:,3)  = angle(exx_ef(:)); % Endfire_Phase
output_exx(:,4)  = real(exx_ef(:));  % Endfire_Real
output_exx(:,5)  = imag(exx_ef(:));  % Endfire_Imag
output_exx(:,6)  = abs(exx_bs(:));   % Broadside_Amplitude
output_exx(:,7)  = angle(exx_bs(:)); % Broadside_Phase
output_exx(:,8)  = real(exx_bs(:));  % Broadside_Real
output_exx(:,9)  = imag(exx_bs(:));  % Broadside_Imag
output_exx(:,10) = abs(exx_ob(:));   % Oblique_Amplitude
output_exx(:,11) = angle(exx_ob(:)); % Oblique_Phase
output_exx(:,12) = real(exx_ob(:));  % Oblique_Real
output_exx(:,13) = imag(exx_ob(:));  % Oblique_Imag

header_exx = {'Offset', ...
    'Endfire_Amplitude','Endfire_Phase','Endfire_Real','Endfire_Imag', ...
    'Broadside_Amplitude','Broadside_Phase','Broadside_Real','Broadside_Imag', ...
    'Oblique_Amplitude','Oblique_Phase','Oblique_Real','Oblique_Imag'};

csv_filename_exx = sprintf('Exx_single_freq_air_%.0fMHz.csv', (freq/1e6));
writecell(header_exx, csv_filename_exx);
writematrix(output_exx, csv_filename_exx, 'Delimiter', ',', 'WriteMode', 'append');

% --- Eyx ---
output_eyx = zeros(nx, 13);
output_eyx(:,1)  = r_ef(:);
output_eyx(:,2)  = abs(eyx_ef(:));
output_eyx(:,3)  = angle(eyx_ef(:));
output_eyx(:,4)  = real(eyx_ef(:));
output_eyx(:,5)  = imag(eyx_ef(:));
output_eyx(:,6)  = abs(eyx_bs(:));
output_eyx(:,7)  = angle(eyx_bs(:));
output_eyx(:,8)  = real(eyx_bs(:));
output_eyx(:,9)  = imag(eyx_bs(:));
output_eyx(:,10) = abs(eyx_ob(:));
output_eyx(:,11) = angle(eyx_ob(:));
output_eyx(:,12) = real(eyx_ob(:));
output_eyx(:,13) = imag(eyx_ob(:));

header_eyx = strrep(header_exx, 'Exx', 'Eyx');
csv_filename_eyx = sprintf('Eyx_single_freq_air_%.0fMHz.csv', (freq/1e6));
writecell(header_eyx, csv_filename_eyx);
writematrix(output_eyx, csv_filename_eyx, 'Delimiter', ',', 'WriteMode', 'append');

% --- Ezx ---
output_ezx = zeros(nx, 13);
output_ezx(:,1)  = r_ef(:);
output_ezx(:,2)  = abs(ezx_ef(:));
output_ezx(:,3)  = angle(ezx_ef(:));
output_ezx(:,4)  = real(ezx_ef(:));
output_ezx(:,5)  = imag(ezx_ef(:));
output_ezx(:,6)  = abs(ezx_bs(:));
output_ezx(:,7)  = angle(ezx_bs(:));
output_ezx(:,8)  = real(ezx_bs(:));
output_ezx(:,9)  = imag(ezx_bs(:));
output_ezx(:,10) = abs(ezx_ob(:));
output_ezx(:,11) = angle(ezx_ob(:));
output_ezx(:,12) = real(ezx_ob(:));
output_ezx(:,13) = imag(ezx_ob(:));

header_ezx = strrep(header_exx, 'Exx', 'Ezx');
csv_filename_ezx = sprintf('Ezx_single_freq_air_%.0fMHz.csv', (freq/1e6));
writecell(header_ezx, csv_filename_ezx);
writematrix(output_ezx, csv_filename_ezx, 'Delimiter', ',', 'WriteMode', 'append');

% Plot the absolute value and phase
figure;

% Absolute value plot
subplot(2, 1, 1);
semilogy(r_ef, output_exx(:, 6), 'b-', 'LineWidth', 1.5); % Broadside Amplitude
hold on;
semilogy(r_ef, output_exx(:, 10), 'r--', 'LineWidth', 1.5); % Oblique (45 deg) Amplitude
semilogy(r_ef, output_exx(:, 2), 'g-.', 'LineWidth', 1.5); % Endfire Amplitude

% Add reference solutions
semilogy(rr, abs(exxe), 'k-', 'LineWidth', 1.2); % Reference Endfire
semilogy(rr, abs(exxb).*ones(size(rr)), 'c--', 'LineWidth', 1.2); % Reference Broadside (constant)
semilogy(rr, abs(exxd), 'm-.', 'LineWidth', 1.2); % Reference Diagonal

hold off;
grid on;
xlabel('Radial Offset (m)');
ylabel('|Exx| (V/m)');
title('Magnitude of Exx at Single Frequency');
legend('Broadside', '45°', 'Endfire', 'Ref Endfire', 'Ref Broadside', 'Ref 45°', 'Location', 'Best');

% Phase plot
subplot(2, 1, 2);
plot(r_ef, output_exx(:, 7), 'b-', 'LineWidth', 1.5); % Broadside Phase (rad)
hold on;
plot(r_ef, output_exx(:, 11), 'r--', 'LineWidth', 1.5); % Oblique (45 deg) Phase (rad)
plot(r_ef, output_exx(:, 3), 'g-.', 'LineWidth', 1.5); % Endfire Phase (rad)

% Add reference solutions
plot(rr, angle(exxe), 'k-', 'LineWidth', 1.2); % Reference Endfire
plot(rr, angle(exxb).*ones(size(rr)), 'c--', 'LineWidth', 1.2); % Reference Broadside (constant)
plot(rr, angle(exxd), 'm-.', 'LineWidth', 1.2); % Reference Diagonal

hold off;
grid on;
xlabel('Radial Offset (m)');
ylabel('Phase (radians)');
title('Phase of Exx at Single Frequency');
legend('Broadside', '45°', 'Endfire', 'Ref Endfire', 'Ref Broadside', 'Ref 45°', 'Location', 'Best');
