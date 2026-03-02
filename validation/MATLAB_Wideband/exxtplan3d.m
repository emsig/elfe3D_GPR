% routine to compute the x-component of the electric field for an
% x-directed electric dipole source, receiver is assumed to be on the
% ground and the source very close to but just above the ground

% if receiver is not on but above the ground, we need to include the direct
% contributions, but I assume we will not use that here
clear variables
clc
% height of source above surface
z=0.025;
% radial offset along surface
r=0.50213;
% number of frequencies
nf=250;
% relative espilon and sigma values of first layer and half space
epsr1=4;
epsr2=9;
sigma1=1e-3;
sigma2=1e-2;
% thickness of first layer
d=1;
% free space velocity and magnetic permeability
c0=299792458;
mu=4e-7*pi;
% center frequency of Ricker wavelet
fc=1e8;
% frequency axis
freq=(1:nf)*6*fc/nf;
disp(freq);
% zero padding factor
mf=8;
% time step and time axis
dt=1/(mf*max(freq));
time=(0:nf)*dt*1e9;
% ricker wavelet
wav=2*(freq./fc).^2.*exp(-(freq./fc).^2)/(fc*sqrt(pi));
% number of angles on ground surface to compute the electric field
nphi=3;
% angle axis
phi=pi*linspace(0,90,nphi)'/180;
% zeta and eta values
smu=2i*pi*freq*mu;
eta0=2i*pi*freq/(c0^2*mu);
eta1=eta0*epsr1+sigma1;
eta2=eta0*epsr2+sigma2;
% preset arrays to zero that occur in loops
% electric fields on surface and reflection response
exx=zeros(nphi,nf+1);
exx_o=zeros(nphi,nf+1);
exxr=zeros(nphi,nf+1);
% electric sub-fields on surface and reflection response
vc1=zeros(1,nf);
vc2=zeros(1,nf);
vc1r=zeros(1,nf);
vc2r=zeros(1,nf);
% set tolerance for integrations
atol=1e-16;
rtol=1e-8;
% loop over frequency values
for jf=1:nf
% bounds for integrals
    a=imag(sqrt(eta0(jf)*smu(jf)));
    b=imag(sqrt(eta1(jf)*smu(jf)));
    mxbnd=sqrt((40/abs(z))^2+abs(eta0(jf)*smu(jf)));
    vc=1; % J0-integral
    vc1(jf)=integral(@(kappa)hfhorx(kappa,eta0(jf),eta1(jf),smu(jf),r,z,vc)...
                ,0,a,'reltol',rtol);
    vc1(jf)=vc1(jf) + integral(@(kappa)hfhorx(kappa,eta0(jf),eta1(jf),smu(jf),r,z,vc)...
                ,a,b,'reltol',rtol);
    vc1(jf)=vc1(jf) + integral(@(kappa)hfhorx(kappa,eta0(jf),eta1(jf),smu(jf),r,z,vc)...
                ,b,mxbnd,'reltol',rtol);
    vc1r(jf)=integral(@(kappa)hfhorxr(kappa,eta0(jf),eta1(jf),eta2(jf),smu(jf),r,z,d,vc)...
                ,0,a,'reltol',rtol);
    vc1r(jf)=vc1r(jf) + integral(@(kappa)hfhorxr(kappa,eta0(jf),eta1(jf),eta2(jf),smu(jf),r,z,d,vc)...
                ,a,b,'reltol',rtol);
    vc1r(jf)=vc1r(jf) + integral(@(kappa)hfhorxr(kappa,eta0(jf),eta1(jf),eta2(jf),smu(jf),r,z,d,vc)...
                ,b,mxbnd,'reltol',rtol);
    vc=2; %J2-integral
    vc2(jf)=integral(@(kappa)hfhorx(kappa,eta0(jf),eta1(jf),smu(jf),r,z,vc)...
                ,0,a,'reltol',rtol);
    vc2(jf)=vc2(jf) + integral(@(kappa)hfhorx(kappa,eta0(jf),eta1(jf),smu(jf),r,z,vc)...
                ,a,b,'reltol',rtol);
    vc2(jf)=vc2(jf) + integral(@(kappa)hfhorx(kappa,eta0(jf),eta1(jf),smu(jf),r,z,vc)...
                ,b,mxbnd,'reltol',rtol);
    vc2r(jf)=integral(@(kappa)hfhorxr(kappa,eta0(jf),eta1(jf),eta2(jf),smu(jf),r,z,d,vc)...
                ,0,a,'reltol',rtol);
    vc2r(jf)=vc2r(jf) + integral(@(kappa)hfhorxr(kappa,eta0(jf),eta1(jf),eta2(jf),smu(jf),r,z,d,vc)...
                ,a,b,'reltol',rtol);
    vc2r(jf)=vc2r(jf) + integral(@(kappa)hfhorxr(kappa,eta0(jf),eta1(jf),eta2(jf),smu(jf),r,z,d,vc)...
                ,b,mxbnd,'reltol',rtol);
end
% compute electric field from the sub-results for the three angles
% note that the multiplication with the exponential function is to delay
% the time pulse such that the direct waves lie entirely on the positive
% time axis, the dealy is T0=1/fc

exx_o(:,2:nf+1)=-(repmat(cos(2*phi),1,nf).*repmat(vc2+vc2r,nphi,1) -...
  repmat(vc1+vc1r,nphi,1))./(8*pi);

exx(:,2:nf+1)=-(repmat(cos(2*phi),1,nf).*repmat(vc2+vc2r,nphi,1) -...
  repmat(vc1+vc1r,nphi,1)).*repmat(wav.*exp(-2i*pi*freq/fc),nphi,1)./(8*pi);

% Define names of output files
fnamet = ["Exxbs", "ExxJ0", "Exxef"];  % Broadside, 45 degrees, and Endfire

% Plot each result and write to file
for ip = 1:3
    figure('Units', 'normalized', 'Position', [0 0 1 0.7]); % Full width, 40% of screen height

    % Real part
    subplot(1,3,1)
    plot(freq(1:nf), real(exx(ip,1:nf)), 'ko-', 'linewidth', 2)
    xlabel('Frequency')
    ylabel(sprintf('Real(%s) (V/m)', fnamet(ip)))
    set(gca, 'fontsize', 18)

    % Imaginary part
    subplot(1,3,2)
    plot(freq(1:nf), imag(exx(ip,1:nf)), 'ro-', 'linewidth', 2)
    xlabel('Frequency')
    ylabel(sprintf('Imag(%s) (V/m)', fnamet(ip)))
    set(gca, 'fontsize', 18)

    % Magnitude
    subplot(1,3,3)
    loglog(freq(1:nf), abs(exx(ip,1:nf)), 'bo-', 'linewidth', 2)
    xlabel('Frequency')
    ylabel(sprintf('|%s| (V/m)', fnamet(ip)))
    set(gca, 'fontsize', 18)

    % Super title for the figure
    sgtitle(['Reflection response, d = ', num2str(d,2), ' cm, \epsilon_{r;1}=', num2str(epsr1,2), ', \epsilon_{r;2}=', num2str(epsr2,2), ', \phi = ', num2str(phi(ip) * 180 / pi, 2), ' degrees'], 'FontSize', 20, 'FontWeight', 'bold')

    % Save the figure with the corresponding filename
    print(fnamet(ip), '-dpng') 
end

% plot the three fields in one plot for comparison
figure('Units', 'normalized', 'Position', [0 0 1 0.7]); % Full width, 40% of screen height

% Real Part
subplot(1,3,1)
plot(freq(1:nf), real(exx(1,1:nf)), 'ko-', 'linewidth', 2)
hold on
plot(freq(1:nf), real(exx(2,1:nf)), 'ro--', 'linewidth', 2)
plot(freq(1:nf), real(exx(3,1:nf)), 'bo:', 'linewidth', 2)
hold off
xlabel('Frequency')
ylabel('Real(Exx) (V/m)')
legend('Broadside', '45 degrees', 'Endfire')
set(gca, 'fontsize', 18)

% Imaginary Part
subplot(1,3,2)
plot(freq(1:nf), imag(exx(1,1:nf)), 'ko-', 'linewidth', 2)
hold on
plot(freq(1:nf), imag(exx(2,1:nf)), 'ro--', 'linewidth', 2)
plot(freq(1:nf), imag(exx(3,1:nf)), 'bo:', 'linewidth', 2)
hold off
xlabel('Frequency')
ylabel('Imag(Exx) (V/m)')
legend('Broadside', '45 degrees', 'Endfire')
set(gca, 'fontsize', 18)

% Magnitude
subplot(1,3,3)
loglog(freq(1:nf), abs(exx(1,1:nf)), 'ko-', 'linewidth', 2)
hold on
loglog(freq(1:nf), abs(exx(2,1:nf)), 'ro--', 'linewidth', 2)
loglog(freq(1:nf), abs(exx(3,1:nf)), 'bo:', 'linewidth', 2)
hold off
xlabel('Frequency')
ylabel('|Exx| (V/m)')
legend('Broadside', '45 degrees', 'Endfire')
set(gca, 'fontsize', 18)

% Super Title
sgtitle(['Reflection response, d = ', num2str(d,2), ' cm, \epsilon_{r;1}=', num2str(epsr1,2), ', \epsilon_{r;2}=', num2str(epsr2,2), ' r=', num2str(r,2)], 'FontSize', 20, 'FontWeight', 'bold')

% Save the figure
print('Exx', '-dpng')

% Output all field components and their absolute values and phase values to csv for all frequencies
output_data = zeros(nf, 13); % Preallocate for frequency, real, imag, abs, phase for each angle
output_data(:, 1) = freq; % Frequency

for ip = 1:3
  % Compute absolute values and phase
  abs_exx = abs(exx(ip, 1:nf));
  phase_exx = angle(exx(ip, 1:nf)) * 180 / pi; % Convert phase to degrees

  % Combine data into the matrix
  output_data(:, 2 + (ip - 1) * 4) = real(exx(ip, 1:nf)); % Real part
  output_data(:, 3 + (ip - 1) * 4) = imag(exx(ip, 1:nf)); % Imaginary part
  output_data(:, 4 + (ip - 1) * 4) = abs_exx; % Magnitude
  output_data(:, 5 + (ip - 1) * 4) = phase_exx; % Phase
end

% Write to CSV file
csv_filename = 'Evert_WB_250.csv';
header = {'Frequency', 'Real_Endfire', 'Imag_Endfire', 'Abs_Endfire', 'Phase_Endfire', ...
      'Real_45deg', 'Imag_45deg', 'Abs_45deg', 'Phase_45deg', ...
      'Real_B', 'Imag_B', 'Abs_B', 'Phase_B'};
writecell(header, csv_filename); % Write the header
writematrix(output_data, csv_filename, 'Delimiter', ',', 'WriteMode', 'append'); % Append the data

% Output all field components and their absolute values and phase values to csv for all frequencies
output_data = zeros(nf, 13); % Preallocate for frequency, real, imag, abs, phase for each angle
output_data(:, 1) = freq; % Frequency

for ip = 1:3
  % Compute absolute values and phase
  abs_exx_o = abs(exx_o(ip, 1:nf));
  phase_exx_o = angle(exx_o(ip, 1:nf)) * 180 / pi; % Convert phase to degrees

  % Combine data into the matrix
  output_data(:, 2 + (ip - 1) * 4) = real(exx_o(ip, 1:nf)); % Real part
  output_data(:, 3 + (ip - 1) * 4) = imag(exx_o(ip, 1:nf)); % Imaginary part
  output_data(:, 4 + (ip - 1) * 4) = abs_exx_o; % Magnitude
  output_data(:, 5 + (ip - 1) * 4) = phase_exx_o; % Phase
end

% Write to CSV file
csv_filename = 'Evert_WB_250_o.csv';
header = {'Frequency', 'Real_Endfire', 'Imag_Endfire', 'Abs_Endfire', 'Phase_Endfire', ...
      'Real_45deg', 'Imag_45deg', 'Abs_45deg', 'Phase_45deg', ...
      'Real_B', 'Imag_B', 'Abs_B', 'Phase_B'};
writecell(header, csv_filename); % Write the header
writematrix(output_data, csv_filename, 'Delimiter', ',', 'WriteMode', 'append'); % Append the data
