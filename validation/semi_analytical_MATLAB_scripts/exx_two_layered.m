% routine to compute the x-component of the electric field for an
% x-directed electric dipole source, receiver is assumed to be on the
% ground and the source very close to but just above the ground

% if receiver is not on but above the ground, we need to include the direct
% contributions, but I assume we will not use that here
clear variables
clc
% height of source above surface
z=0.025;
% number of offsets
nr=76;
% radial offset along surface
rr=linspace(0.1, 1, nr);
% relative espilon and sigma values of first layer and half space
epsr1=4;
sigma1=1e-4;
epsr2=9;
sigma2=1e-3;
% thickness of layer 1
d = 1;
% free space velocity and magnetic permeability
c0=299792458;
mu=4e-7*pi;
% frequency value
freq=1e8;
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
exxi=zeros(nphi,nr);
exxr=zeros(nphi,nr);
% Preallocate results
output = zeros(length(rr), 1 + 4 * length(phi)); % Columns: r, Real, Imag, Abs, Phase for each angle
% set tolerance for integrations
atol=1e-16;
rtol=1e-8;
% loop over offsets
for ir=1:nr
    r=rr(ir);
% bounds for integrals
    a=imag(sqrt(eta0*smu));
    b=imag(sqrt(eta1*smu));
    mxbnd=sqrt((40/abs(z))^2+abs(eta1*smu));
    vc=1; % J0-integral
    vc1=integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,0,a,'reltol',rtol);
    vc1=vc1 + integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,a,b,'reltol',rtol);
    vc1=vc1 + integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,b,mxbnd,'reltol',rtol);
    vcr1=integral(@(kappa)hfhorxr(kappa,eta0,eta1,eta2,smu,r,z,d,vc)...
                ,0,a,'reltol',rtol);
    vcr1=vcr1 + integral(@(kappa)hfhorxr(kappa,eta0,eta1,eta2,smu,r,z,d,vc)...
                ,a,mxbnd,'reltol',rtol);
    vc=2; %J2-integral
    vc2=integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,0,a,'reltol',rtol);
    vc2=vc2 + integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,a,b,'reltol',rtol);
    vc2=vc2 + integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,b,mxbnd,'reltol',rtol);
    vcr2=integral(@(kappa)hfhorxr(kappa,eta0,eta1,eta2,smu,r,z,d,vc)...
                ,0,a,'reltol',rtol);
    vcr2=vcr2 + integral(@(kappa)hfhorxr(kappa,eta0,eta1,eta2,smu,r,z,d,vc)...
                ,a,mxbnd,'reltol',rtol);
% compute electric field from the sub-results for the three angles
    exxi(:,ir)=-(vc1-cos(2*phi).*vc2)/(4*pi);
    exxr(:,ir)=(vcr1-cos(2*phi).*vcr2)/(8*pi);
end
exx=exxi+exxr;
% compute maximum amplitude
mx1=max(max(abs(exx)));

% Saving to CSV using output
csv_filename = sprintf('Exx_4_9_%.0fMHz_NR_s.csv', (freq/1e6));
header = {'Offset', ...
    'Endfire_Amplitude','Endfire_Phase','Endfire_Real','Endfire_Imag', ...
    'Broadside_Amplitude','Broadside_Phase','Broadside_Real','Broadside_Imag', ...
    'Oblique_Amplitude','Oblique_Phase','Oblique_Real','Oblique_Imag'};

output(1:nr, 1)  = rr;                % Offset
% Endfire (phi=0, index 1)
output(1:nr, 2)  = abs(exx(1,1:nr));      % Endfire_Amplitude
output(1:nr, 3)  = angle(exx(1,1:nr)); % Endfire_Phase (radians)
output(1:nr, 4)  = real(exx(1,1:nr));     % Endfire_Real
output(1:nr, 5)  = imag(exx(1,1:nr));     % Endfire_Imag
% Broadside (phi=90, index 3)
output(1:nr, 6)  = abs(exx(3,1:nr));      % Broadside_Amplitude
output(1:nr, 7)  = angle(exx(3,1:nr)); % Broadside_Phase (radians)
output(1:nr, 8)  = real(exx(3,1:nr));     % Broadside_Real
output(1:nr, 9)  = imag(exx(3,1:nr));     % Broadside_Imag
% Oblique (phi=45, index 2)
output(1:nr,10)  = abs(exx(2,1:nr));      % Oblique_Amplitude
output(1:nr,11)  = angle(exx(2,1:nr)); % Oblique_Phase (radians)
output(1:nr,12)  = real(exx(2,1:nr));     % Oblique_Real
output(1:nr,13)  = imag(exx(2,1:nr));     % Oblique_Imag

writecell(header, csv_filename); % Write the header
writematrix(output, csv_filename, 'Delimiter', ',', 'WriteMode', 'append'); % Append the data

figure(1)
semilogy(rr,abs(exx(1,1:nr)),'r',rr,abs(exx(2,1:nr)),'g',rr,abs(exx(3,1:nr)),'b','linewidth',2)
axis([0 5 mx1/1e5 1e3])
xlabel('offset (m)')
ylabel('electric field amplitude (V/m)')
legend('endfire','45 degrees','broadside')
title(['\epsilon_{r;1,2}=(',num2str(epsr1,2),',',num2str(epsr2,2),'), \sigma_{1,2}=(',num2str(sigma1,4),',',num2str(sigma2,3),'), d=2 m'])
set(gca,'fontsize',18)

figure(2)
plot(rr,angle(exx(1,1:nr)),'r',rr,angle(exx(2,1:nr)),'g',rr,angle(exx(3,1:nr)),'b','linewidth',2)
axis([0 5 -pi pi])
xlabel('offset (m)')
ylabel('electric field amplitude (V/m)')
legend('endfire','45 degrees','broadside')
title(['\epsilon_{r;1,2}=(',num2str(epsr1,2),',',num2str(epsr2,2),'), \sigma_{1,2}=(',num2str(sigma1,4),',',num2str(sigma2,3),'), d=2 m'])
set(gca,'fontsize',18)
