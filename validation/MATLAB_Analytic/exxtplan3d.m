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
nr=100;
% radial offset along surface
rr=0.05*(1:nr);
% relative espilon and sigma values of first layer and half space
epsr1=4;
sigma1=1e-4;
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
% modify eta0 for half space test
eta0=eta1;
% preset arrays to zero that occur in loops
% electric fields on surface and reflection response
exx=zeros(nphi,nr);
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
    % vc1=vc1 + integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
    %             ,a,b,'reltol',rtol);
    vc1=vc1 + integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,a,mxbnd,'reltol',rtol);
    vc=2; %J2-integral
    vc2=integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,0,a,'reltol',rtol);
    % vc2=vc2 + integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
    %             ,a,b,'reltol',rtol);
    vc2=vc2 + integral(@(kappa)hfhorx(kappa,eta0,eta1,smu,r,z,vc)...
                ,a,mxbnd,'reltol',rtol);
% compute electric field from the sub-results for the three angles
exx(:,ir)=-(vc1-cos(2*phi).*vc2)/(4*pi);
end
% compute maximum amplitude
mx1=max(max(abs(exx)));

% whole space test
eta0=eta1;
grg0=sqrt(smu.*eta0);
% endfire
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
figure(1)
semilogy(rr,abs(exx(1,1:nr)),'r',rr,abs(exx(2,1:nr)),'g',rr,abs(exx(3,1:nr)),'b',...
    rr,abs(exxe),'k--',rr,abs(exxd),'r--',rr,abs(exxb),'y--','linewidth',2)
axis([0 5 mx1/1e5 mx1])
legend('endfire-numerical','45 degrees-numerical','broadside-numerical','endfire-exact','45 degrees-exact','broadside-exact')
xlabel('offset (m)')
ylabel('electric field amplitude (V/m)')
title(['Reflection response, \epsilon_{r;1}=',num2str(epsr1,2)])
set(gca,'fontsize',18)
nrmsee=abs(exx(1,:)-exxe)./abs(exxe);
nrmsed=abs(exx(2,:)-exxd)./abs(exxd);
nrmseb=abs(exx(3,:)-exxb)./abs(exxb);creat
figure(2)
semilogy(rr,nrmsee,'r',rr,nrmsed,'g',rr,nrmseb,'b','linewidth',2)
axis([0 5 1e-10 1e-2])
legend('endfire','45 degrees','broadside')
xlabel('offset (m)')
ylabel('normalized error in amplitude')
title('amplitude of normalized difference in electric field')
set(gca,'fontsize',18)
figure(3)
plot(rr,angle(exx(1,:))-angle(exxe),'r',rr,angle(exx(2,:))-angle(exxd),'g',rr,angle(exx(3,:))-angle(exxb),'b','linewidth',2)
legend('endfire','45 degrees','broadside')
xlabel('offset (m)')
ylabel('normalized error in amplitude')
title('difference in electric field phases')
set(gca,'fontsize',18)
