function exyz=hfhorx(kappa,eta0,eta1,smu,r,a,vc)
grg0=sqrt(smu*eta0+kappa.^2);
grg1=sqrt(smu*eta1+kappa.^2);
J0kc=besselj(0,kappa*r).*kappa;
J2kc=besselj(2,kappa*r).*kappa.^3;
if vc == 2  % J2-integral
        exyz=J2kc.*exp(-grg0*a)./(eta1.*grg0+eta0.*grg1);
else % J0-integral
        exyz=J0kc.*exp(-grg0*a).*(2*smu./(grg0+grg1) + ...
            kappa.^2./(eta1.*grg0+eta0.*grg1));
end
