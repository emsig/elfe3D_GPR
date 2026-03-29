function exyz=hfhorxr(kappa,eta0,eta1,eta2,smu,r,a,d,vc)
% vertical wave numbers in air and the two subsurface layers
grg0=sqrt(smu*eta0+kappa.^2);
grg1=sqrt(smu*eta1+kappa.^2);
grg2=sqrt(smu*eta2+kappa.^2);
% compute local reflection coefficients of the TE and TM modes, r0 is
% surface reflection and r1 is the subsurface boundary reflection
r0te=(grg0-grg1)./(grg0+grg1);
r1te=(grg1-grg2)./(grg1+grg2);
r0tm=(eta1*grg0-eta0*grg1)./(eta1*grg0+eta0*grg1);
r1tm=(eta2*grg1-eta1*grg2)./(eta2*grg1+eta1*grg2);
% compute subsurface reflection response including mutliples against ground
% surface
glre=(1-r0te.^2).*r1te.*exp(-2*grg1*d)./(1+r0te.*r1te.*exp(-2*grg1*d));
glrm=(1-r0tm.^2).*r1tm.*exp(-2*grg1*d)./(1+r0tm.*r1tm.*exp(-2*grg1*d));
% Bessel function of order 0 and 2
J0kc=besselj(0,kappa*r).*kappa;
J2kc=besselj(2,kappa*r).*kappa;
% electric field part depending on the contribution connected to the J0 or
% the J2 integrals.
if vc == 2  % J2-integral
        exyz=J2kc.*exp(-grg0*a).*(grg0.*glrm/eta0+smu*glre./grg0);
else % J0-integral
        exyz=J0kc.*exp(-grg0*a).*(grg0.*glrm/eta0-smu*glre./grg0);
end
% done
% 
% for a more general model with N+1 layers we use this recursive scheme
% in that case we need eta0 and eta(1:nl+1) and we replace lines 3 to 15
% above with the lines below
% 
% we start on the bottom boundary where the total reflection response is
% equal to the local reflection response
% grg1=sqrt(smu*eta(nl)+kappa.^2);
% grg2=sqrt(smu*eta(nl+1)+kappa.^2);
% grtm=(eta(nl+1)*grg1-eta(nl)*grg2)./(eta(nl+1)*grg1+eta(nl)*grg2);
% grte=(grg1-grg2)./(grg1+grg2);
% % loop over all subsurface boundaries
% for il=nl-1:-1:1
%     grg2=grg1;
%     grg1=sqrt(eta(il)*smu+kappa.^2);
%     rtm=(eta(il+1)*grg1-eta(il)*grg2)./(eta(il+1)*grg1+eta(il)*grg2);
%     rte=(grg1-grg2)./(grg1+grg2);
%     grtm=(rtm+grtm.*exp(-2*grg2*d(il+1)))./...
%       (1+rtm.*grtm.*exp(-2*grg2*d(il+1)));
%     grte=(rte+grte.*exp(-2*grg2*d(il+1)))./...
%       (1+rte.*grte.*exp(-2*grg2*d(il+1)));
% end
% % now we are just above the top subsurface boundary and move to surface
% grg0=sqrt(eta0*smu+kappa.^2);
% r0te=(grg0-grg1)./(grg0+grg1);
% r0tm=(eta(1)*grg0-eta0*grg1)./(eta(1)*grg0+eta0*grg1);
% glre=(1-r0te.^2).*r1te.*exp(-2*grg1*d(1))./(1+r0te.*rte.*exp(-2*grg1*d(1)));
% glrm=(1-r0tm.^2).*r1tm.*exp(-2*grg1*d(1))./(1+r0tm.*rtm.*exp(-2*grg1*d(1)));
