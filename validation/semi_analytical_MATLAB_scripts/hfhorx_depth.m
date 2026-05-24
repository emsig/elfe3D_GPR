function exyz = hfhorx_depth(kappa, eta0, eta1, smu, r, a, z1, vc)
% HFHORX_DEPTH  Hankel integrand for HED field at depth z1 in a half-space.
%
%   Computes the Ex kernel for an x-directed HED at height a in air,
%   evaluated at depth z1 > 0 inside layer 1 (half-space model).
%
% Usage:
%   exyz = hfhorx_depth(kappa, eta0, eta1, smu, r, a, z1, vc)
%
%   kappa  : horizontal wavenumber [1/m]
%   eta0   : admittivity of air    = 2i*pi*freq/(c0^2*mu)
%   eta1   : admittivity of layer1 = eta0*epsr1 + sigma1
%   smu    : j*omega*mu
%   r      : horizontal source-receiver distance [m]
%   a      : source height above surface [m]
%   z1     : evaluation depth below surface [m]  (positive)
%   vc     : 1 = J0 integral,  2 = J2 integral
%
% Derivation:
%   The tangential E-field is continuous at the air/layer-1 interface.
%   At z = 0- (just below interface), E_tan equals E_tan at z = 0+
%   (surface field from hfhorx). Below the interface, the field propagates
%   downward as exp(-grg1 * z1) — the only wave in an infinite half-space
%   is the transmitted downward wave, so no upgoing component exists.
%
%   Therefore:
%       kernel_depth = hfhorx(kappa, ..., a, vc) * exp(-grg1 * z1)
%
%   This is exact for a half-space (single interface, no reflections below).
%
% Scattered field from buried dipole to surface (reciprocity):
%   For receivers AT the surface, set a = 0 (surface source in reciprocal
%   problem) and z1 = sphere_depth. By source-receiver reciprocity:
%       G_xx(r_rx, z=0; r_sphere, z=-z1)  =  G_xx(r_sphere, z=-z1; r_rx, z=0)
%   The right-hand side is hfhorx_depth with a=0, z1=sphere_depth.

grg1  = sqrt(smu * eta1 + kappa.^2);
exyz  = hfhorx(kappa, eta0, eta1, smu, r, a, vc) .* exp(-grg1 .* z1);
