function exyz = hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, a, d, z1, vc)
% HFHORXR_DEPTH  Hankel integrand for HED field at depth z1 in layer 1
%                of a two-layer earth model.
%
% Usage:
%   exyz = hfhorxr_depth(kappa, eta0, eta1, eta2, smu, r, a, d, z1, vc)
%
%   kappa  : horizontal wavenumber [1/m]
%   eta0   : admittivity of air
%   eta1   : admittivity of layer 1
%   eta2   : admittivity of layer 2
%   smu    : j*omega*mu
%   r      : horizontal distance [m]
%   a      : source height above surface [m]
%   d      : thickness of layer 1 [m]
%   z1     : evaluation depth in layer 1 [m],  0 < z1 < d
%   vc     : 1 = J0 integral,  2 = J2 integral
%
% Approximation:
%   Uses the same exp(-grg1*z1) extension as the half-space case:
%       kernel = hfhorxr(kappa, ..., d, vc) * exp(-grg1 * z1)
%
%   This is the additional REFLECTION contribution from layer 2.
%   The total depth field is: hfhorx_depth + hfhorxr_depth.
%
% %% UNCERTAIN — this approximation treats the layer-2 reflection as
%   propagating to the surface and then back down to depth z1 via
%   exp(-grg1*z1). In reality, the downgoing reflection field at depth z1
%   is modified by multiple reflections between the surface and layer2 that
%   accumulate extra exp(-grg1*z1) phase shifts. The error grows as z1/d
%   increases. The approximation is most reliable when z1 << d.
%   For the nominal parameters (z1=0.7 m, d=1.0 m, z1/d=0.7), treat
%   results with caution; validate against FEM output when available.
%
% The derivation would be exact if layer 2 were absent (hfhorxr → 0),
% reducing to hfhorx_depth. The layer-2 correction adds the reflected
% wave contribution which modifies the standing-wave pattern in layer 1.

grg1  = sqrt(smu * eta1 + kappa.^2);
exyz  = hfhorxr(kappa, eta0, eta1, eta2, smu, r, a, d, vc) .* exp(-grg1 .* z1);
