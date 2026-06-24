function ch = test_synthetic_channel(n_tor, f_kHz, m_tor, ...
                                     psi_grid, psi_peak, psi_width, amplitude)
% TEST_SYNTHETIC_CHANNEL  Build a single-m analytical antenna channel for
% testing. The radial profile is a real-valued Gaussian centered at
% psi_peak with width psi_width; only the imaginary part is zero. This is
% a stripped-down stand-in for an eigenmode, fine for verifying multi-
% channel infrastructure and conservation behavior, not for physics
% validation.
%
% Inputs
%   n_tor      : toroidal mode number (scalar).
%   f_kHz      : wave frequency in kHz.
%   m_tor      : poloidal mode number (scalar; only one m harmonic).
%   psi_grid   : radial grid (column vector).
%   psi_peak   : envelope center.
%   psi_width  : envelope sigma.
%   amplitude  : peak amplitude for delta-phi profile (use ~1e-4 to start).
%
% Output ch: channel struct compatible with eval_perturb_field.

global plasma unit %#ok<GVMIS>
if isempty(unit) || ~isfield(unit, 'gtc_utime')
    error('test_synthetic_channel:no_unit', ...
        'physics_unit must run first so unit.gtc_utime is available.');
end

psi_grid = psi_grid(:);
N        = numel(psi_grid);

env = amplitude .* exp(-0.5*((psi_grid - psi_peak)./psi_width).^2);
denv_dpsi = -((psi_grid - psi_peak)./(psi_width^2)) .* env;

% Single m harmonic: real profile, zero imaginary
phi_real_cell    = { env.' };
phi_imag_cell    = { zeros(1, N) };
apara_real_cell  = { 0.1*env.' };     % small parallel A_|| for completeness
apara_imag_cell  = { zeros(1, N) };
dptdp_real_cell  = { denv_dpsi.' };
dptdp_imag_cell  = { zeros(1, N) };
dapdp_real_cell  = { 0.1*denv_dpsi.' };
dapdp_imag_cell  = { zeros(1, N) };

omega_mosaic = 2*pi*(f_kHz*1000)*unit.gtc_utime;

ch = ant_make_channel(omega_mosaic, n_tor, m_tor, psi_grid, ...
    phi_real_cell, phi_imag_cell, apara_real_cell, apara_imag_cell, ...
    dptdp_real_cell, dptdp_imag_cell, dapdp_real_cell, dapdp_imag_cell);
end
