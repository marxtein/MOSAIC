function ch = ant_make_channel(omega, n_tor, m_modes, psi_grid, ...
    phi_real, phi_imag, apara_real, apara_imag, ...
    dptdp_real, dptdp_imag, dapdp_real, dapdp_imag)
% ANT_MAKE_CHANNEL  Pack one antenna channel for multi-frequency runs.
%
% A channel = a single (omega, n_tor) coherent perturbation with a list
% of poloidal harmonics m_modes and their radial profiles, plus
% radial derivatives. Field structure mirrors the legacy ant.* fields so
% existing inline loops can be replaced by eval_perturb_field with no
% numerical change at N_channels = 1.

ch.omega       = omega;
ch.n           = n_tor;
ch.m_modes     = m_modes(:).';
ch.num_modes   = numel(ch.m_modes);
ch.psi         = psi_grid(:);

ch.phi_real    = phi_real;
ch.phi_imag    = phi_imag;
ch.apara_real  = apara_real;
ch.apara_imag  = apara_imag;
ch.dptdp_real  = dptdp_real;
ch.dptdp_imag  = dptdp_imag;
ch.dapdp_real  = dapdp_real;
ch.dapdp_imag  = dapdp_imag;
end
