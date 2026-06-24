function ant = ant_attach_channel(ant, n_tor)
% ANT_ATTACH_CHANNEL  Wrap the legacy ant.* fields into ant.channels(1).
%
% Call this after ant_interface(). Optionally pass n_tor (toroidal mode
% number) if not already set on the loaded mat file.
%
% After this call, ant.channels(1) carries all data the multi-channel
% evaluator eval_perturb_field needs. Legacy ant.* top-level fields are
% kept intact, so existing callers (poincare_perturb, cal_char_freq_perturb)
% continue to work without modification.

if nargin >= 2 && ~isempty(n_tor)
    ant.n = n_tor;
end
if ~isfield(ant, 'n') || isempty(ant.n)
    error('ant_attach_channel:missing_n', ...
        'ant.n (toroidal mode number) must be set before attaching channel.');
end
if ~isfield(ant, 'omega') || isempty(ant.omega)
    error('ant_attach_channel:missing_omega', ...
        'ant.omega must be set before attaching channel.');
end

ant.channels(1) = ant_make_channel( ...
    ant.omega, ant.n, ant.m_modes, ant.psi, ...
    ant.phi_real,   ant.phi_imag,   ant.apara_real, ant.apara_imag, ...
    ant.dptdp_real, ant.dptdp_imag, ant.dapdp_real, ant.dapdp_imag);
end
