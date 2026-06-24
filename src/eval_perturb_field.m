function F = eval_perturb_field(channels, pdum, tdum, zdum, t_now, amp_mod)
% EVAL_PERTURB_FIELD  Multi-channel antenna perturbation evaluator.
%
% Inputs
%   channels : struct array, channels(k) carries (omega, n, m_modes,
%              num_modes, psi, phi_{real,imag}, apara_{real,imag},
%              dptdp_{real,imag}, dapdp_{real,imag}).
%   pdum     : psi (scalar).
%   tdum     : poloidal angle theta (scalar).
%   zdum     : toroidal angle zeta (scalar).
%   t_now    : absolute time entering the wave phase (= tstep*(istep+0.5*irk)
%              in the legacy poincare_perturb convention).
%   amp_mod  : scalar global amplitude modulator (matches legacy amp_mod).
%
% Output struct F (sum over channels):
%   F.dptdp, F.dptdt, F.dptdz       (delta-phi and its derivatives)
%   F.apara, F.dapdp, F.dapdt, F.dapdz, F.paparapt
%
% Plus F.per_channel(k).{dptdp,dptdt,dptdz,apara,dapdp,dapdt,dapdz,paparapt}
% which carries the same eight quantities split by channel - this is what
% downstream diagnostics use to compute per-channel conserved quantities
% like E'_k = E - (omega_k/n_k) * P_zeta.
%
% At numel(channels) == 1 this function is bit-identical to the legacy
% inline loops in poincare_perturb.m / cal_char_freq_perturb.m provided
% the caller passes t_now = tstep*(istep+0.5*irk).

nc = numel(channels);
zero_pc = struct('dptdp', 0, 'dptdt', 0, 'dptdz', 0, ...
                 'apara', 0, 'dapdp', 0, 'dapdt', 0, 'dapdz', 0, ...
                 'paparapt', 0);
F.per_channel = repmat(zero_pc, 1, nc);

F.dptdp    = 0;  F.dptdt = 0;  F.dptdz = 0;
F.apara    = 0;  F.dapdp = 0;  F.dapdt = 0;  F.dapdz = 0;
F.paparapt = 0;

for k = 1:nc
    ch = channels(k);

    % radial index: legacy nearest-neighbor lookup on ch.psi
    if pdum > ch.psi(end)
        ind = numel(ch.psi);
    elseif pdum < ch.psi(1)
        ind = 1;
    else
        [~, ind] = min(abs(pdum - ch.psi));
    end

    dptdp_k = 0; dptdt_k = 0; dptdz_k = 0;
    apara_k = 0; dapdp_k = 0; dapdt_k = 0; dapdz_k = 0; paparapt_k = 0;

    omega_k = ch.omega;
    n_k     = ch.n;

    for i = 1:ch.num_modes
        m_i   = ch.m_modes(i);
        phase = -m_i*tdum + n_k*zdum - omega_k*t_now;
        cph   = cos(phase);
        sph   = sin(phase);

        phi_r    = amp_mod * ch.phi_real{i}(ind);
        phi_i    = amp_mod * ch.phi_imag{i}(ind);
        dphi_r   = amp_mod * ch.dptdp_real{i}(ind);
        dphi_i   = amp_mod * ch.dptdp_imag{i}(ind);
        apara_r  = amp_mod * ch.apara_real{i}(ind);
        apara_i  = amp_mod * ch.apara_imag{i}(ind);
        dapara_r = amp_mod * ch.dapdp_real{i}(ind);
        dapara_i = amp_mod * ch.dapdp_imag{i}(ind);

        dptdp_k = dptdp_k + dphi_r*cph - dphi_i*sph;
        dptdt_k = dptdt_k + m_i*phi_r*sph + m_i*phi_i*cph;
        dptdz_k = dptdz_k - n_k*phi_r*sph - n_k*phi_i*cph;

        apara_k = apara_k + apara_r*cph - apara_i*sph;
        dapdp_k = dapdp_k + dapara_r*cph - dapara_i*sph;
        dapdt_k = dapdt_k + m_i*apara_r*sph + m_i*apara_i*cph;
        dapdz_k = dapdz_k - n_k*apara_r*sph - n_k*apara_i*cph;
        paparapt_k = paparapt_k + omega_k*apara_r*sph + omega_k*apara_i*cph;
    end

    F.per_channel(k).dptdp    = dptdp_k;
    F.per_channel(k).dptdt    = dptdt_k;
    F.per_channel(k).dptdz    = dptdz_k;
    F.per_channel(k).apara    = apara_k;
    F.per_channel(k).dapdp    = dapdp_k;
    F.per_channel(k).dapdt    = dapdt_k;
    F.per_channel(k).dapdz    = dapdz_k;
    F.per_channel(k).paparapt = paparapt_k;

    F.dptdp    = F.dptdp    + dptdp_k;
    F.dptdt    = F.dptdt    + dptdt_k;
    F.dptdz    = F.dptdz    + dptdz_k;
    F.apara    = F.apara    + apara_k;
    F.dapdp    = F.dapdp    + dapdp_k;
    F.dapdt    = F.dapdt    + dapdt_k;
    F.dapdz    = F.dapdz    + dapdz_k;
    F.paparapt = F.paparapt + paparapt_k;
end
end
