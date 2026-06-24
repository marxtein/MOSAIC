% TEST_MULTICHANNEL_CONSERVATION
% Multi-channel antenna E'_k conservation test.
%
% Three coherent antenna channels are turned on simultaneously (chosen to
% mirror the past GTC multi-n benchmark; see memory transport-study):
%   k = 1   n = 3, RSAE-like   f = 59.46 kHz
%   k = 2   n = 4, RSAE-like   f = 66.07 kHz
%   k = 3   n = 6, TAE-like    f = 97.50 kHz
%
% Each channel uses a single m harmonic with a Gaussian radial envelope
% built by test_synthetic_channel (this is not a real eigenmode - it is a
% clean stand-in for verifying the multi-channel infrastructure and the
% conservation diagnostic).
%
% Conserved quantity diagnosed per channel:
%   E'_k(t) = E(t) - (omega_k / n_k) * qpart * P_zeta(t)
%
% Expectation:
%   - Single channel on:        E'_k drift ~ machine precision.
%   - All 3 channels on:        E'_k drifts because no single (omega, n)
%                                rotating frame removes all wave phases.
%                                Drift amplitude indicates how strongly the
%                                particle is influenced by competing channels.
%
% Output: ../output/test_multichannel_conservation.{fig,png}.

clc; clear; clear global; close all;

global eq particle unit plasma ps_option passing_option ant_option

% Use analytical equilibrium so the test has no external data dependency.
tokamak         = 'Analytical equilibrium';
ps_option       = 'PLam';
passing_option  = 'co-passing';
ant_option      = 'on';

% Set up equilibrium + units (same prologue as run_paper.m for the
% analytical branch).
plasma.R0      = 83.5;     % cm
plasma.b0      = 20125;    % gauss
particle.charge = 1.0;
particle.mass   = 1.0;
eq.qcoef       = [0.82; 1.1; 1.0];
eq.psiw        = 0.0375;
run construct_analytic_eq;
run physics_unit;

% Build 3 synthetic channels.
psi_grid  = linspace(0.0, eq.psiw, 200).';
psi_peak  = 0.55*eq.psiw;
psi_width = 0.10*eq.psiw;
amp_phi   = 2e-5;          % calibrated: perturbation visible above RK2 integration noise; orbits remain bounded over mstep=4000.

ch1 = test_synthetic_channel(3,  59.46, 5,  psi_grid, psi_peak, psi_width, amp_phi);
ch2 = test_synthetic_channel(4,  66.07, 7,  psi_grid, psi_peak, psi_width, amp_phi);
ch3 = test_synthetic_channel(6,  97.50, 11, psi_grid, psi_peak, psi_width, amp_phi);

% Wrap each into ant struct so ant_merge_channels can stack them.
ant_wrap = @(c) struct('channels', c);
ant_multi = ant_merge_channels(ant_wrap(ch1), ant_wrap(ch2), ant_wrap(ch3));

fprintf('Built multi-channel antenna: %d channels.\n', numel(ant_multi.channels));
for k = 1:numel(ant_multi.channels)
    ck = ant_multi.channels(k);
    fprintf('  ch%d: n=%d m=%d omega=%.4e (=%.2f kHz) num_modes=%d\n', ...
        k, ck.n, ck.m_modes(1), ck.omega, ...
        ck.omega/(2*pi*unit.gtc_utime)/1000, ck.num_modes);
end

% Initial conditions: a handful of representative orbits.
np_test = 6;
qpart = particle.charge; apart = particle.mass;

% Place 6 particles at different (psi, pitch) to span topology types.
psi0    = eq.psiw * [0.35 0.45 0.55 0.60 0.65 0.70];
theta0  = [0       pi/4    pi/2    3*pi/4  pi       5*pi/4];
zeta0   = zeros(1, np_test);
E0_keV  = [30      30      40      30      45       40];
lam0    = [0.4     0.5     0.6     0.7     0.5      0.3];

zpart = zeros(14, np_test);
zpart(1, :) = psi0;
zpart(2, :) = theta0;
zpart(3, :) = zeta0;
% Build u_para from E and lambda at theta0 (use b at that location).
for m = 1:np_test
    bm = qdspline.spline2d(0, psi0(m), theta0(m), eq.lsp, eq.lst, ...
                           eq.dpsi, eq.dtheta, eq.bsp);
    E_norm = E0_keV(m)*1000/unit.energy_norm;
    mu     = lam0(m)*E_norm/bm;
    zpart(6, m) = sqrt(mu);
    upara  = sqrt(2*(E_norm - mu*bm)/apart);
    % rho_|| stored in slot 4: u_para = rho_||*b*cmratio  =>  rho_|| = upara/(b*cmratio)
    cmratio = qpart/apart;
    zpart(4, m) = upara/(bm*cmratio);
end

mstep   = 4000;
tstep   = 10;     % omega_cp^-1
eq_curr = 1;
eq_gradB = 1;
amp_mod  = 1.0;   % multi-channel amplitudes already baked in via amp_phi

fprintf('Pushing %d orbits, mstep=%d ...\n', np_test, mstep);
tic;
B = push_orbit_multichannel(zpart, mstep, tstep, qpart, apart, ...
                            eq_curr, eq_gradB, ant_multi, amp_mod);
fprintf('Push took %.2f s.\n', toc);

% Diagnostics: relative drift of E'_k per channel per particle.
nc = numel(ant_multi.channels);
t_unit_s = unit.gtc_utime;           % omega_cp^-1 in seconds
t_ms = B.t * t_unit_s * 1e3;          % time in ms

% Energy in keV
E_keV     = B.E      * unit.energy_norm / 1000;
Eprime_keV = B.Eprime * unit.energy_norm / 1000;

fh = figure('name','multi-channel E_k conservation', ...
            'unit','normalized','position',[0.05 0.05 0.85 0.85], ...
            'color','w');

% Top: E(t) and P_zeta(t) for all particles (sanity)
subplot(3, 1, 1); hold on; box on; grid on;
for m = 1:np_test
    valid = ~isnan(B.E(:, m));
    plot(t_ms(valid), E_keV(valid, m), 'LineWidth', 1.2);
end
ylabel('$E$ [keV]','Interpreter','latex','fontsize',16);
title('Multi-channel antenna test: E(t), P_\zeta(t), and E''_k(t) drift', ...
      'Interpreter','tex','fontsize',16);
set(gca,'fontsize',14);

subplot(3, 1, 2); hold on; box on; grid on;
for m = 1:np_test
    valid = ~isnan(B.E(:, m));
    plot(t_ms(valid), B.Pzeta(valid, m)/eq.psiw, 'LineWidth', 1.2);
end
ylabel('$P_\zeta / \psi_w$','Interpreter','latex','fontsize',16);
set(gca,'fontsize',14);

% Bottom: per-channel E'_k drift (relative to its initial value), all
% particles, separated by channel via line color, by particle via line
% style group.
subplot(3, 1, 3); hold on; box on; grid on;
colors = {[0.85 0.10 0.10], [0.10 0.55 0.85], [0.20 0.65 0.20]};
legend_entries = {};
for k = 1:nc
    Eprime0 = Eprime_keV(1, :, k);   % 1 x np
    for m = 1:np_test
        valid = ~isnan(B.E(:, m));
        drift = Eprime_keV(valid, m, k) - Eprime0(m);
        plot(t_ms(valid), drift, 'Color', colors{k}, 'LineWidth', 1.0);
    end
    legend_entries{end+1} = sprintf('ch%d  n=%d  f=%.2f kHz', ...
        k, ant_multi.channels(k).n, ...
        ant_multi.channels(k).omega/(2*pi*unit.gtc_utime)/1000);
end
xlabel('$t$ [ms]','Interpreter','latex','fontsize',16);
ylabel('$E''_k(t) - E''_k(0)$ [keV]','Interpreter','latex','fontsize',16);
set(gca,'fontsize',14);

% one legend handle per channel (first particle of that channel)
hL = gobjects(1, nc);
for k = 1:nc
    hL(k) = plot(nan, nan, 'Color', colors{k}, 'LineWidth', 2.0);
end
legend(hL, legend_entries, 'Location', 'best', 'fontsize', 12);

% Summary numbers
fprintf('\nE''_k drift summary (keV):\n');
for k = 1:nc
    drift_max = 0;
    for m = 1:np_test
        valid = ~isnan(B.E(:, m));
        d = max(abs(Eprime_keV(valid, m, k) - Eprime_keV(1, m, k)));
        drift_max = max(drift_max, d);
    end
    fprintf('  ch%d (n=%d, f=%.2f kHz):  max |E''_k - E''_k(0)| = %.3e keV\n', ...
        k, ant_multi.channels(k).n, ...
        ant_multi.channels(k).omega/(2*pi*unit.gtc_utime)/1000, drift_max);
end

out_dir = '../output';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
savefig(fh, fullfile(out_dir, 'test_multichannel_conservation.fig'));
exportgraphics(fh, fullfile(out_dir, 'test_multichannel_conservation.png'), ...
               'Resolution', 200);
fprintf('Saved figure to %s/test_multichannel_conservation.{fig,png}\n', out_dir);
