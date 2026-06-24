% TEST_MULTI_N_RESONANCE
% No-antenna test: overlay resonance lines for multiple toroidal mode
% numbers (n) and frequencies on a single Pzeta-lambda plane.
%
% Reference frequencies are drawn from a past GTC multi-n benchmark:
%   n = 3  RSAE   f = 59.46 kHz
%   n = 4  RSAE   f = 66.07 kHz
%   n = 6  TAE    f = 97.50 kHz
% (see memory: transport-study, IOP server 2026-06-07 audit)
%
% Procedure:
%   1. Run an unperturbed PLam scan via run_paper.m  -> char_freq grid.
%   2. For each (n_k, omega_k), compute l_res = (omega_k - n_k*omega_phi)
%      ./ omega_b on the Pzeta-lambda grid.
%   3. Draw integer-p resonance contours (p = -3..3) colored per channel
%      on one figure.
%
% Output: ../output/test_multi_n_resonance.fig + .png

clc; clear; clear global; close all;

global ps_option passing_option ant_option

tokamak         = 'Analytical equilibrium';
ps_option       = 'PLam';
passing_option  = 'co-passing';
ant_option      = 'off';

% Run unperturbed pipeline; this populates char_freq in the caller scope.
run run_paper.m;

% Channel definitions: (label, n_tor, f_kHz, plot color)
chs = struct( ...
    'label', {'n=3 RSAE', 'n=4 RSAE', 'n=6 TAE'}, ...
    'n_tor', {3, 4, 6}, ...
    'f_kHz', {59.46, 66.07, 97.50}, ...
    'color', {[0.85 0.10 0.10], [0.10 0.55 0.85], [0.20 0.65 0.20]} );

% Poloidal harmonic indices p to draw per channel.
p_levels = -3:3;

% Figure
fh = figure('name','multi-n resonance lines (no antenna)', ...
            'unit','normalized', ...
            'position',[0.1 0.1 0.7 0.7], ...
            'color','w');
hold on; box on; grid on;
set(gca,'fontsize',16);

legend_handles = gobjects(1, numel(chs));

for k = 1:numel(chs)
    omega_k = 2*pi*chs(k).f_kHz*1000;       % rad/s, same unit as char_freq.omega_*
    n_k     = chs(k).n_tor;
    l_res   = (omega_k - n_k*char_freq.omega_phi) ./ char_freq.omega_b;

    % Draw each integer p as a contour at the level value p (Eq. 70 Bao 24)
    for p = p_levels
        [~, hC] = contour(char_freq.Pzeta_norm, char_freq.lambda, ...
                          l_res, [p p], ...
                          'LineColor', chs(k).color, 'LineWidth', 1.4);
        if p == p_levels(1)
            legend_handles(k) = hC;
        end
    end
end

xlabel('$\bar{P}_\zeta / \psi_w$','Interpreter','latex','fontsize',20);
ylabel('$\lambda = \mu B_a / E$','Interpreter','latex','fontsize',20);
if exist('energy_in','var')
    E_label = sprintf('E = %g keV', energy_in);
else
    E_label = '';
end
title(sprintf('Multi-n resonance lines  (%s, %s, %s)', ...
              E_label, passing_option, tokamak), ...
      'Interpreter','none','fontsize',16);

legend(legend_handles, {chs.label}, 'Location', 'northeast', 'fontsize', 14);

% Save
out_dir = '../output';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
savefig(fh, fullfile(out_dir, 'test_multi_n_resonance.fig'));
exportgraphics(fh, fullfile(out_dir, 'test_multi_n_resonance.png'), ...
               'Resolution', 200);

fprintf('Saved figure to %s/test_multi_n_resonance.{fig,png}\n', out_dir);
