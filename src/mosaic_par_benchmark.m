% This file is part of MOSAIC version 1.0
% MOSAIC is released under the GNU General Public License v3.0 (GPLv3):
%
% Copyright (c) 2026 Jian Bao (jbao@iphy.ac.cn)
% Institute of Physics, Chinese Academy of Sciences
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================== Parameter Setup ==================
del_t = tstep;
del_t_gtc = 0.5767/2;
xylabelsize = 18;

% Two particles to compare
m_list = [22, 19];   % MOSAIC particle index
col_offset = [0, 2]; % subplot column offset

% Figure setup
figure('position',[200, 100, 1800, 800],...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','normal',...
    'DefaultAxesLineWidth',1,...
    'DefaultAxesTickLength',[0.01,0.01]);
labels = {'(a1)','(a2)','(a3)','(a4)', ...
          '(b1)','(b2)','(b3)','(b4)'};
type={'passing','trapped'};
%% ================== Loop Over Particles ==================
for kk = 1:2
    if kk==1
    path = '..\input\gtc_benchmark\m=22\';
    filename = 'poincare.txt';
    run read_para_w.m
    P = read_poincare_particles(path, filename);
    del_t_gtc = 0.2883;
    else
    path = '..\input\gtc_benchmark\m=19\';
    filename = 'poincare.txt';
    run read_para_w.m
    P = read_poincare_particles(path, filename);
    del_t_gtc = 0.5767;
    end

    m = m_list(kk); % current MOSAIC particle
    offset =col_offset(kk)*4/2;

    %% ----- Find corresponding particle in P (GTC) -----
    target = max(max(A.R(:,m)));
    max_R_per_element = arrayfun(@(x) max(x.R(:)), P);
    differences = abs(max_R_per_element - target);
    [~, i] = min(differences); % matched particle index in P

    %% ----- Extract A (MOSAIC) data -----
    A_R = squeeze(A.R(:, m) * spdata.rmaj/100);
    A_Z = squeeze(A.Z(:, m) * spdata.rmaj/100);
    t_A = (1:length(A_R)) * del_t;

    % Extract toroidal angle (adjust index if needed)
    zeta_col_index = 3;
    A_zeta = squeeze(A.part_pos(:, zeta_col_index, m));

    %% ================== 1. R-Z Orbit ==================
    subplot(2,4,1+offset)
    daspect([1 1 1]);
    %text(0.02,0.05,labels{1+offset}, ...
    %    'Units','normalized','FontSize',17);hold on
    % Plot equilibrium boundary (no legend)
    [~,ind] = min(abs(A.init.psi(m) - eq.psi));
    plot(squeeze(eq.xsp(1,ind,:)*spdata.rmaj/100),...
        squeeze(eq.zsp(1,ind,:)*spdata.rmaj/100),...
        'k--','linewidth',1.5,'HandleVisibility','off'); hold on;

    plot(squeeze(eq.xsp(1,end,:)*spdata.rmaj/100),...
        squeeze(eq.zsp(1,end,:)*spdata.rmaj/100),...
        'k-','linewidth',2,'HandleVisibility','off'); hold on;

    plot(eq.xsp(1,1,1)*spdata.rmaj/100,0,'k.','markersize',15,'HandleVisibility','off');hold on;
    % Plot MOSAIC vs GTC
    plot(A_R, A_Z, 'r-', 'LineWidth', 2, ...
        'DisplayName', ['MOSAIC ', type{kk}]); hold on;

    h=plot(P(i).R * spdata.rmaj/100, ...
         P(i).Z * spdata.rmaj/100, ...
         'b--', 'LineWidth', 2, ...
         'DisplayName', ['GTC ', type{kk}]);
    xlabel('$X$(m)','FontSize',xylabelsize,'Interpreter','latex');
    ylabel('$Z$(m)','FontSize',xylabelsize,'Interpreter','latex');
    axis equal; grid on;
    xlim([0.495,2.522])
    legend('Location','northwest')

    %% ================== 2. R vs Time ==================
    subplot(2,4,2+offset)
    %text(0.02,0.05,labels{2+offset}, ...
    %    'Units','normalized','FontSize',17);hold on

    plot(t_A, A_R, 'r-', 'LineWidth', 2, 'DisplayName','MOSAIC'); hold on;

    plot(P(i).time_step * del_t_gtc, ...
        P(i).R * spdata.rmaj/100, ...
        'b--', 'LineWidth', 2, 'DisplayName','GTC');

    xlabel('$t*\Omega_{ci0}$','FontSize',xylabelsize,'Interpreter','latex');
    ylabel('$X$(m)','FontSize',xylabelsize,'Interpreter','latex');
    xlim([0,max(P(i).time_step * del_t_gtc)])
    grid on;

    %% ================== 3. Z vs Time ==================
    subplot(2,4,3+offset)
    %text(0.02,0.05,labels{3+offset}, ...
    %    'Units','normalized','FontSize',17);hold on
    plot(t_A, A_Z, 'r-', 'LineWidth', 2, 'DisplayName','MOSAIC'); hold on;

    plot(P(i).time_step * del_t_gtc, ...
        P(i).Z * spdata.rmaj/100, ...
        'b--', 'LineWidth', 2, 'DisplayName','GTC');

    xlabel('$t*\Omega_{ci0}$','FontSize',xylabelsize,'Interpreter','latex');
    ylabel('$Z$(m)','FontSize',xylabelsize,'Interpreter','latex');
    xlim([0,max(P(i).time_step * del_t_gtc)])
    ylim([-0.61,0.61])
    grid on;

    %% ================== 4. Zeta vs Time ==================
    subplot(2,4,4+offset)
    %text(0.02,0.05,labels{4+offset}, ...
    %    'Units','normalized','FontSize',17);hold on

    plot(t_A, A_zeta/pi, 'r-', 'LineWidth', 2, 'DisplayName','MOSAIC'); hold on;

    plot(P(i).time_step * del_t_gtc, ...
        P(i).zeta/pi, ...
        'b--', 'LineWidth', 2, 'DisplayName','GTC');

    xlabel('$t*\Omega_{ci0}$','FontSize',xylabelsize,'Interpreter','latex');
    ylabel('$\zeta/\pi$','FontSize',xylabelsize,'Interpreter','latex');
    xlim([0,max(P(i).time_step * del_t_gtc)])
    
    grid on;

end
saveas(gca,'../output/gtc_benchmark','png');
print(gcf, '-depsc', '-r600', ['../output/','gtc_benchmark','.eps']);
