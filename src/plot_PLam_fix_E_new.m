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
iplot = 1;

n = 1;
omega0 = 82.03*1000*2*pi; %EHL-2: TAE
q_mode = 3;

if iplot == 1
    %%
    figure('unit','normalized',...
        'Position',[0.0 0.0 0.8 0.5],...
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.02]);
    
    p1 = subplot(141)
    set(p1,'position',[0.07 0.13 0.18 0.8]);
    omega_b_all = char_freq.omega_b;
    omega_b_all(find(strcmp(char_freq.type,'counter-passing'))) = nan;
    contourf(char_freq.Pzeta_norm, char_freq.lambda, omega_b_all,50,'edgecolor','none');colormap('jet');hold on;
    cc=colorbar('northoutside')
    set(cc, 'linewidth',2,'Ticklength',[0.02,0.015],'TickDirection','out');
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    
    legend({'$\omega_\theta(rad/s)$'},'location','northeast','Interpreter','latex','fontsize',28);
    set(gca,'TickDir','out','box','off');
    xlabel('$P_\zeta/\psi_w$','Interpreter','latex','fontsize',30);
    ylabel('$\lambda=\mu B_a/E$','Interpreter','latex','fontsize',30);
    
    
    p2 = subplot(142)
    set(p2,'position',[0.3 0.13 0.18 0.8]);
    omega_phi_trap = char_freq.omega_phi;
    omega_phi_trap(find(~(strcmp(char_freq.type,'trapped')|strcmp(char_freq.type,'potato')))) = nan;
    contourf(char_freq.Pzeta_norm, char_freq.lambda, omega_phi_trap,50,'edgecolor','none');hold on;
    cc=colorbar('northoutside')
    set(cc, 'linewidth',2,'Ticklength',[0.02,0.015],'TickDirection','out');
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    
    legend({'$\omega_\phi(rad/s)$'},'location','northeast','Interpreter','latex','fontsize',28);
    set(gca,'TickDir','out','box','off');
    xlabel('$P_\zeta/\psi_w$','Interpreter','latex','fontsize',30);
    
    
    p3 = subplot(143)
    set(p3,'position',[0.53 0.13 0.18 0.8]);
    omega_phi_copassing = char_freq.omega_phi;
    omega_phi_copassing(find(~strcmp(char_freq.type,'co-passing'))) = nan;
    contourf(char_freq.Pzeta_norm, char_freq.lambda, omega_phi_copassing,50,'edgecolor','none');hold on;
    cc=colorbar('northoutside')
    set(cc, 'linewidth',2,'Ticklength',[0.02,0.015],'TickDirection','out');
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    
    legend({'$\omega_\phi(rad/s)$'},'location','northeast','Interpreter','latex','fontsize',28);
    set(gca,'TickDir','out','box','off');
    xlabel('$P_\zeta/\psi_w$','Interpreter','latex','fontsize',30);
    
    
    p4 = subplot(144)
    set(p4,'position',[0.76 0.13 0.18 0.8]);
    [C,h] = contourf(char_freq.Pzeta_norm, char_freq.lambda, char_freq.q_avrg,50,'edgecolor','none');hold on;
    [C,h2] = contourf(char_freq.Pzeta_norm, char_freq.lambda, char_freq.q_avrg,[q_mode,q_mode],'facecolor','none','linewidth',4,'linecolor','m','linestyle',':');hold on;
    cc=colorbar('northoutside')
    set(cc, 'linewidth',2,'Ticklength',[0.02,0.015],'TickDirection','out');
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    
    hl=legend([h,h2],'$\langle q\rangle$',['$\langle q\rangle=$',num2str(q_mode)]);
    set(hl,'Interpreter','latex','location','northeast','fontsize',28);
    set(gca,'TickDir','out','box','off');
    xlabel('$P_\zeta/\psi_w$','Interpreter','latex','fontsize',30);
    
    
    
    if ~exist('../output','dir');
        mkdir('../output');
    end
    saveas(gca,'../output/char_freq','png');
    print(gcf, '-depsc', '-r600', ['../output/','char_freq.eps']);
    
elseif iplot == 2
    %%
    figure('unit','normalized',...
        'Position',[0.0 0.0 0.6 0.5],...
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.02]);
    
    p1 = subplot(131)
    set(p1,'position',[0.08 0.13 0.25 0.8]);
    omega_b = char_freq.omega_b;
    omega_b(find(strcmp(char_freq.type,'co-passing'))) = nan;
    contourf(char_freq.Pzeta_norm, char_freq.lambda, omega_b,50,'edgecolor','none');colormap('jet');hold on;
    cc=colorbar('northoutside')
    set(cc, 'linewidth',2,'Ticklength',[0.02,0.015],'TickDirection','out');
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    
    legend({'$\omega_\theta(rad/s)$'},'location','northeast','Interpreter','latex','fontsize',28);
    set(gca,'TickDir','out','box','off');
    xlabel('$P_\zeta/\psi_w$','Interpreter','latex','fontsize',30);
    ylabel('$\lambda=\mu B_a/E$','Interpreter','latex','fontsize',30);
    
    
    p2 = subplot(132)
    set(p2,'position',[0.4 0.13 0.25 0.8]);
    omega_phi = char_freq.omega_phi;
    omega_phi(find(~strcmp(char_freq.type,'counter-passing'))) = nan;
    contourf(char_freq.Pzeta_norm, char_freq.lambda, omega_phi,50,'edgecolor','none');hold on;
    cc=colorbar('northoutside')
    set(cc, 'linewidth',2,'Ticklength',[0.02,0.015],'TickDirection','out');
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    
    legend({'$\omega_\phi(rad/s)$'},'location','northeast','Interpreter','latex','fontsize',28);
    set(gca,'TickDir','out','box','off');
    xlabel('$P_\zeta/\psi_w$','Interpreter','latex','fontsize',30);
    
    
    p3 = subplot(133)
    set(p3,'position',[0.72 0.13 0.25 0.8]);
    [C,h] = contourf(char_freq.Pzeta_norm, char_freq.lambda, char_freq.q_avrg,50,'edgecolor','none');hold on;
    [C,h2] = contourf(char_freq.Pzeta_norm, char_freq.lambda, char_freq.q_avrg,[q_mode,q_mode],'facecolor','none','linewidth',4,'linecolor','m','linestyle',':');hold on;
    cc=colorbar('northoutside')
    set(cc, 'linewidth',2,'Ticklength',[0.02,0.015],'TickDirection','out');
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    
    hl=legend([h,h2],'$\langle q\rangle$',['$\langle q\rangle=$',num2str(q_mode)]);
    set(hl,'Interpreter','latex','location','northeast','fontsize',28);
    set(gca,'TickDir','out','box','off');
    xlabel('$P_\zeta/\psi_w$','Interpreter','latex','fontsize',30);
    
    
    
    if ~exist('../output','dir');
        mkdir('../output');
    end
    saveas(gca,'../output/char_freq2','png');
    print(gcf, '-depsc', '-r600', ['../output/','char_freq2.eps']);
    
elseif iplot == 3
    %%
    figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.4,0.5],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.015]);
    
    p1 = subplot(121)
    set(p1,'position',[0.1 0.13 0.4 0.8]);
    p_harmonic = (n*char_freq.omega_phi - omega0)./char_freq.omega_b; % Eq. (70) in Bao et al 24 NF
    p_harmonic(find(strcmp(char_freq.type,'co-passing')|strcmp(char_freq.type,'stagnation'))) = nan;
    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,p_harmonic,100,'edgecolor','none');colorbar;colormap('jet');hold on;
    [C2,h2] = contourf(char_freq.Pzeta_norm,char_freq.lambda,p_harmonic,floor(min(min(p_harmonic))):ceil(max(max(p_harmonic))),'facecolor','none','linewidth',2,'linecolor','k');hold on;
    clabel(C2,h2,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',400);
    set(h2,'LineWidth',2,'linecolor','k','linestyle','-')
    
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    legend([h],{'p\_harmonic'},'location','southeast','fontsize',25);
    set(gca,'ytick',0:0.2:1.6)
    xlabel('$P_\zeta/\psi_w$','fontsize',30,'interpreter','latex')
    ylabel('$\lambda=\mu B_a/E$','fontsize',30,'interpreter','latex')
    title('Trapped');
    
    p2 = subplot(122)
    set(p2,'position',[0.57 0.13 0.4 0.8]);
    p_harmonic = (n*char_freq.omega_phi - omega0)./char_freq.omega_b; % omega = n*omega_phi - p*omega_theta
    p_harmonic(find(strcmp(char_freq.type,'trapped')|strcmp(char_freq.type,'potato'))) = nan;
    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,p_harmonic,100,'edgecolor','none');colorbar;colormap('jet');hold on;
    [C2,h2] = contourf(char_freq.Pzeta_norm,char_freq.lambda,p_harmonic,floor(min(min(p_harmonic))):ceil(max(max(p_harmonic))),'facecolor','none','linewidth',2,'linecolor','k');hold on;
    clabel(C2,h2,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',400);
    set(h2,'LineWidth',2,'linecolor','k','linestyle','-')
    
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    legend([h],{'p\_harmonic'},'location','northeast','fontsize',25);
    set(gca,'ytick',0:0.2:1.6)
    xlabel('$P_\zeta/\psi_w$','fontsize',30,'interpreter','latex')
    title('Co-passing');
    
    if ~exist('../output','dir');
        mkdir('../output');
    end
    saveas(gca,'../output/resonance','png');
    print(gcf, '-depsc', '-r600', ['../output/','resonance.eps']);
    
elseif iplot == 4
    %%
    figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.24,0.5],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.015]);
    
    
    set(gca,'position',[0.18 0.13 0.73 0.8]);
    p_harmonic = (n*char_freq.omega_phi - omega0)./char_freq.omega_b; % omega = n*omega_phi - p*omega_theta
    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,p_harmonic,100,'edgecolor','none');colorbar;colormap('jet');hold on;
    [C2,h2] = contourf(char_freq.Pzeta_norm,char_freq.lambda,p_harmonic,floor(min(min(p_harmonic))):ceil(max(max(p_harmonic))),'facecolor','none','linewidth',2,'linecolor','k');hold on;
    clabel(C2,h2,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',400);
    set(h2,'LineWidth',2,'linecolor','k','linestyle','-')
    
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    legend([h],{'p\_harmonic'},'location','northeast','fontsize',25);
    set(gca,'ytick',0:0.2:1.6)
    xlabel('$P_\zeta/\psi_w$','fontsize',30,'interpreter','latex');
    ylabel('$\lambda=\mu B_a/E$','fontsize',30,'interpreter','latex');
    title('Counter-passing');
    
    if ~exist('../output','dir');
        mkdir('../output');
    end
    saveas(gca,'../output/resonance2','png');
    print(gcf, '-depsc', '-r600', ['../output/','resonance2.eps']);
    
end