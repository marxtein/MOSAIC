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
global particle

%iplot = 2;

if iplot == 1 % l_num*P_num aligned particle trajectory
    %%
    figure('name','orbit',...
        'unit','normalized',...
        'position',[0.0,0.0,0.68,1],... % figure position
        'DefaultAxesFontSize',17,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',1,...
        'DefaultAxesTickLength',[0.01,0.01]);
    
    for m=1:np
        p(m)=subplot(l_num,P_num,m);
        %set(p(m),'position',[0.04+mod((m-1),5)*0.193    0.825-floor((m-1)/5)*0.193    0.15    0.15]);
set(p(m),'position',[0.05+mod((m-1),5)*0.155   0.83-floor((m-1)/5)*0.195   0.17   0.15]);
    
        plot(squeeze(eq.xsp(1,end,:)*spdata.rmaj/100),squeeze(eq.zsp(1,end,:)*spdata.rmaj/100),'k-','linewidth',2);hold on;
        plot(eq.xsp(1,1,1)*spdata.rmaj/100,0,'k.','markersize',15);hold on;
        daspect([1 1 1]);
        xlim([0.45 1.45]*spdata.rmaj/100);
        
        if mod(m,P_num)==1
            ylabel('$Z$(m)','Interpreter','latex')
        end
        if m>np-P_num
            xlabel('$X$(m)','Interpreter','latex')
        end
        set(get(gca,'Title'),'Position',[1.6, 1.03, 0]);
        if strcmp(char_freq.exist{m},'no')
            title(['$l$=',num2str(m),' (',char_freq.type{m},')'],'Interpreter','latex');
            
            continue;
        end
        plot(A.R(:,m)*spdata.rmaj/100,A.Z(:,m)*spdata.rmaj/100,'r.');hold on;
        if m<np
            [~,ind] = min(abs(A.init.psi(m) - eq.psi));
        end
        plot(squeeze(eq.xsp(1,ind,:)*spdata.rmaj/100),squeeze(eq.zsp(1,ind,:)*spdata.rmaj/100),'k--','linewidth',2);hold on;
        plot(squeeze(eq.xsp(1,end,:)*spdata.rmaj/100),squeeze(eq.zsp(1,end,:)*spdata.rmaj/100),'k-','linewidth',2);hold on;
        %plot(eq.xsp(1,1,1)*spdata.rmaj/100,eq.zsp(1,1,1)*spdata.rmaj,'*','linewidth',1.5,'markersize',10);hold on;
        %daspect([1 1 1]);
        %xlim([0.3 1.5]*spdata.rmaj/100);
        %set(gca,'xtick',(0.6:0.4:1.4)*spdata.rmaj/100)
        if strcmp(char_freq.confine{m},'lost')
            title(['$l$=',num2str(m),' (',char_freq.confine{m},')'],'Interpreter','latex');
            
        else
            title(['$l$=',num2str(m),' (',char_freq.type{m},')'],'Interpreter','latex');
            
        end
        %ax = subplot(l_num,P_num,np); % Last subplot
hold on;


        set(gca,'Fontsize',16)
    end
  % ==========================

% ==========================
ax = axes('Position',[0 0 0.995 1.85],'Visible','off','HandleVisibility','off');
hold(ax,'on');

% Draw three invisible lines, specifically for the Legend
h1 = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2);                 % Particle trajectory
h2 = plot(ax, NaN, NaN, 'k--', 'LineWidth', 2);               % Initial flux surface
h3 = plot(ax, NaN, NaN, 'k-', 'LineWidth', 2);                % Outermost flux surface
h4 = plot(ax, NaN, NaN, 'k.', 'MarkerSize', 15);             % Magnetic axis (black dot)

% Full Legend
legend(ax, [h1,h2,h3,h4], ...
    'Particle trajectory', ...
    'Initial flux surface', ...
    'Outermost flux surface', ...
    'Magnetic axis', ...
    'Location','eastoutside',...
    'FontSize',20,...
    'Interpreter','latex');

hold(ax,'off');
    if ~exist('../output','dir');
        mkdir('../output');
    end
    saveas(gca,'../output/orbit','png');
    print(gcf, '-depsc', '-r600', ['../output/','orbit','.eps']);

    if strcmp(ant_option,'on')
        saveas(gca,'../output/orbit_ant','png');
        print(gcf, '-depsc', '-r600', ['../output/','orbit_ant','.eps']);
    else
        saveas(gca,'../output/orbit','png');
        print(gcf, '-depsc', '-r600', ['../output/','orbit','.eps']);
    end
    %%
    % figure('name','energy',...
    %     'unit','normalized',...
    %     'position',[0.0,0.0,0.6,1],... % figure position
    %     'DefaultAxesFontSize',15,...
    %     'DefaultAxesFontWeight','normal',...
    %     'DefaultAxesLineWidth',1,...
    %     'DefaultAxesTickLength',[0.02,0.015]);
    % 
    % for m=1:np
    %     p(m)=subplot(l_num,P_num,m);
    %     set(p(m),'position',[0.04+mod((m-1),5)*0.196    0.825-floor((m-1)/5)*0.196    0.15    0.15]);
    %     if strcmp(char_freq.exist{m},'no')
    %         title(['$l$=',num2str(m),' (',char_freq.type{m},')'],'Interpreter','latex');
    %         if (m-1)/5 >= 4
    %             xlabel('time step')
    %         end
    %         continue;
    %     end
    %     if strcmp(ant_option,'on')
    %         plot(A.E_motion(:,m),'r-','linewidth',2);hold on;
    %         plot(A.E_power(:,m),'k--','linewidth',1);hold on;
    %         %plot(A.E_prime(:,m),'c--','linewidth',1);hold on;
    %         legend('$E$ (Equation of motion)','$E$ (Power transfer)','Interpreter','latex')
    %     else
    %         plot(A.E_motion(:,m),'r-','linewidth',2);hold on;
    %         legend('E (Equation of motion)','Interpreter','latex')
    %     end
    %     if strcmp(char_freq.confine{m},'lost')
    %         title(['m=',num2str(m),' (',char_freq.confine{m},')']);
    %     else
    %         title(['m=',num2str(m),' (',char_freq.type{m},')']);
    %     end
    %     if (m-1)/5 >= 4
    %         xlabel('time step')
    %     end
    % end
    % if ~exist('../output','dir');
    %     mkdir('../output');
    % end
    % saveas(gca,'../output/energy','png');
    % print(gcf, '-depsc', '-r600', ['../output/','energy','.eps']);
 %%    
    
    hFig = figure ;
    set(hFig, 'Position',[1,1,2080,350]);
    mm=0;
    for m=[11 12 19 23 24]
        mm=mm+1;
        subplot(1,5,mm);
         if strcmp(char_freq.exist{m},'no')
            title(['m=',num2str(m),' (',char_freq.type{m},')']);
          
            continue;
        end
        if strcmp(ant_option,'on')
            plot((1:size(A.R,1))*tstep,A.E_motion(:,m),'r-','linewidth',2);hold on;
            plot((1:size(A.R,1))*tstep,A.E_power(:,m),'b--','linewidth',2);hold on;
            legend('$E$ from Eq.(38)','$E$ from Eq.(40)', 'Location', 'southeast','Interpreter','latex')
            %legend('E (Equation of motion)','E (Equation of power transfer)', 'Location', 'southeast')
        else
            plot((1:size(A.R,1))*tstep,A.E_motion(:,m),'r-','linewidth',2);hold on;
            legend('$E$ (Equation of motion)', 'Location', 'southeast','Interpreter','latex')
        end
        if strcmp(char_freq.confine{m},'lost')
            title(['$l$=',num2str(m),' (',char_freq.confine{m},')'],'Interpreter','latex');
        else
            title(['$l$=',num2str(m),' (',char_freq.type{m},')'],'Interpreter','latex');
        end
        grid on;
        set(gca,'Fontsize',16)
        if (m-1)/5 >= 4
            xlabel('$t*\Omega_{ci0}$','Interpreter','latex','Fontsize',20)
        end
        xlabel('$t*\Omega_{ci0}$','Interpreter','latex','Fontsize',20)
        if m==11
                ylabel('$E$(keV)','Interpreter','latex','Fontsize',20)
        end
        
    end
    if strcmp(ant_option,'on')
        saveas(gca,'../output/energy_paper_ant','epsc');
        saveas(gca,'../output/energy_paper_ant','png');
    else
        saveas(gca,'../output/energy_paper_ant','epsc');
        saveas(gca,'../output/energy_paper_ant','png');
    end

    % %%
    % figure('name','Pzeta',...
    %     'unit','normalized',...
    %     'position',[0.0,0.0,0.6,1],... % figure position
    %     'DefaultAxesFontSize',15,...
    %     'DefaultAxesFontWeight','normal',...
    %     'DefaultAxesLineWidth',0.5,...
    %     'DefaultAxesTickLength',[0.02,0.015]);
    % 
    % for m=1:np
    %     p(m)=subplot(l_num,P_num,m);
    %     set(p(m),'position',[0.04+mod((m-1),5)*0.196    0.825-floor((m-1)/5)*0.196    0.15    0.15]);
    %     if strcmp(char_freq.exist{m},'no')
    %         title(['m=',num2str(m),' (',char_freq.type{m},')']);
    %         if (m-1)/5 >= 4
    %             xlabel('time step')
    %         end
    %         continue;
    %     end
    % 
    %     plot(A.Pzeta_norm(:,m),'b-','linewidth',2);hold on;
    %     legend({'$P_{\zeta}/\psi_w$'},'interpreter','latex','fontsize',15)
    % 
    %     if strcmp(char_freq.confine{m},'lost')
    %         title(['m=',num2str(m),' (',char_freq.confine{m},')']);
    %     else
    %         title(['m=',num2str(m),' (',char_freq.type{m},')']);
    %     end
    %     if (m-1)/5 >= 4
    %         xlabel('time step')
    %     end
    % end
    % if ~exist('../output','dir');
    %     mkdir('../output');
    % end
    % saveas(gca,'../output/Pzeta','png');
    % print(gcf, '-depsc', '-r600', ['../output/','Pzeta','.eps']);
    
    
elseif iplot == 2 % Pzeta-lambda phase space markered with particle trajectory classification
    %%
    figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.45,0.59],... % figure position
        'DefaultAxesFontSize',25,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',1,...
        'DefaultAxesTickLength',[0.005,0.005]);
    
    % plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta_bound_right,PLam_2D.lam_uni_grid,'m--','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'color',[1,0.7,0],'linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'m','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    plot(PLam_2D.Pzeta_bound_right,PLam_2D.lam_uni_grid,'c-.','linewidth',2);hold on;

   
    run trapped_region.m
    h3 = fill(x_fill, y_fill, [1,1,0]);  
    set(h3,'facealpha',0.3,'edgecolor','none');

    run copass_region1.m
    h4 = fill(x_fill, y_fill, 'r'); 
    set(h4,'facealpha',0.15,'edgecolor','none');

    run copass_region2.m
    h5 = fill(x_fill, y_fill, 'r');  
    set(h5,'facealpha',0.15,'edgecolor','none');

    run potato_region.m
    h6 = fill(x_fill, y_fill, [0,0,1]);  
    set(h6,'facealpha',0.3,'edgecolor','none');

    run stagnation_region.m
    h7 = fill(x_fill, y_fill, [0,1,0]);  
    set(h7,'facealpha',0.3,'edgecolor','none');

    
    h = fill([PLam_2D.Pzeta_bound_right fliplr(PLam_2D.Pzeta_bound_left1)],[PLam_2D.lam_uni_grid fliplr(PLam_2D.lam_uni_grid)],'r');
    set(h,'edgealpha',0,'facealpha',0.0)
    h2 = fill([PLam_2D.Pzeta_bound_counter_passing_left fliplr(PLam_2D.Pzeta_bound_counter_passing_right)],[PLam_2D.lam_uni_grid_counter_passing fliplr(PLam_2D.lam_uni_grid_counter_passing)],'c');
    set(h2,'edgealpha',0,'facealpha',0.2)
    
    fsize=15;%'FontWeight', 'bold'
   % annotation('textbox',[0.7 0.74 0.9 0.03],'String','$CP$','fontsize',fsize,'edgecolor','none','interpreter','latex')
    % text(-1, 0.3, '$Ctr-P$', 'FontSize', 20, 'Color', 'k','interpreter','latex');
    % text(-0.19, 0.5, '$Co-P$', 'FontSize', 20, 'Color', 'k','interpreter','latex');
    % text(-0.4, 1.05, '$T$', 'FontSize', 20, 'Color', 'k','interpreter','latex');
    % text(-0.15, 0.9, '$P$', 'FontSize', 20, 'Color', 'k','interpreter','latex');
    % text(-0.0, 1.03, '$S$', 'FontSize', 20, 'Color', 'k','interpreter','latex');
    
    grid on;
    
    l_num = 5;
    P_num = 5;
    lam_tmp = linspace(traj.lam0,traj.lam1,l_num);
    Pzeta_norm_tmp = linspace(traj.Pzeta_norm0,traj.Pzeta_norm1,P_num);
    [Pzeta_2D,lam_2D] = meshgrid(Pzeta_norm_tmp,lam_tmp);
    
    
    PLam_2D.Pzeta_1D = reshape(Pzeta_2D,l_num*P_num,1);
    PLam_2D.lam_1D = reshape(lam_2D,l_num*P_num,1);
    PLam_2D.E_1D = PLam_2D.energy_out*ones(l_num*P_num,1);
    h_confined = [];
    h_lost = [];
    for i = 1:l_num*P_num
        if strcmp(char_freq.confine{i},'confined')
            h=plot(PLam_2D.Pzeta_1D(i),PLam_2D.lam_1D(i),'mo','markersize',12,'markerfacecolor','m');hold on;
            h_confined = [h_confined, h];
        elseif strcmp(char_freq.confine{i},'lost')
            h=plot(PLam_2D.Pzeta_1D(i),PLam_2D.lam_1D(i),'mo','markersize',12);hold on;
            h_lost = [h_lost, h];
        elseif strcmp(char_freq.confine{i},'non-defined')
            plot(PLam_2D.Pzeta_1D(i),PLam_2D.lam_1D(i),'ko','markersize',12);hold on;
        end
        text(PLam_2D.Pzeta_1D(i),PLam_2D.lam_1D(i),num2str(i),'Color',[0.9290 0.6940 0.1250],'FontSize',24,'FontWeight','bold',...
            'HorizontalAlignment','center','VerticalAlignment','top');
    end
    legend_handles = [h_confined(1), h_lost(1)];
    legend_labels = {'Confined particles', 'Lost particles'};
    set(gca,'ytick',0:0.2:1.6)
    grid on;
    axis([floor(min(PLam_2D.Pzeta_bound_left2))    ceil(max(PLam_2D.Pzeta_bound_right))    0    1.5000]);
    xlabel('$P_\zeta/\psi_w$','fontsize',35,'interpreter','latex')
    ylabel('$\lambda=\mu B_a/E$','fontsize',35,'interpreter','latex')
    %strings = ['$E$=',num2str(energy_in),'keV',',','$\ m/m_p$=',num2str(particle.mass),',','$\ Z/e$=',num2str(particle.charge)];
    %title(strings,'fontsize',40,'Interpreter','latex');

    %strings = ['$E$=',num2str(energy_in),'keV'];
    %annotation('textbox',[0.35 0.84 0.9 0.03],'String',strings,'fontsize',35,'edgecolor','none','Interpreter','latex')

    % strings = ['$\ m/m_p$=',num2str(particle.mass)];
    % annotation('textbox',[0.6 0.81 0.9 0.03],'String',strings,'fontsize',35,'edgecolor','none','Interpreter','latex')
    % strings = ['$\ Z/e$=',num2str(particle.charge)];
    % annotation('textbox',[0.53 0.74 0.9 0.03],'String',strings,'fontsize',35,'edgecolor','none','Interpreter','latex')
    %hl = legend([h h2],'Co-passing/Trapped/Potato/Stagnation','Counter-passing');
     hl = legend([h2 h4 h3 h6 h7 legend_handles],'Counter-passing','Co-passing','Trapped','Potato','Stagnation','Confined particles', ...
    'Lost particles');
   
    set(hl,'fontsize',19)
    
    set(gca,'ytick',0:0.2:1.6)
    ylim([0,1.6]);
    if ~exist('../output','dir');
        mkdir('../output');
    end
    saveas(gca,'../output/CoM_2','png');
    print(gcf, '-depsc', '-r600', ['../output/','CoM_2','.eps']);
    
elseif iplot == 3
    %%
    figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.4,0.5],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.015]);
    
    % plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta_bound_right,PLam_2D.lam_uni_grid,'m--','linewidth',2);hold on;

    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'color',[1,0.7,0],'linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'m','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    plot(PLam_2D.Pzeta_bound_right,PLam_2D.lam_uni_grid,'c-.','linewidth',2);hold on;

%test

% figure
% plot(xA_seg,yA_seg);hold on;
% plot(xB_seg,yB_seg);hold on;
% plot(xC_seg,yC_seg);hold on;
% %


    run trapped_region.m
    h3 = fill(x_fill, y_fill, [1,1,0]);  % Yellow
    set(h3,'facealpha',0.0,'edgecolor','none');

    run copass_region1.m
    h4 = fill(x_fill, y_fill, [0,1,0]);  % Yellow
    set(h4,'facealpha',0.0,'edgecolor','none');

    run copass_region2.m
    h5 = fill(x_fill, y_fill, [0,1,0]);  % Yellow
    set(h5,'facealpha',0.0,'edgecolor','none');

    run potato_region.m
    h6 = fill(x_fill, y_fill, [0,0,1]);  % Yellow
    set(h6,'facealpha',0.0,'edgecolor','none');

    run stagnation_region.m
    h7 = fill(x_fill, y_fill, [0,1,1]);  % Yellow
    set(h7,'facealpha',0.0,'edgecolor','none');

    
    h = fill([PLam_2D.Pzeta_bound_right fliplr(PLam_2D.Pzeta_bound_left1)],[PLam_2D.lam_uni_grid fliplr(PLam_2D.lam_uni_grid)],'r');
    set(h,'edgealpha',0,'facealpha',0.0)
    h2 = fill([PLam_2D.Pzeta_bound_counter_passing_left fliplr(PLam_2D.Pzeta_bound_counter_passing_right)],[PLam_2D.lam_uni_grid_counter_passing fliplr(PLam_2D.lam_uni_grid_counter_passing)],'c');
    set(h2,'edgealpha',0,'facealpha',1)
    hl = legend([h2 h4 h3 h6 h7],'Counter-passing','Co-passing','Trapped','Potato','Stagnation');
    set(hl,'fontsize',25)
    set(gca,'ytick',0:0.2:1.6)
    grid on;
    axis([floor(min(PLam_2D.Pzeta_bound_left2))    ceil(max(PLam_2D.Pzeta_bound_right))    0    1.5000]);
    xlabel('$P_\zeta/\psi_w$','fontsize',30,'interpreter','latex')
    ylabel('$\lambda=\mu B_a/E$','fontsize',30,'interpreter','latex')
    strings = ['$E$=',num2str(energy_in),'keV',',','$\ m/m_p$=',num2str(particle.mass),',','$\ Z/e$=',num2str(particle.charge)];
    title(strings,'fontsize',40,'Interpreter','latex');
    
    if ~exist('../output','dir');
        mkdir('../output');
    end
    saveas(gca,'../output/CoM','png');
    print(gcf, '-depsc', '-r600', ['../output/','CoM','.eps']);
    
end