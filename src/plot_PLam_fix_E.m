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
iplot = 3;

%n = 4;
%omega0 = 82.03*1000*2*pi;

if iplot == 1

    figure;
    subplot(141)
    omega_b_all = char_freq.omega_b;
    omega_b_all(find(strcmp(char_freq.type,'counter-passing'))) = nan;
    contourf(char_freq.Pzeta_norm,char_freq.lambda,omega_b_all,50);title('\omega\_b');colorbar;hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k--','linewidth',2);hold on;
    xlabel('P_\zeta/\psi_w');
    ylabel('\lambda');
    subplot(142)
    omega_phi_trap = char_freq.omega_phi;
    omega_phi_trap(find(~(strcmp(char_freq.type,'trapped')|strcmp(char_freq.type,'potato')))) = nan;
    contourf(char_freq.Pzeta_norm,char_freq.lambda,omega_phi_trap,50);title('\omega\_d');colorbar;hold on;
    contour(char_freq.Pzeta_norm,char_freq.lambda,char_freq.omega_d,[0 0],'linewidth',3,'color','c');hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k--','linewidth',2);hold on;
    xlabel('P_\zeta/\psi_w');
    ylabel('\lambda');
    subplot(143)
    omega_phi_copassing = char_freq.omega_phi;
    omega_phi_copassing(find(~strcmp(char_freq.type,'co-passing'))) = nan;
    contourf(char_freq.Pzeta_norm,char_freq.lambda,omega_phi_copassing,50);title('\omega\_\phi');colorbar;hold on;
    contour(char_freq.Pzeta_norm,char_freq.lambda,char_freq.omega_phi,[0 0],'linewidth',3,'color','c');hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k--','linewidth',2);hold on;
    xlabel('P_\zeta/\psi_w');
    ylabel('\lambda');
    subplot(144)
    contourf(char_freq.Pzeta_norm,char_freq.lambda,char_freq.q_avrg,50);title('q\_avrg');colorbar;hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k--','linewidth',2);hold on;
    xlabel('P_\zeta/\psi_w');
    ylabel('\lambda');

elseif iplot == 2
    colors = {'parula','jet','hot','turbo'}; % 不同 colormap
    fields = {'omega_b','omega_d','omega_phi','q_avrg'};
    figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.25,0.4],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.015]);

    l_res = (omega0 - n*char_freq.omega_phi)./char_freq.omega_b; % Eq. (70) in Bao et al 24 NF

    %l_res(find(strcmp(char_freq.type,'trapped'))) = nan;
    %l_res(find(strcmp(char_freq.type,'potato'))) = nan;
    %l_res(find(strcmp(char_freq.type,'stagnation'))) = nan;
    %l_res(find(strcmp(char_freq.type,'co-passing'))) = nan;
    %l_res(find(strcmp(char_freq.type,'counter-passing'))) = nan;
    l_res(find(strcmp(char_freq.exist,'no'))) = nan;
    l_res(find(strcmp(char_freq.confine,'lost'))) = nan;

    if strcmp(passing_option,'counter-passing')
        v = 12:1:23
    elseif strcmp(passing_option,'co-passing')
        v = -2:1:14
    end

    %v = 12:1:24
    contourf(char_freq.Pzeta_norm,char_freq.lambda,l_res,100,'edgecolor','none');colorbar('northoutside');colormap('jet');hold on;caxis([v(1) v(end)])
    %title('l\_res\_d');
    l_passing = l_res;
    v = 7:1:19;
    l_passing(find(~(strcmp(char_freq.type,'co-passing')|strcmp(char_freq.type,'counter-passing')))) = nan;
    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,l_passing,v,'facecolor','none','linewidth',1,'linecolor','k');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',400);
    set(h,'LineWidth',1.5,'linecolor','k','linestyle',':')
    v = -0:1:4;
    l_trap_potato = l_res;
    l_trap_potato(find(~(strcmp(char_freq.type,'trapped')|strcmp(char_freq.type,'potato')))) = nan;
    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,l_trap_potato,v,'facecolor','none','linewidth',1,'linecolor','k');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',400);
    set(h,'LineWidth',1.5,'linecolor','k','linestyle',':')

    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'color',[1,0.7,0],'linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'m','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    plot(PLam_2D.Pzeta_bound_right,PLam_2D.lam_uni_grid,'c-.','linewidth',2);hold on;
    %q=3等值线
    [C,h] =contour(char_freq.Pzeta_norm, char_freq.lambda, ...
        char_freq.(fields{4}), [0 -3], 'LineWidth', 2, 'Color','k');
    h.LineStyle = '--';
    %grid on;
    % axis([round(min(PLam_2D.Pzeta(1,:)))    ceil(max(PLam_2D.Pzeta_bound_right))    0    1.4000]);
    set(gca,'LineWidth',2,'FontSize',20,'ticklength',[0.02 0.0])
    xlabel('P_\zeta/\psi_w','fontsize',25)
    ylabel('\lambda=\muB_a/E','fontsize',25)
    strings = ['E = ',num2str(PLam_2D.energy_out),'keV'];
    annotation('textbox',[0.64 0.79 0.9 0.03],'String',strings,'fontsize',25,'edgecolor','none')

    %---------------------------------------------------------------------------------------------图例
    if strcmp(passing_option,'counter-passing')
        posx = min(xlim) + 0.05*range(xlim);
        posy = min(ylim) + 0.77*range(ylim);
    elseif strcmp(passing_option,'co-passing')
        posx = min(xlim) + 0.05*range(xlim);
        posy = min(ylim) + 0.12*range(ylim);
    end

    w = 0.45*range(xlim);
    h = 0.19*range(ylim);
    fill([posx posx+w posx+w posx],[posy-h/2 posy-h/2 posy+h posy+h],'w',...
        'EdgeColor','k','LineWidth',1.0);

    % 小椭圆渐变
    x0 = posx + 0.08*range(xlim);
    y0 = posy + 0.15*range(ylim);
    a = 0.03 * range(xlim);
    b = 0.03 * range(ylim);
    N = 200; theta = linspace(0,2*pi,N);
    cmap = colormap(colors{1}); nC = size(cmap,1);
    for k = 1:nC
        r = 1 - (k-1)/(nC-1);
        fill(x0 + a*r*cos(theta), y0 + b*r*sin(theta), cmap(k,:), ...
            'EdgeColor','none');
    end
    plot(x0 + a*cos(theta), y0 + b*sin(theta), 'k', 'LineWidth', 1.2);

    % 图例文字
    text(x0 + 0.06*range(xlim), y0, '\it{l}', ...
        'FontSize',23,'FontName','Times New Roman', ...
        'FontWeight','bold','Interpreter','tex', ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
    text(x0 - 0.045*range(xlim), y0-h/1.9, ' \cdot\cdot\cdot \it{l} \in ℤ', ...
        'FontSize',23,'FontName','Times New Roman', ...
        'FontWeight','bold','Interpreter','tex', ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
    text(x0 - 0.05*range(xlim), y0-h/1, '  --- \langle q \rangle = 3', ...
        'FontSize',23,'FontName','Times New Roman', ...
        'FontWeight','bold','Interpreter','tex', ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
    %--------------------------------------------------------------------------------------------------
    if strcmp(passing_option,'counter-passing')
        saveas(gca,'../output/pzeta_lambda_2d_phase_space_counter','epsc');
    elseif strcmp(passing_option,'co-passing')
        saveas(gca,'../output/pzeta_lambda_2d_phase_space_co','epsc');
    end

elseif iplot==3
    figure('Position',[100 100 1400 400],'DefaultAxesLineWidth',1);
% figure('name','region',...
%         'unit','normalized',...
%         'position',[100 100 1400 400],... % figure position
%         'DefaultAxesFontSize',20,...
%         'DefaultAxesFontWeight','normal',...
%         'DefaultAxesLineWidth',2,...
%         'DefaultAxesTickLength',[0.02,0.015]);
    titles = {'$\omega_\theta$','$\omega_d$','$\omega_\zeta$','$\langle q \rangle$'};
    fields = {'omega_b','omega_d','omega_phi','q_avrg'};
    colors = {'parula','jet','hot','turbo'}; % 不同 colormap
    n = numel(fields);

    for i = 1:n
        subplot(1,n,i)
        omega_va=4.3889e+06/sqrt(2);
        % ===== 主图绘制 =====
        if i==4
            contourf(char_freq.Pzeta_norm, char_freq.lambda, ...
                -char_freq.(fields{i}), 30, 'LineColor','none');hold on
            delta_fre=max(-char_freq.(fields{i})(:))-min(-char_freq.(fields{i})(:));
            caxis([min(-char_freq.(fields{i})(:))-0.1*delta_fre max(-char_freq.(fields{i})(:))]);
        else
            if i==1 && strcmp(passing_option,'co-passing')
                omega_ff = char_freq.omega_b;
                omega_ff(find(strcmp(char_freq.type,'counter-passing'))) = nan;
            elseif i==2 && strcmp(passing_option,'co-passing')
                    omega_ff = char_freq.omega_phi;
                    omega_ff(find(~(strcmp(char_freq.type,'trapped')|strcmp(char_freq.type,'potato')))) = nan;

            elseif i==3 && strcmp(passing_option,'co-passing')
                    omega_ff = char_freq.omega_phi;
                    omega_ff(find(~strcmp(char_freq.type,'co-passing'))) = nan;
            elseif strcmp(passing_option,'counter-passing')
                omega_ff =char_freq.(fields{i});
            end
            contourf(char_freq.Pzeta_norm, char_freq.lambda,omega_ff, 30, 'LineColor','none');hold on
            %delta_fre=max(char_freq.(fields{i})(:)/omega_va)-min(char_freq.(fields{i})(:)/omega_va);
            %caxis([min(omega_ff) max(omega_ff)]);
            end
            % contour(char_freq.Pzeta_norm, char_freq.lambda, ...
            %         char_freq.(fields{i}), [0 0],'linewidth',3,'color','c');hold on
            %colormap(colors{i});
            set(gca,'FontSize',11);
            axis tight;
            xlabel('$P_\zeta / \psi_w$','FontSize',15,'Interpreter','latex');
        ylabel('$\lambda$','FontSize',15,'Interpreter','latex');
        %title(titles{i},'Interpreter','tex','FontSize',14);

        % 上方 colorbar
        cb = colorbar('northoutside');
        %cb.Label.String = titles{i};
        cb.FontSize = 10;
        cb.Label.FontName = 'Times New Roman';
cb.FontWeight = 'bold';          
cb.LineWidth = 1; 
        hold on;
        % 轨迹线
        plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
        plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
        plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
        plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'color',[1,0.7,0],'linewidth',2);hold on;
        plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'m','linewidth',2);hold on;
        %plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
        plot(PLam_2D.Pzeta_bound_right,PLam_2D.lam_uni_grid,'c-.','linewidth',2);hold on;
        if strcmp(passing_option,'counter-passing')
            posx = min(xlim) + 0.47*range(xlim);
            posy = min(ylim) + 0.87*range(ylim);
            xlim([-1.32,0.0]);
        elseif strcmp(passing_option,'co-passing')
            posx = min(xlim) + 0.325*range(xlim);
            posy = min(ylim) + 0.08*range(ylim);
            xlim([-1,0.537]);
        end

        % q=3 等值线
        if i <= 4
            [C,h] =contour(char_freq.Pzeta_norm, char_freq.lambda, ...
                char_freq.(fields{4}), [0 -3], 'LineWidth', 1.5, 'Color','k');
            h.LineStyle = '--';
        end

        % ====== 添加椭圆渐变色小图例 ======
        % ===== 图例框部分 =====
        % 白底矩形框
        %---------------------------------------------------------------------------------------------图例
        if i<= 3
            %posx = min(xlim) + 0.05*range(xlim);
            %posy = min(ylim) + 0.07*range(ylim);
            w = 0.46*range(xlim);
            h = 0.12*range(ylim);
            fill([posx posx+w posx+w posx],[posy-h/2 posy-h/2 posy+h posy+h],'w',...
                'EdgeColor','k','LineWidth',1.0);

            % 小椭圆渐变
            x0 = posx + 0.06*range(xlim);
            y0 = posy + 0.07*range(ylim);
            a = 0.04 * range(xlim);
            b = 0.02 * range(ylim);
            N = 200; theta = linspace(0,2*pi,N);
            cmap = colormap(colors{2}); nC = size(cmap,1);
            for k = 1:nC
                r = 1 - (k-1)/(nC-1);
                fill(x0 + a*r*cos(theta), y0 + b*r*sin(theta), cmap(k,:), ...
                    'EdgeColor','none');
            end
            plot(x0 + a*cos(theta), y0 + b*sin(theta), 'k', 'LineWidth', 1.2);

            % 图例文字
            text(x0 + 0.06*range(xlim), y0, ['  ' titles{i} '{(rad/s)}'], ...
                'FontSize',13,'FontName','Times New Roman', ...
                'FontWeight','normal','Interpreter','latex', ...
                'HorizontalAlignment','left','VerticalAlignment','middle');
            text(x0 - 0.03*range(xlim), y0-h/1.5, '- -$\langle q \rangle$ = 3', ...
                'FontSize',13,'FontName','Times New Roman', ...
                'FontWeight','normal','Interpreter','latex', ...
                'HorizontalAlignment','left','VerticalAlignment','middle');
        else
            % posx = min(xlim) + 0.05*range(xlim);
            %posy = min(ylim) + 0.07*range(ylim);
            w = 0.46*range(xlim);
            h = 0.12*range(ylim);
            fill([posx posx+w posx+w posx],[posy-h/2 posy-h/2 posy+h posy+h],'w',...
                'EdgeColor','k','LineWidth',1.0);

            % 小椭圆渐变
            x0 = posx + 0.06*range(xlim);
            y0 = posy + 0.07*range(ylim);
            a = 0.04 * range(xlim);
            b = 0.02 * range(ylim);
            N = 200; theta = linspace(0,2*pi,N);
            cmap = colormap(colors{2}); nC = size(cmap,1);
            for k = 1:nC
                r = 1 - (k-1)/(nC-1);
                fill(x0 + a*r*cos(theta), y0 + b*r*sin(theta), cmap(k,:), ...
                    'EdgeColor','none');
            end
            plot(x0 + a*cos(theta), y0 + b*sin(theta), 'k', 'LineWidth', 1.2);

            % 图例文字
            text(x0 + 0.05*range(xlim), y0, [titles{i}], ...
                'FontSize',13,'FontName','Times New Roman', ...
                'FontWeight','normal','Interpreter','latex', ...
                'HorizontalAlignment','left','VerticalAlignment','middle');
            text(x0 - 0.03*range(xlim), y0-h/1.5, '- -$\langle q \rangle$ = 3', ...
                'FontSize',13,'FontName','Times New Roman', ...
                'FontWeight','normal','Interpreter','latex', ...
                'HorizontalAlignment','left','VerticalAlignment','middle');
        end
        %---------------------------------------------------------------------------------------------图例
    end
    if strcmp(passing_option,'counter-passing')
        saveas(gca,'../output/pzeta_lambda_freq_counter','epsc');
    elseif strcmp(passing_option,'co-passing')
        saveas(gca,'../output/pzeta_lambda_freq_co','epsc');
    end
elseif iplot==4
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
end