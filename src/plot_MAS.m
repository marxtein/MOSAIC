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
if strcmp(ps_option,'MG')
    iplot = 1;
elseif strcmp(ps_option,'GB')
    iplot = 2;
elseif strcmp(ps_option,'PLam')
    iplot = 3;
end

n = 4;
omega0 = 82.03*1000*2*pi;


if iplot == 1

    l_res = (omega0 - n*char_freq.omega_phi)./char_freq.omega_b;

    if strcmp(passing_option,'counter-passing')
        v = 12:1:23;
        %l_res = l_res1;
    elseif strcmp(passing_option,'co-passing')
        v = 4:1:14;
        %l_res = l_res;
        for pp=1:size(char_freq.Pzeta_norm,1)
            for ppp=1:size(char_freq.Pzeta_norm,2)
                fdp(pp,ppp)=min(min(((char_freq.Pzeta_norm(pp,ppp)-char_freq1.Pzeta_norm).^2+(char_freq.E(pp,ppp)-char_freq1.E).^2).^0.5));
            end
        end
    end

    %contourf(char_freq.Pzeta_norm,char_freq.E,l_res,100,'edgecolor','none');colormap('parula');hold on;caxis([v(1) v(end)])
    %cc=colorbar('northoutside');

    %contourf(char_freq.Pzeta_norm,char_freq.E,l_res-round(l_res),100,'edgecolor','none');title('l\_res\_d');colorbar;colormap('jet');hold on;caxis([-0.5 0.5])
    l_passing = l_res;
    l_passing(find(~(strcmp(char_freq.type,'co-passing')|strcmp(char_freq.type,'counter-passing')))) = nan;
    l_passing1 = l_passing;
    l_passing2 = l_passing;
    if strcmp(passing_option,'counter-passing')
        v = 12:1:17;
    elseif strcmp(passing_option,'co-passing')
        v = 8:1:14;
        l_passing1(find(fdp<0.015))=nan;
        l_passing2(find(fdp>=0.015)) =nan;
        [C,h] = contourf(char_freq.Pzeta_norm,char_freq.E,l_passing2,v,'facecolor','none','linewidth',1,'linecolor','C','linestyle',':');hold on;
        %clabel(C,h,'FontSize',16,'Color','k','labelspacing',700);
    end


    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.E,l_passing1,v,'facecolor','none','linewidth',2,'linecolor','k','linestyle',':');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',700);
    v = 5:1:8;
    l_trap_potato = l_res;
    l_trap_potato(find(~(strcmp(char_freq.type,'trapped')|strcmp(char_freq.type,'potato')))) = nan;

    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.E,l_trap_potato,v,'facecolor','none','linewidth',2,'linecolor','k','linestyle',':');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',800);

    %title('l\_res\_d');
    E_convert = 1; % when using ps_cls.Pzeta_E_fixed_mu()
    h1 = plot(MG.pzeta_norm,MG.E_pzeta(1,:)*E_convert,'b','linewidth',2);hold on;
    h2 = plot(MG.pzeta_norm,MG.E_pzeta(2,:)*E_convert,'g','linewidth',2);hold on;
    h3 = plot(MG.pzeta_norm,MG.E_pzeta(3,:)*E_convert,'r','linewidth',2);hold on;
    h4 = plot(MG.pzeta_norm,MG.E_pzeta(4,:)*E_convert,'color',[1,0.7,0],'linewidth',2);hold on;
    h5 = plot(MG.pzeta_norm,MG.E_pzeta(5,:)*E_convert,'m','linewidth',2);hold on;
    h7 = plot(MG.Pzeta_bound5,MG.E_uni_grid1*E_convert,'c--','linewidth',2);hold on;

    %h6 = plot(MG.pzeta_norm,MG.E_pzeta(6,:)*E_convert,'k--','linewidth',2);hold on;
    [C,h21] =contour(char_freq.Pzeta_norm, char_freq.E,char_freq.q_avrg, [0 -3], 'LineWidth', 2, 'Color','k','linestyle','--');hold on
    axis([MG.pzeta_norm0 MG.pzeta_norm1 14 E_max]);
    % axis([-1 MG.pzeta_norm1 14 E_max]);

    %---------------------------------------------------------------------------------------------图例
    % if strcmp(passing_option,'counter-passing')
    %     posx = min(xlim) + 0.03*range(xlim);
    %     posy = min(ylim) + 0.13*range(ylim);
    % elseif strcmp(passing_option,'co-passing')
    %     posx = min(xlim) + 0.05*range(xlim);
    %     posy = min(ylim) + 0.77*range(ylim);
    % end
    %
    % w = 0.45*range(xlim);
    % h = 0.19*range(ylim);
    % fill([posx posx+w posx+w posx],[posy-h/2 posy-h/2 posy+h posy+h],'w',...
    %     'EdgeColor','k','LineWidth',1.0);
    %
    % % 小椭圆渐变
    % x0 = posx + 0.08*range(xlim);
    % y0 = posy + 0.15*range(ylim);
    % a = 0.03 * range(xlim);
    % b = 0.03 * range(ylim);
    % N = 200; theta = linspace(0,2*pi,N);
    % cmap = colormap('parula'); nC = size(cmap,1);
    % for k = 1:nC
    %     r = 1 - (k-1)/(nC-1);
    %     fill(x0 + a*r*cos(theta), y0 + b*r*sin(theta), cmap(k,:), ...
    %         'EdgeColor','none');
    % end
    % plot(x0 + a*cos(theta), y0 + b*sin(theta), 'k', 'LineWidth', 1.2);
    %
    % % 图例文字
    % text(x0 + 0.06*range(xlim), y0, '\it{p}', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    % text(x0 - 0.045*range(xlim), y0-h/1.9, ' \cdot\cdot\cdot \it{p} \in ℤ', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    % text(x0 - 0.05*range(xlim), y0-h/1, '  --- \langle q \rangle = 3', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    %--------------------------------------------------------------------------------------------------
    %set(gca,'TickDir','out','box','off');
    set(gca,'LineWidth',1.3,'FontSize',20,'ticklength',[0.008 0.008])

    %grid on;
    
xlabel('$P_\zeta/\psi_w$','fontsize',23,'Interpreter','latex')
    ylabel('$E$(keV)','fontsize',23,'Interpreter','latex')
    if strcmp(passing_option,'counter-passing')
        %saveas(gca,'../output/co&counter_MG','epsc');
    elseif strcmp(passing_option,'co-passing')
        %saveas(gca,'../output/pzeta_E_fix_mu_co','epsc');
    end

elseif iplot == 2


    l_res = (omega0 - n*char_freq.omega_phi)./char_freq.omega_b;
    if strcmp(passing_option,'counter-passing')
        v = 12:1:23;
        %l_res = l_res1;
    elseif strcmp(passing_option,'co-passing')
        v = 4:1:14;
       % l_res = l_res;
        for pp=1:size(char_freq.Pzeta_norm,1)
            for ppp=1:size(char_freq.Pzeta_norm,2)
                fdp(pp,ppp)=min(min(((char_freq.Pzeta_norm(pp,ppp)-char_freq1.Pzeta_norm).^2+(char_freq.lambda(pp,ppp)-char_freq1.lambda).^2).^0.5));
            end
        end
    end

    %contourf(char_freq.Pzeta_norm,char_freq.lambda,l_res,100,'edgecolor','none');colormap('parula');hold on;caxis([v(1) v(end)])
    %cc=colorbar('northoutside');
    % title('l\_res\_d');
    %contourf(char_freq.Pzeta_norm,char_freq.lambda,l_res-round(l_res),100,'edgecolor','none');title('l\_res\_d');colorbar;colormap('jet');hold on;caxis([-0.5 0.5])
    l_passing = l_res;
    l_passing(find(~(strcmp(char_freq.type,'co-passing')|strcmp(char_freq.type,'counter-passing')))) = nan;
    l_passing1 = l_passing;
    l_passing2 = l_passing;
    if strcmp(passing_option,'counter-passing')
        v = 12:1:17;
    elseif strcmp(passing_option,'co-passing')
        v = 7:1:12;

        l_passing1(find(fdp<0.015))=nan;
        l_passing2(find(fdp>=0.015)) =nan;
        [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,l_passing2,v,'facecolor','none','linewidth',1,'linecolor','C','linestyle',':');hold on;
        %clabel(C,h,'FontSize',16,'Color','k','labelspacing',700);
    end


    %l_passing1(find(char_freq.lambda>GB.lambda_pzeta(1,:)) = nan;
    %l_passing1(find(~(strcmp(char_freq.type,'co-passing')|strcmp(char_freq.type,'counter-passing')))) = nan;


    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,l_passing1,v,'facecolor','none','linewidth',2,'linecolor','k');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',800);
    set(h,'LineWidth',2,'linecolor','k','linestyle',':')


    v = 1:1:8;
    l_trap_potato = l_res;
    l_trap_potato(find(~(strcmp(char_freq.type,'trapped')|strcmp(char_freq.type,'potato')))) = nan;

    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,l_trap_potato,v,'facecolor','none','linewidth',2,'linecolor','k');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',800);
    set(h,'LineWidth',2,'linecolor','k','linestyle',':')
    [C,h21] =contour(char_freq.Pzeta_norm, char_freq.lambda,char_freq.q_avrg, [0 -3], 'LineWidth', 2, 'Color','k','linestyle','--');hold on

    h1 = plot(GB.pzeta_norm,GB.lambda_pzeta(1,:),'b','linewidth',2);hold on;
    h2 = plot(GB.pzeta_norm,GB.lambda_pzeta(2,:),'g','linewidth',2);hold on;
    h3 = plot(GB.pzeta_norm,GB.lambda_pzeta(3,:),'r','linewidth',2);hold on;
    h4 = plot(GB.pzeta_norm,GB.lambda_pzeta(4,:),'color',[1,0.7,0],'linewidth',2);hold on;
    h5 = plot(GB.pzeta_norm,GB.lambda_pzeta(5,:),'m','linewidth',2);hold on;
    % plot(GB.Pzeta_bound0,GB.lambda_uni_grid,'k','linewidth',3);hold on;
    plot(GB.Pzeta_bound1,GB.lambda_uni_grid,'c--','linewidth',2);hold on;

    % legend([h1 h2 h3 h4 h5],'right-wall','left-wall','axis','barely-trap','deep-trap')
    %grid on;
    % strings = ['\muB_a = ',num2str(GB.E_perp_axis),' keV'];
    %annotation('textbox',[0.2 0.28 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none')
    axis([GB.pzeta_norm0 GB.pzeta_norm1 0 max(GB.lambda_pzeta(1,:))]);

    %---------------------------------------------------------------------------------------------图例
    % if strcmp(passing_option,'counter-passing')
    %     posx = min(xlim) + 0.03*range(xlim);
    %     posy = min(ylim) + 0.77*range(ylim);
    % elseif strcmp(passing_option,'co-passing')
    %     posx = min(xlim) + 0.03*range(xlim);
    %     posy = min(ylim) + 0.13*range(ylim);
    % end
    % 
    % w = 0.45*range(xlim);
    % h = 0.19*range(ylim);
    % fill([posx posx+w posx+w posx],[posy-h/2 posy-h/2 posy+h posy+h],'w',...
    %     'EdgeColor','k','LineWidth',1.0);
    % 
    % % 小椭圆渐变
    % x0 = posx + 0.08*range(xlim);
    % y0 = posy + 0.15*range(ylim);
    % a = 0.03 * range(xlim);
    % b = 0.03 * range(ylim);
    % N = 200; theta = linspace(0,2*pi,N);
    % cmap = colormap('parula'); nC = size(cmap,1);
    % for k = 1:nC
    %     r = 1 - (k-1)/(nC-1);
    %     fill(x0 + a*r*cos(theta), y0 + b*r*sin(theta), cmap(k,:), ...
    %         'EdgeColor','none');
    % end
    % plot(x0 + a*cos(theta), y0 + b*sin(theta), 'k', 'LineWidth', 1.2);
    % 
    % % 图例文字
    % text(x0 + 0.06*range(xlim), y0, '\it{p}', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    % text(x0 - 0.045*range(xlim), y0-h/1.9, ' \cdot\cdot\cdot \it{p} \in ℤ', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    % text(x0 - 0.05*range(xlim), y0-h/1, '  --- \langle q \rangle = 3', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    %--------------------------------------------------------------------------------------------------
        set(gca,'LineWidth',1.3,'FontSize',20,'ticklength',[0.008 0.008])
    xlabel('$P_\zeta/\psi_w$','fontsize',23,'Interpreter','latex')
    ylabel('$\lambda$','fontsize',23,'Interpreter','latex')

    if strcmp(passing_option,'counter-passing')
        %saveas(gca,'../output/pzeta_lambda_fix_mu_counter','epsc');
    elseif strcmp(passing_option,'co-passing')
        %saveas(gca,'../output/pzeta_lambda_fix_mu_co','epsc');
    end
elseif iplot == 3
    colors = {'parula','jet','hot','turbo'}; % 不同 colormap
    fields = {'omega_b','omega_d','omega_phi','q_avrg'};
  
    l_res = (omega0 - n*char_freq.omega_phi)./char_freq.omega_b; % Eq. (70) in Bao et al 24 NF

    %l_res(find(strcmp(char_freq.type,'trapped'))) = nan;
    %l_res(find(strcmp(char_freq.type,'potato'))) = nan;
    %l_res(find(strcmp(char_freq.type,'stagnation'))) = nan;
    %l_res(find(strcmp(char_freq.type,'co-passing'))) = nan;
    %l_res(find(strcmp(char_freq.type,'counter-passing'))) = nan;
    l_res(find(strcmp(char_freq.exist,'no'))) = nan;
    l_res(find(strcmp(char_freq.confine,'lost'))) = nan;

    if strcmp(passing_option,'counter-passing')
        v = 12:1:23;
    elseif strcmp(passing_option,'co-passing')
        v = -2:1:14;
         for pp=1:size(char_freq.Pzeta_norm,1)
            for ppp=1:size(char_freq.Pzeta_norm,2)
                fdp(pp,ppp)=min(min(((char_freq.Pzeta_norm(pp,ppp)-char_freq1.Pzeta_norm).^2+(char_freq.lambda(pp,ppp)-char_freq1.lambda).^2).^0.5));
            end
        end
    end

    l_passing = l_res;
    l_passing(find(~(strcmp(char_freq.type,'co-passing')|strcmp(char_freq.type,'counter-passing')))) = nan;
    l_passing1 = l_passing;
    l_passing2 = l_passing;
    if strcmp(passing_option,'counter-passing')
        v = 12:1:17;
    elseif strcmp(passing_option,'co-passing')
        v = 7:1:19;
        l_passing1(find(fdp<0.015))=nan;
        l_passing2(find(fdp>=0.015)) =nan;
        [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,l_passing2,v,'facecolor','none','linewidth',1,'linecolor','C','linestyle',':');hold on;
        %clabel(C,h,'FontSize',16,'Color','k','labelspacing',700);
    end
    %v = 12:1:24
    %contourf(char_freq.Pzeta_norm,char_freq.lambda,l_res,100,'edgecolor','none');colorbar('northoutside');colormap('jet');hold on;caxis([v(1) v(end)])
    %title('l\_res\_d');
    
    v = 7:1:17;
    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,l_passing1,v,'facecolor','none','linewidth',2,'linecolor','k');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',400);
    set(h,'LineWidth',2,'linecolor','k','linestyle',':')
    v = -0:1:4;
    l_trap_potato = l_res;
    l_trap_potato(find(~(strcmp(char_freq.type,'trapped')|strcmp(char_freq.type,'potato')))) = nan;
    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.lambda,l_trap_potato,v,'facecolor','none','linewidth',2,'linecolor','k');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',400);
    set(h,'LineWidth',2,'linecolor','k','linestyle',':')

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
        set(gca,'LineWidth',1.3,'FontSize',20,'ticklength',[0.008 0.008])

    xlabel('$P_\zeta/\psi_w$','fontsize',23,'Interpreter','latex')
    ylabel('$\lambda$','fontsize',23,'Interpreter','latex')
   % strings = ['E = ',num2str(PLam_2D.energy_out),'keV'];
    %annotation('textbox',[0.64 0.79 0.9 0.03],'String',strings,'fontsize',25,'edgecolor','none')

    %---------------------------------------------------------------------------------------------图例
    % if strcmp(passing_option,'counter-passing')
    %     posx = min(xlim) + 0.05*range(xlim);
    %     posy = min(ylim) + 0.77*range(ylim);
    % elseif strcmp(passing_option,'co-passing')
    %     posx = min(xlim) + 0.05*range(xlim);
    %     posy = min(ylim) + 0.12*range(ylim);
    % end
    % 
    % w = 0.45*range(xlim);
    % h = 0.19*range(ylim);
    % fill([posx posx+w posx+w posx],[posy-h/2 posy-h/2 posy+h posy+h],'w',...
    %     'EdgeColor','k','LineWidth',1.0);
    % 
    % % 小椭圆渐变
    % x0 = posx + 0.08*range(xlim);
    % y0 = posy + 0.15*range(ylim);
    % a = 0.03 * range(xlim);
    % b = 0.03 * range(ylim);
    % N = 200; theta = linspace(0,2*pi,N);
    % cmap = colormap(colors{1}); nC = size(cmap,1);
    % for k = 1:nC
    %     r = 1 - (k-1)/(nC-1);
    %     fill(x0 + a*r*cos(theta), y0 + b*r*sin(theta), cmap(k,:), ...
    %         'EdgeColor','none');
    % end
    % plot(x0 + a*cos(theta), y0 + b*sin(theta), 'k', 'LineWidth', 1.2);
    % 
    % % 图例文字
    % text(x0 + 0.06*range(xlim), y0, '\it{p}', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    % text(x0 - 0.045*range(xlim), y0-h/1.9, ' \cdot\cdot\cdot \it{p} \in ℤ', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    % text(x0 - 0.05*range(xlim), y0-h/1, '  --- \langle q \rangle = 3', ...
    %     'FontSize',23,'FontName','Times New Roman', ...
    %     'FontWeight','bold','Interpreter','tex', ...
    %     'HorizontalAlignment','left','VerticalAlignment','middle');
    %--------------------------------------------------------------------------------------------------
    % if strcmp(passing_option,'counter-passing')
    %     saveas(gca,'../output/pzeta_lambda_2d_phase_space_counter','epsc');
    % elseif strcmp(passing_option,'co-passing')
    %     saveas(gca,'../output/pzeta_lambda_2d_phase_space_co','epsc');
    % end


end
