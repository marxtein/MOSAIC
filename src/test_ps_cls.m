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
global particle eq PE_2D PLam_2D MG GB

if strcmp(ps_option,'PLam') | strcmp(ps_option,'PLam_traj')  | strcmp(ps_option,'PE') | strcmp(ps_option,'PLam_traj2')

    %if strcmp(ps_option,'PLam_traj')
        by=0;%0:Not drawing the sixth boundary to save time
    %end
moduleName = '3D phase space boundary construction ';
fprintf('Starting module: %s...\n', moduleName);
tic;
A = ps_cls.surf_Pzeta_E_lambda(particle.charge,particle.mass,E0,E1,E_grid,Pzeta_grid,by);
elapsedTime = toc;
fprintf('=> Module [%s] execution time: %.4f seconds\n\n', moduleName, elapsedTime);
%% Pzeta-E-lambda 3D phase space boundary
%figure;
 figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.33,0.5],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',1,...
        'DefaultAxesTickLength',[0.02,0.015]);
%contourf(squeeze(pzeta(1,:,:)),squeeze(lambda(1,:,:)),squeeze(energy_3D(1,:,:)),'edgecolor','none')

h1 = surf(squeeze(A.Pzeta(1,:,:)),squeeze(A.energy(1,:,:)),squeeze(A.lambda(1,:,:)),'EdgeColor','none','FaceColor','b','FaceAlpha',0.4);hold on;
h2 = surf(squeeze(A.Pzeta(2,:,:)),squeeze(A.energy(2,:,:)),squeeze(A.lambda(2,:,:)),'EdgeColor','none','FaceColor','g','FaceAlpha',0.6);hold on;
h3 = surf(squeeze(A.Pzeta(3,:,:)),squeeze(A.energy(3,:,:)),squeeze(A.lambda(3,:,:)),'EdgeColor','none','FaceColor','r','FaceAlpha',0.15);hold on;
h4 = surf(squeeze(A.Pzeta(4,:,:)),squeeze(A.energy(4,:,:)),squeeze(A.lambda(4,:,:)),'EdgeColor','none','FaceColor','y','FaceAlpha',0.4);hold on;
h5 = surf(squeeze(A.Pzeta(5,:,:)),squeeze(A.energy(5,:,:)),squeeze(A.lambda(5,:,:)),'EdgeColor','none','FaceColor','m','FaceAlpha',0.5);hold on;
%h6 = surf(squeeze(A.Pzeta(6,:,:)),squeeze(A.energy(6,:,:)),squeeze(A.lambda(6,:,:)),'EdgeColor','none','FaceColor','k','FaceAlpha',0.3);hold on;
h7 = surf(squeeze(A.Pzeta(7,:,:)),squeeze(A.energy(7,:,:)),squeeze(A.lambda(7,:,:)),'EdgeColor','none','FaceColor','c','FaceAlpha',0.4);hold on;


% legend([h1,h2,h3,h4,h5,h7],...
%     'Boundary 1',...
%     'Boundary 2',...
%     'Boundary 3',...
%     'Boundary 4',...
%     'Boundary 5',...
%     'Boundary 6',...
%     'FontSize',18,'Location','northeast');
axis([floor(min(A.Pzeta(1,end,:)))    (max(A.Pzeta(3,end,:)))+0.1    E0   E1   0    1.4000]);
set(gca,'LineWidth',1,'FontSize',19,'ticklength',[0.02 0.0])
xlabel('$P_\zeta/\psi_w$','FontSize',25,'Interpreter','latex');
ylabel('$E$(keV)','FontSize',25,'Interpreter','latex');
zlabel('$\lambda=\mu B_a/E$','FontSize',25,'Interpreter','latex');
%strings = '3D phase space'
%annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none')
saveas(gca,'../output/3d_phase_space','epsc')
elseif strcmp(ps_option,'MG')
    A =ps_cls2.surf_Pzeta_E_mu(particle.charge,particle.mass,E0,E1,E_grid,Pzeta0,Pzeta1,Pzeta_grid,E_max,E_grid_num,psi_diag_norm)
           
    %% Pzeta-E-lambda 3D phase space boundary
    %figure;

     figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.33,0.5],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.015]);
    %contourf(squeeze(pzeta(1,:,:)),squeeze(lambda(1,:,:)),squeeze(energy_3D(1,:,:)),'edgecolor','none')
    h1=surf(squeeze(A.Pzeta(1,:,:)),squeeze(A.energy(1,:,:)),squeeze(A.mu(1,:,:)),'EdgeColor','none','FaceColor','b','FaceAlpha',0.4);hold on;
    h2=surf(squeeze(A.Pzeta(2,:,:)),squeeze(A.energy(2,:,:)),squeeze(A.mu(2,:,:)),'EdgeColor','none','FaceColor','g','FaceAlpha',0.4);hold on;
    h3=surf(squeeze(A.Pzeta(3,:,:)),squeeze(A.energy(3,:,:)),squeeze(A.mu(3,:,:)),'EdgeColor','none','FaceColor','r','FaceAlpha',0.15);hold on;
    h4=surf(squeeze(A.Pzeta(4,:,:)),squeeze(A.energy(4,:,:)),squeeze(A.mu(4,:,:)),'EdgeColor','none','FaceColor','y','FaceAlpha',0.4);hold on;
    h5=surf(squeeze(A.Pzeta(5,:,:)),squeeze(A.energy(5,:,:)),squeeze(A.mu(5,:,:)),'EdgeColor','none','FaceColor','m','FaceAlpha',0.5);hold on;
    %h6=surf(squeeze(A.Pzeta(6,:,:)),squeeze(A.energy(6,:,:)),squeeze(A.mu(6,:,:)),'EdgeColor','none','FaceColor','k','FaceAlpha',0.5);hold on;
    h7=surf(squeeze(A.Pzeta(7,:,:)),squeeze(A.energy(7,:,:)),squeeze(A.mu(7,:,:)),'EdgeColor','none','FaceColor','c','FaceAlpha',0.3);hold on;
    legend([h1,h2,h3,h4,h5,h7],...
        'Boundary 1',...
        'Boundary 2',...
        'Boundary 3',...
        'Boundary 4',...
        'Boundary 5',...
        'Boundary 6',...
        'FontSize',18,'Location','northeast');
    axis([floor(min(A.Pzeta(1,end,:)))    ceil(max(A.Pzeta(3,end,:)))    E0   E1   0    92]);
    set(gca,'LineWidth',2,'FontSize',19,'ticklength',[0.02 0.0])
    xlabel('$P_\zeta/\psi_w$','fontsize',25,'Interpreter','latex')
    ylabel('$\mu B_a$(keV)','fontsize',25,'Interpreter','latex')
    zlabel('$E$(keV)','fontsize',25,'Interpreter','latex')
    %title('3D phase space','fontsize',30)
    %strings = '3D phase space'
    %annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',22,'edgecolor','none')
    saveas(gca,'../output/PEmu','epsc')

elseif strcmp(ps_option,'GB')
    %A = ps_cls.surf_Pzeta_E_lambda(particle.charge,particle.mass,E0,E1,E_grid,Pzeta_grid);
   A =ps_cls2.surf_Pzeta_lambda_mu(particle.charge,particle.mass,E0,E1,E_grid,Pzeta0,Pzeta1,Pzeta_grid,lam_min,lam_grid_num)
           
   %(qpart,apart,E0_keV,E1_keV,E_grid,Pzeta0,Pzeta1,Pzeta_grid,lam_min,lam_grid_num)
          
    %% Pzeta-E-lambda 3D phase space boundary
    %figure;
    figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.33,0.5],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.015]);
    %contourf(squeeze(pzeta(1,:,:)),squeeze(lambda(1,:,:)),squeeze(energy_3D(1,:,:)),'edgecolor','none')
    h1=surf(squeeze(A.Pzeta(1,:,:)),squeeze(A.energy(1,:,:)),squeeze(A.lambda(1,:,:)),'EdgeColor','none','FaceColor','b','FaceAlpha',0.5);hold on;
    h2=surf(squeeze(A.Pzeta(2,:,:)),squeeze(A.energy(2,:,:)),squeeze(A.lambda(2,:,:)),'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);hold on;
    h3=surf(squeeze(A.Pzeta(3,:,:)),squeeze(A.energy(3,:,:)),squeeze(A.lambda(3,:,:)),'EdgeColor','none','FaceColor','r','FaceAlpha',0.5);hold on;
    h4=surf(squeeze(A.Pzeta(4,:,:)),squeeze(A.energy(4,:,:)),squeeze(A.lambda(4,:,:)),'EdgeColor','none','FaceColor','y','FaceAlpha',0.5);hold on;
    h5=surf(squeeze(A.Pzeta(5,:,:)),squeeze(A.energy(5,:,:)),squeeze(A.lambda(5,:,:)),'EdgeColor','none','FaceColor','m','FaceAlpha',0.5);hold on;
    %h6=surf(squeeze(A.Pzeta(6,:,:)),squeeze(A.energy(6,:,:)),squeeze(A.lambda(6,:,:)),'EdgeColor','none','FaceColor','k','FaceAlpha',0.5);hold on;
    h7=surf(squeeze(A.Pzeta(7,:,:)),squeeze(A.energy(7,:,:)),squeeze(A.lambda(7,:,:)),'EdgeColor','none','FaceColor','c','FaceAlpha',0.3);hold on;
    
    legend([h1,h2,h3,h4,h5,h7],...
        'Boundary 1',...
        'Boundary 2',...
        'Boundary 3',...
        'Boundary 4',...
        'Boundary 5',...
        'Boundary 6',...
        'FontSize',18,'Location','northeast');
    axis([floor(min(A.Pzeta(1,end,:)))    ceil(max(A.Pzeta(3,end,:)))    E0   E1   0    1.4000]);
    %axis([floor(min(A.Pzeta(1,end,:)))    ceil(max(A.Pzeta(3,end,:)))    E0   E1   0    92]);
    set(gca,'LineWidth',2,'FontSize',19,'ticklength',[0.02 0.0])
    xlabel('$P_\zeta/\psi_w$','fontsize',25,'Interpreter','latex')
    ylabel('$\mu B_a$(keV)','fontsize',25,'Interpreter','latex')
    zlabel('$\lambda$','fontsize',25,'Interpreter','latex')
    %title('3D phase space','fontsize',30)
    %strings = '3D phase space'
    %annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none')
    saveas(gca,'../output/Plammu','epsc')


end
if strcmp(ps_option,'PE')
    
    %% Pzeta-E 2D phase space boundary
    ps_3D = A;
    PE_2D = ps_cls.cro_sec_Pzeta_E(lam_in,E_grid_num,ps_3D);
    
    figure;
    plot(PE_2D.Pzeta(1,:),PE_2D.energy(1,:),'b.');hold on;
    plot(PE_2D.Pzeta(2,:),PE_2D.energy(2,:),'g.');hold on;
    plot(PE_2D.Pzeta(3,:),PE_2D.energy(3,:),'r.');hold on;
    plot(PE_2D.Pzeta(4,:),PE_2D.energy(4,:),'y.');hold on;
    plot(PE_2D.Pzeta(5,:),PE_2D.energy(5,:),'m.');hold on;
    plot(PE_2D.Pzeta(6,:),PE_2D.energy(6,:),'k.');hold on;
    plot(PE_2D.Pzeta_bound_right(1,:),PE_2D.E_uni_grid,'k--','linewidth',3);hold on;
    plot(PE_2D.Pzeta_bound_right(3,:),PE_2D.E_uni_grid,'k--','linewidth',3);hold on;
    plot(PE_2D.Pzeta_bound_left(2,:),PE_2D.E_uni_grid,'c--','linewidth',3);hold on;
    plot(PE_2D.Pzeta_bound_left(3,:),PE_2D.E_uni_grid,'c--','linewidth',3);hold on;
    
    %plot(PE_2D.Pzeta_bound_left(4,:),PE_2D.E_uni_grid,'y--','linewidth',3);hold on;
    %plot(PE_2D.Pzeta_bound_left(5,:),PE_2D.E_uni_grid,'m--','linewidth',3);hold on;
    %plot(PE_2D.Pzeta_bound_left(6,:),PE_2D.E_uni_grid,'k--','linewidth',3);hold on;
    %plot(PE_2D.Pzeta_bound0,PE_2D.E_uni_grid,'linewidth',3);hold on;
    %plot(PE_2D.Pzeta_bound1,PE_2D.E_uni_grid,'linewidth',3);hold on;
    
    grid on;
    set(gca,'LineWidth',2,'FontSize',16,'ticklength',[0.02 0.0])
    xlabel('P_\zeta/\psi_w','fontsize',30);
    ylabel('Energy (keV)','fontsize',30);
    strings = ['\lambda = ',num2str(PE_2D.lambda0)];
    annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none')
    %saveas(gca,'./pzeta_E_2d_phase_space','epsc')
    
elseif strcmp(ps_option,'PLam') | strcmp(ps_option,'PLam_traj')  | strcmp(ps_option,'PLam_traj2')
    
    %% Pzeta-lambda 2D phase space boundary
    tic;
    if eq.gpsi(1,1)>0
        sign_gpsi = 'positive';
    else
        sign_gpsi = 'negative';
    end
    
    PLam_2D = ps_cls.cro_sec_Pzeta_lambda(energy_in,A,sign_gpsi);
    toc;
    figure;
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'color',[1,0.7,0],'linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'m','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    plot(PLam_2D.Pzeta_bound_right,PLam_2D.lam_uni_grid,'c-.','linewidth',2);hold on;

   
    % plot(PLam_2D.Pzeta_bound_left1,PLam_2D.lam_uni_grid,'r-.','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta_bound_left2,PLam_2D.lam_uni_grid,'r-.','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta_bound_right_co_EE,PLam_2D.lam_uni_grid,'r-.','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta_bound_counter_passing_left,PLam_2D.lam_uni_grid_counter_passing,'c--','linewidth',2);hold on;
    % plot(PLam_2D.Pzeta_bound_counter_passing_right,PLam_2D.lam_uni_grid_counter_passing,'c--','linewidth',2);hold on;
    grid on;
    axis([floor(min(PLam_2D.Pzeta_bound_left2))    ceil(max(PLam_2D.Pzeta_bound_right))    0    1.4000]);
    set(gca,'LineWidth',2,'FontSize',16,'ticklength',[0.02 0.0])
    xlabel('P_\zeta/\psi_w','fontsize',30)
    ylabel('\lambda=\muB_a/E','fontsize',30)
    strings = ['E = ',num2str(energy_in),'keV'];
    annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none')
    %saveas(gca,'./pzeta_lambda_2d_phase_space','epsc');
    
elseif strcmp(ps_option,'MG')| strcmp(ps_option,'Poincare')
    
    %% MG: Pzeta-E at fixed mu 2D phase space boundary
    MG = ps_cls2.Pzeta_E_fixed_mu(particle.charge,particle.mass,Eperp_in,Pzeta0,Pzeta1,Pzeta_grid_num,E_max,E_grid_num,psi_diag_norm);
      figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.25,0.4],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.015]);
    % h1 = plot(MG.pzeta_norm,MG.E_pzeta(1,:),'b','linewidth',2);hold on;
    % h2 = plot(MG.pzeta_norm,MG.E_pzeta(2,:),'g','linewidth',2);hold on;
    % h3 = plot(MG.pzeta_norm,MG.E_pzeta(3,:),'r','linewidth',2);hold on;
    % h4 = plot(MG.pzeta_norm,MG.E_pzeta(4,:),'c','linewidth',2);hold on;
    % h5 = plot(MG.pzeta_norm,MG.E_pzeta(5,:),'c--','linewidth',2);hold on;

    h1 = plot(MG.pzeta_norm,MG.E_pzeta(1,:),'b','linewidth',2);hold on;
    h2 = plot(MG.pzeta_norm,MG.E_pzeta(2,:),'g','linewidth',2);hold on;
    h3 = plot(MG.pzeta_norm,MG.E_pzeta(3,:),'r','linewidth',2);hold on;
    h4 = plot(MG.pzeta_norm,MG.E_pzeta(4,:),'color',[1,0.7,0],'linewidth',2);hold on;
    h5 = plot(MG.pzeta_norm,MG.E_pzeta(5,:),'m','linewidth',2);hold on;
    %h6 = plot(MG.pzeta_norm,MG.E_pzeta(6,:),'k--','linewidth',2);hold on;
    %plot(MG.Pzeta_bound0,MG.E_uni_grid,'k','linewidth',3);hold on;
    %plot(MG.Pzeta_bound1,MG.E_uni_grid,'k','linewidth',3);hold on;
    %plot(MG.Pzeta_bound2,MG.E_uni_grid,'k','linewidth',3);hold on;
    %plot(MG.Pzeta_bound3,MG.E_uni_grid,'k','linewidth',3);hold on;
     plot(MG.Pzeta_bound5,MG.E_uni_grid1,'c--','linewidth',3);hold on;
    %legend([h1 h2 h3 h4 h5 h6],'right-wall','left-wall','axis','barely-trap','deep-trap','psi-diag-zero-theta')
    grid on;
    strings = ['\muB_a = ',num2str(MG.E_perp_axis),' keV'];
    annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none')
    axis([MG.pzeta_norm0 MG.pzeta_norm1 10 E_max]);
    xlabel('P_\zeta/\psi_w','fontsize',30);
    ylabel('E (keV)','fontsize',30);
    set(gca,'fontsize',16);
    
elseif strcmp(ps_option,'GB')
    
    %% GB: Pzeta-lambda at fixed mu 2D phase space boundary
    if eq.gpsi(1,1)>0
        sign_gpsi = 'positive';
    else
        sign_gpsi = 'negative';
    end
    GB = ps_cls2.Pzeta_lambda_fixed_mu(particle.charge,particle.mass,Eperp_in,Pzeta0,Pzeta1,Pzeta_grid_num,lam_min,lam_grid_num,sign_gpsi);
    
    figure;
    h1 = plot(GB.pzeta_norm,GB.lambda_pzeta(1,:),'b','linewidth',2);hold on;
    h2 = plot(GB.pzeta_norm,GB.lambda_pzeta(2,:),'g','linewidth',2);hold on;
    h3 = plot(GB.pzeta_norm,GB.lambda_pzeta(3,:),'r','linewidth',2);hold on;
    h4 = plot(GB.pzeta_norm,GB.lambda_pzeta(4,:),'color',[1,0.7,0],'linewidth',2);hold on;
    h5 = plot(GB.pzeta_norm,GB.lambda_pzeta(5,:),'m','linewidth',2);hold on;
    %plot(GB.Pzeta_bound0,GB.lambda_uni_grid,'k','linewidth',3);hold on;
    %plot(GB.Pzeta_bound1,GB.lambda_uni_grid,'k','linewidth',3);hold on;
   % plot(GB.Pzeta_bound2,GB.lambda_uni_grid,'k','linewidth',3);hold on;
   % plot(GB.Pzeta_bound3,GB.lambda_uni_grid,'k','linewidth',3);hold on;
    %plot(GB.Pzeta_bound0_counter_passing,GB.lambda_uni_grid_counter_passing,'m','linewidth',3);hold on;
    %plot(GB.Pzeta_bound1_counter_passing,GB.lambda_uni_grid_counter_passing,'m','linewidth',3);hold on;
    plot(GB.Pzeta_bound1,GB.lambda_uni_grid,'c--','linewidth',3);hold on;
    %legend([h1 h2 h3 h4 h5],'right-wall','left-wall','axis','barely-trap','deep-trap')
    grid on;
    strings = ['\muB_a = ',num2str(GB.E_perp_axis),' keV'];
    annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none')
    axis([GB.pzeta_norm0 GB.pzeta_norm1 0 max(GB.lambda_pzeta(1,:))]);
    xlabel('P_\zeta/\psi_w','fontsize',30);
    ylabel('\lambda','fontsize',30);
    set(gca,'fontsize',16);
    
end






