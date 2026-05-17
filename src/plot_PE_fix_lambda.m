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

n = 4;
omega0 = 82.03*1000*2*pi;



if iplot == 1
    
    figure;
    subplot(141)
    contourf(char_freq.Pzeta_norm,char_freq.E,char_freq.omega_b);title('\omega\_b');colorbar;
    xlabel('P_\zeta/\psi_w');
    ylabel('E (keV)');
    subplot(142)
    contourf(char_freq.Pzeta_norm,char_freq.E,char_freq.omega_d);title('\omega\_d');colorbar;
    xlabel('P_\zeta/\psi_w');
    ylabel('E (keV)');
    subplot(143)
    contourf(char_freq.Pzeta_norm,char_freq.E,char_freq.omega_phi);title('\omega\_\phi');colorbar;
    xlabel('P_\zeta/\psi_w');
    ylabel('E (keV)');
    subplot(144)
    contourf(char_freq.Pzeta_norm,char_freq.E,char_freq.q_avrg);title('q\_avrg');colorbar;
    xlabel('P_\zeta/\psi_w');
    ylabel('E (keV)');
    
    
elseif iplot == 2
    
    figure;
    
    l_res = (omega0 - n*char_freq.omega_phi)./char_freq.omega_b;
    
    v = 5:1:14 % co-passing
    %v = 10:1:30 % counter-passing
    contourf(char_freq.Pzeta_norm,char_freq.E,l_res,100,'edgecolor','none');title('l\_res\_d');colorbar;colormap('jet');hold on;caxis([v(1) v(end)])
    [C,h] = contourf(char_freq.Pzeta_norm,char_freq.E,l_res,v,'facecolor','none','linewidth',2,'linecolor','k');hold on;
    clabel(C,h,'FontSize',16,'Color','k','FontWeight','bold','labelspacing',800);
    set(h,'LineWidth',2,'linecolor','k','linestyle','-')
    
    %{
    contourf(char_freq.Pzeta_norm,char_freq.E,l_res);title('l\_res\_d');colorbar;hold on;
    contour(char_freq.Pzeta_norm,char_freq.E,l_res,[-5 -5],'k','linewidth',3);hold on;
    contour(char_freq.Pzeta_norm,char_freq.E,l_res,[-6 -6],'g','linewidth',3);hold on;
    contour(char_freq.Pzeta_norm,char_freq.E,l_res,[-7 -7],'r','linewidth',3);hold on;
    contour(char_freq.Pzeta_norm,char_freq.E,l_res,[-8 -8],'y','linewidth',3);hold on;
    contour(char_freq.Pzeta_norm,char_freq.E,l_res,[-9 -9],'c','linewidth',3);hold on;
    contour(char_freq.Pzeta_norm,char_freq.E,l_res,[-14 -14],'m','linewidth',3);hold on;
    %}
    %xlabel('P_\zeta/\psi_w');
    %ylabel('E (keV)');
    
    plot(PE_2D.Pzeta(1,:),PE_2D.energy(1,:),'b.');hold on;
    plot(PE_2D.Pzeta(2,:),PE_2D.energy(2,:),'g.');hold on;
    plot(PE_2D.Pzeta(3,:),PE_2D.energy(3,:),'r.');hold on;
    
    plot(PE_2D.Pzeta_bound_right(1,:),PE_2D.E_uni_grid,'k--','linewidth',3);hold on;
    plot(PE_2D.Pzeta_bound_right(3,:),PE_2D.E_uni_grid,'k--','linewidth',3);hold on;
    plot(PE_2D.Pzeta_bound_left(2,:),PE_2D.E_uni_grid,'c--','linewidth',3);hold on;
    plot(PE_2D.Pzeta_bound_left(3,:),PE_2D.E_uni_grid,'c--','linewidth',3);hold on;
    
    
    grid on;
    set(gca,'LineWidth',2,'FontSize',16,'ticklength',[0.02 0.0])
    xlabel('P_\zeta/\psi_w','fontsize',30);
    ylabel('Energy (keV)','fontsize',30);
    strings = ['\lambda = ',num2str(PE_2D.lambda0)];
    annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none')
    %saveas(gca,'./pzeta_E_2d_phase_space','epsc')
    

end

    
