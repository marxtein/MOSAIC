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
test_opt = 3

if test_opt == 1
    
    Rtmp=spdata.xsp(1,end,125)/(spdata.rmaj/100);
    Ztmp=spdata.zsp(1,end,125)/(spdata.rmaj/100);
    
    psitmp = qdspline.spline2d(0, Ztmp + 0.5*eqRZ.zdim - eqRZ.zmid, Rtmp - eqRZ.rleft, eqRZ.nz, eqRZ.nr, eqRZ.delZ, eqRZ.delR, eqRZ.psisp)
    psitmp*eqRZ.psi_norm
    spdata.psiw
    B0_amp = qdspline.spline2d(0, Ztmp + 0.5*eqRZ.zdim - eqRZ.zmid, Rtmp - eqRZ.rleft, eqRZ.nz, eqRZ.nr, eqRZ.delZ, eqRZ.delR, eqRZ.b0_amp_sp)
    
elseif test_opt == 2
    
    Rtmp = squeeze(spdata.xsp(1,:,:))/(spdata.rmaj/100);
    Ztmp = squeeze(spdata.zsp(1,:,:))/(spdata.rmaj/100);
    psitmp = zeros(spdata.lsp,spdata.lst);
    B0_tmp = zeros(spdata.lsp,spdata.lst);
    
    for i = 1:spdata.lsp
        for j = 1:spdata.lst
            psitmp(i,j) = qdspline.spline2d(0, Ztmp(i,j) + 0.5*eqRZ.zdim - eqRZ.zmid, Rtmp(i,j) - eqRZ.rleft, eqRZ.nz, eqRZ.nr, eqRZ.delZ, eqRZ.delR, eqRZ.psisp);
            fpoltmp(i,j) = qdspline.spline2d(0, Ztmp(i,j) + 0.5*eqRZ.zdim - eqRZ.zmid, Rtmp(i,j) - eqRZ.rleft, eqRZ.nz, eqRZ.nr, eqRZ.delZ, eqRZ.delR, eqRZ.fpolsp);
            B0_tmp(i,j) = qdspline.spline2d(0, Ztmp(i,j) + 0.5*eqRZ.zdim - eqRZ.zmid, Rtmp(i,j) - eqRZ.rleft, eqRZ.nz, eqRZ.nr, eqRZ.delZ, eqRZ.delR, eqRZ.b0_amp_sp);
        end
    end
    
    if 0
        close all;
        figure('name','phase_space',...
            'unit','normalized',...
            'position',[0.0,0.0,0.6,0.6],... % figure position
            'DefaultAxesFontSize',15,...
            'DefaultAxesFontWeight','normal',...
            'DefaultAxesLineWidth',2,...
            'DefaultAxesTickLength',[0.02,0.05]);
        
        subplot(231)
        contourf(Rtmp,Ztmp,psitmp,'edgecolor','none');colorbar;
        xlabel('R/R0');
        ylabel('Z/R0');
        title('normalized \psi/(B_aR_0^2), ireg\_RZ')
        daspect([1 1 1]);
        
        subplot(232)
        contourf(Rtmp,Ztmp,fpoltmp,'edgecolor','none');colorbar;
        xlabel('R/R0');
        ylabel('Z/R0');
        title('normalized fpol/(B_aR_0), ireg\_RZ')
        daspect([1 1 1]);
        
        subplot(233)
        contourf(Rtmp,Ztmp,B0_tmp,'edgecolor','none');colorbar;
        xlabel('R/R0');
        ylabel('Z/R0');
        title('normalized B_0/B_a, ireg\_RZ')
        daspect([1 1 1]);
        
        subplot(234)
        a=eqRZ.psi2drz;
        a(find(eqRZ.psi2drz>eqRZ.psiw))=nan;
        contourf(eqRZ.R2drz,eqRZ.Z2drz,a,'edgecolor','none');colorbar;hold on;
        plot(efit.rbbbs/eqRZ.xnorm, efit.zbbbs/eqRZ.xnorm, 'k-','linewidth',2);hold on;
        plot(efit.rlim/eqRZ.xnorm, efit.zlim/eqRZ.xnorm, 'r-','linewidth',2);hold on;
        plot(efit.rmaxis/eqRZ.xnorm, efit.zmaxis/eqRZ.xnorm, 'go','markersize',5,'MarkerFaceColor','g');
        xlabel('R/R0');
        ylabel('Z/R0');
        title('normalized \psi/(B_aR_0^2), reg\_RZ')
        daspect([1 1 1]);
        
        subplot(235)
        a=eqRZ.fpol2drz;
        a(find(eqRZ.psi2drz>eqRZ.psiw))=nan;
        contourf(eqRZ.R2drz,eqRZ.Z2drz,a,'edgecolor','none');colorbar;hold on;
        plot(efit.rbbbs/eqRZ.xnorm, efit.zbbbs/eqRZ.xnorm, 'k-','linewidth',2);hold on;
        plot(efit.rlim/eqRZ.xnorm, efit.zlim/eqRZ.xnorm, 'r-','linewidth',2);hold on;
        plot(efit.rmaxis/eqRZ.xnorm, efit.zmaxis/eqRZ.xnorm, 'go','markersize',5,'MarkerFaceColor','g');
        xlabel('R/R0');
        ylabel('Z/R0');
        title('normalized fpol/(B_aR_0), reg\_RZ')
        daspect([1 1 1]);
        
        subplot(236)
        a=eqRZ.B0_amp;
        a(find(eqRZ.psi2drz>eqRZ.psiw))=nan;
        contourf(eqRZ.R2drz,eqRZ.Z2drz,a,'edgecolor','none');colorbar;hold on;
        plot(efit.rbbbs/eqRZ.xnorm, efit.zbbbs/eqRZ.xnorm, 'k-','linewidth',2);hold on;
        plot(efit.rlim/eqRZ.xnorm, efit.zlim/eqRZ.xnorm, 'r-','linewidth',2);hold on;
        plot(efit.rmaxis/eqRZ.xnorm, efit.zmaxis/eqRZ.xnorm, 'go','markersize',5,'MarkerFaceColor','g');
        xlabel('R/R0');
        ylabel('Z/R0');
        title('normalized B_0/B_a, reg\_RZ')
        daspect([1 1 1]);
    end
    
elseif test_opt == 3
    
    
    Booz = struct;
    
    Booz.psi = char_freq.psi_init;
    Booz.theta = char_freq.theta_init;
    Booz.E = char_freq.E;
    Booz.lambda = char_freq.lambda;
    Booz.Pzeta = char_freq.Pzeta_norm*eq.psiw;
    Booz.R = zeros(size(Booz.psi));
    Booz.Z = zeros(size(Booz.psi));
    Booz.B0 = zeros(size(Booz.psi));
    
    for i = 1:np
        if strcmp(char_freq.exist{i},'yes')
            Booz.R(i) = qdspline.spline2d(0, Booz.psi(i), Booz.theta(i), eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.xsp);
            Booz.Z(i) = qdspline.spline2d(0, Booz.psi(i), Booz.theta(i), eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.zsp);
            Booz.B0(i)= qdspline.spline2d(0, Booz.psi(i), Booz.theta(i), eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
        else
            Booz.R(i) = nan;
            Booz.Z(i) = nan;
            Booz.B0(i) = nan;
        end
    end
    
    RZ = struct;
    
    RZ.E = Booz.E;
    RZ.C_vpar = sign(char_freq.rhopara_init).*sign(particle.charge).*sqrt(1-Booz.lambda.*Booz.B0);
    RZ.R = Booz.R;
    RZ.Z = Booz.Z;
    RZ.psi = zeros(size(RZ.R));
    RZ.fpol = zeros(size(RZ.R));
    RZ.B0 = zeros(size(RZ.R));
    
    for i = 1:np
        RZ.psi(i) = qdspline.spline2d(0, RZ.Z(i) + 0.5*eqRZ.zdim - eqRZ.zmid, RZ.R(i) - eqRZ.rleft, eqRZ.nz, eqRZ.nr, eqRZ.delZ, eqRZ.delR, eqRZ.psisp);
        RZ.fpol(i) = qdspline.spline2d(0, RZ.Z(i) + 0.5*eqRZ.zdim - eqRZ.zmid, RZ.R(i) - eqRZ.rleft, eqRZ.nz, eqRZ.nr, eqRZ.delZ, eqRZ.delR, eqRZ.fpolsp);
        RZ.B0(i) = qdspline.spline2d(0, RZ.Z(i) + 0.5*eqRZ.zdim - eqRZ.zmid, RZ.R(i) - eqRZ.rleft, eqRZ.nz, eqRZ.nr, eqRZ.delZ, eqRZ.delR, eqRZ.b0_amp_sp);
    end
    
    RZ.E_gtc_unit = RZ.E*1000/unit.energy_norm;
    RZ.vpara = sqrt(2*RZ.E_gtc_unit/particle.mass).*RZ.C_vpar;
    RZ.rhopara = particle.mass./(particle.charge*RZ.B0).*RZ.vpara;
    RZ.Pzeta = RZ.fpol.*RZ.rhopara - RZ.psi;
    RZ.lambda = (1 - RZ.C_vpar.^2)./RZ.B0;
    
    
    
    figure;
    plot(PLam_2D.Pzeta(1,:),PLam_2D.lambda(1,:),'b','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(2,:),PLam_2D.lambda(2,:),'g','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(3,:),PLam_2D.lambda(3,:),'r','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),'k','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(5,:),PLam_2D.lambda(5,:),'k--','linewidth',2);hold on;
    plot(PLam_2D.Pzeta(6,:),PLam_2D.lambda(6,:),'k:','linewidth',2);hold on;
    grid on;
    
    for i = 1:np
          plot(Booz.Pzeta/eq.psiw, Booz.lambda,'mo','markersize',10,'markerfacecolor','m');hold on;
          plot(RZ.Pzeta/eq.psiw, RZ.lambda,'x','markersize',10,'markeredgecolor','k','linewidth',2);hold on;
    end
    h1=plot(nan,nan,'mo','markersize',10,'markerfacecolor','m');hold on;
    h2=plot(nan,nan,'x','markersize',10,'markeredgecolor','k','linewidth',2);hold on;
    
    
    axis([-2.0000    1.0000    0    1.4000]);
    set(gca,'LineWidth',2,'FontSize',16,'ticklength',[0.02 0.0])
    xlabel('P_\zeta/\psi_w','fontsize',30)
    ylabel('\lambda=\muB_a/E','fontsize',30)
    legend([h1 h2],'raw','from RZ')
    strings = ['E = ',num2str(PLam_2D.energy_out),'keV'];
    annotation('textbox',[0.2 0.88 0.9 0.03],'String',strings,'fontsize',30,'edgecolor','none');
end






