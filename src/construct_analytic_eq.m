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
global eq grids 
eqplot=0;

%eq.qcoef=[1.797;0.8;-0.2]; % q=q(1) + q(2)*psi_norm +q(3)*psi_norm^2

eq.lsp = 1000;
eq.lst = 1280;
%eq.psiw = 0.026;
eq.psia = 0.0001*eq.psiw;
eq.dpsi = (eq.psiw-eq.psia)/(eq.lsp-1);

%geometry 1d profile
eq.psi=eq.psia+eq.dpsi*(linspace(1,eq.lsp,eq.lsp)-1);
eq.psi_norm=eq.psi/eq.psiw;

% q=q(1) + q(2)*psi_norm +q(3)*psi_norm^2
eq.qpsi(1,:) = eq.qcoef(1) + eq.qcoef(2)*eq.psi_norm + eq.qcoef(3)*eq.psi_norm.^2;

eq.torpsi(1,:) = eq.psiw*eq.qcoef(1)*eq.psi_norm ...
    + (1/2)*eq.psiw*eq.qcoef(2)*eq.psi_norm.^2 ...
    + (1/3)*eq.psiw*eq.qcoef(3)*eq.psi_norm.^3;

eq.eps(1,:) = sqrt(2*eq.torpsi(1,:));

eq.torpsiw = max(eq.torpsi(1,:));

eq.r(1,:) = sqrt(eq.torpsi(1,:)/eq.torpsiw);

eq.rminor=max(eq.eps(1,:));

eq.gpsi(1,:)=ones(1,eq.lsp);

eq.ipsi(1,:)=eq.eps.^2./eq.qpsi;

eq.dtheta = 2*pi/(eq.lst-1);
eq.theta = eq.dtheta*(linspace(1,eq.lst,eq.lst)-1);

%x=rcos(theta)
eq.xsp(1,:,:)= 1 + eq.eps(1,:)'.*cos(eq.theta);% - eq.eps(1,:)'.^2.*(sin(eq.theta)).^2;
%z=rsin(theta)
eq.zsp(1,:,:)= eq.eps(1,:)'.*sin(eq.theta);% + eq.eps(1,:)'.^2.*(sin(eq.theta)).*(cos(eq.theta));
%b = 1 - epsilon*cos(theta) + epsilon^2*cos(theta)^2
%eq.bsp(1,:,:)= 1./(1 +( eq.eps(1,:)'.*cos(eq.theta) - eq.eps(1,:)'.^2.*sin(eq.theta).^2));%+ eq.rpsi(1,:).^2*cos(eq.theta(k))^2;
eq.bsp(1,:,:)= 1./(1 + eq.eps(1,:)'.*cos(eq.theta));

%%
% construct spline for 2d quantities
eq.bsp = qdspline.constructSpline2d(0, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
eq.zsp = qdspline.constructSpline2d(0, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.zsp);
eq.xsp = qdspline.constructSpline2d(0, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.xsp);

% construct spline for 1d quantities
eq.qpsi = qdspline.constructSpline1d(0, eq.lsp, eq.dpsi, eq.qpsi);
eq.torpsi = qdspline.constructSpline1d(0, eq.lsp, eq.dpsi, eq.torpsi);
eq.eps = qdspline.constructSpline1d(0, eq.lsp, eq.dpsi, eq.eps);
eq.r = qdspline.constructSpline1d(0, eq.lsp, eq.dpsi, eq.r);
eq.gpsi = qdspline.constructSpline1d(0, eq.lsp, eq.dpsi, eq.gpsi);
eq.ipsi = qdspline.constructSpline1d(0, eq.lsp, eq.dpsi, eq.ipsi);

%% fine grid
grids.psi0=0.01;
grids.psi1=0.98;

grids.mpsi=20;
grids.mtheta=200;

grids.psia=eq.psiw*grids.psi0; % simulation domain inner boundary
grids.psiw=eq.psiw*grids.psi1; % simulation domain outer boundary

grids.spdpsi = (grids.psiw-grids.psia)/(grids.mpsi-1);
grids.psi=grids.psia+grids.spdpsi*(linspace(1,grids.mpsi,grids.mpsi)-1);

grids.spdtheta = 2*pi/(grids.mtheta-1);
grids.theta = grids.spdtheta*(linspace(1,grids.mtheta,grids.mtheta)-1);

grids.X = qdspline.matrix2d_spline2d(0, grids.psi, grids.theta, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.xsp);
grids.Z = qdspline.matrix2d_spline2d(0, grids.psi, grids.theta, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.zsp);


%%  plasma profile
for i=1:grids.mpsi
    grids.qpsi(i) = qdspline.spline1d(0, grids.psi(i), eq.lsp, eq.dpsi, eq.qpsi);
    
    grids.torpsi(i) = qdspline.spline1d(0, grids.psi(i), eq.lsp, eq.dpsi, eq.torpsi);
    
    grids.eps(i) = qdspline.spline1d(0, grids.psi(i), eq.lsp, eq.dpsi, eq.eps);
    
    grids.r(i) = qdspline.spline1d(0, grids.psi(i), eq.lsp, eq.dpsi, eq.r);
    
    grids.gpsi(i) = qdspline.spline1d(0, grids.psi(i), eq.lsp, eq.dpsi, eq.gpsi);
    
    grids.ipsi(i) = qdspline.spline1d(0, grids.psi(i), eq.lsp, eq.dpsi, eq.ipsi);
end

%%
if eqplot == 1
    %% plot magnetic field
    figure('unit','normalized',...
        'DefaultAxesFontSize',15, ...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',1,...
        'DefaultAxesTickLength',[0.02,0.05]); % 2D view length versus 3D view length
    [C,h]=contourf(squeeze(eq.xsp(1,:,:)),squeeze(eq.zsp(1,:,:)),squeeze(eq.bsp(1,:,:)),20,'edgecolor','none');
    hold on;
    %plot(squeeze(eq.xsp(1,:,:)),squeeze(eq.zsp(1,:,:)),'k.-');
    hold on;
    %plot(squeeze(eq.xsp(1,:,:))',squeeze(eq.zsp(1,:,:))','k.-');
    colorbar;
    daspect([1 1 1]);
    title('Equilibrium magnetic field B_0/B_a','fontsize',20);
    xlabel('R/R_0','fontsize',20);
    ylabel('Z/R_0','fontsize',20);
end

