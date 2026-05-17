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
function  f = eqdataRZ_interface

global efit cocos

% Note: 'rz' as the suffix in variable name indicates horizontal and vertical
% labels (x,y) from the view of plotting. However, the array dimension obeys
% column first in matlab, so the length of vertical label is in front of the
% horizontal label for declaring the arrays, for example: *2drz = zeros(nz,nr).

f.bnorm = abs(efit.fpol(1)/efit.rmaxis);  % normalization is positive only.
f.xnorm = efit.rmaxis;  % normalization is positive only.
f.psi_norm = f.bnorm*f.xnorm^2;  % normalization is positive only.
f.psiw = efit.psiw/f.psi_norm;

f.nr = efit.nr;
f.nz = efit.nz;
f.R = efit.R/f.xnorm;
f.Z = efit.Z/f.xnorm;
f.R2drz = efit.R2drz/f.xnorm;
f.Z2drz = efit.Z2drz/f.xnorm;
f.rmaxis = efit.rmaxis/f.xnorm;  % R coordinate of magnetic axis
f.zmaxis = efit.zmaxis/f.xnorm;  % Z coordinate of magnetic axis
f.raxis_ind = efit.raxis_ind;  % R grid index on magnetic axis
f.zaxis_ind = efit.zaxis_ind;  % Z grid index on magnetic axis
f.fpol2d = griddedInterpolant(efit.psi1d/f.psi_norm,efit.fpol/(f.bnorm*f.xnorm),'cubic'); % possible postive or negative
f.psi2drz = efit.psi2drz/f.psi_norm;
f.fpol2drz = f.fpol2d(f.psi2drz);  % with sign dependence (positive or negative).

f.zdim = efit.zdim/f.xnorm;
f.zmid = efit.zmid/f.xnorm;
f.rleft = efit.rleft/f.xnorm;

f.rbbbs = efit.rbbbs/f.xnorm;
f.zbbbs = efit.zbbbs/f.xnorm;

f.delR = (f.R(end)-f.R(1))/(f.nr-1);
f.delZ = (f.Z(end)-f.Z(1))/(f.nz-1);

%% Psi spline
f.psisp = zeros(9,f.nz,f.nr);
f.psisp(1,:,:) = f.psi2drz;
f.psisp = qdspline.constructSpline2d(0, 0, f.nz, f.nr, f.delZ, f.delR, f.psisp);

%% fpol spline
f.fpolsp = zeros(9,f.nz,f.nr);
f.fpolsp(1,:,:) = f.fpol2drz;
f.fpolsp = qdspline.constructSpline2d(0, 0, f.nz, f.nr, f.delZ, f.delR, f.fpolsp);

%% B0 & BT spline
f.dpsidZ = squeeze(f.psisp(2,:,:));
f.dpsidR = squeeze(f.psisp(4,:,:));

f.BR = cocos.RpZ_Ip*(f.dpsidZ./f.R);  % with sign dependence (positive or negative).
f.BZ = cocos.RpZ_Ip*(-f.dpsidR./f.R);  % with sign dependence (positive or negative).
f.BT = f.fpol2drz./f.R2drz;  % with sign dependence (positive or negative).
f.B0_amp = sqrt(f.BR.^2 + f.BZ.^2 + f.BT.^2);  % B0 amplitude, positive only.

f.b0_amp_sp = zeros(9,f.nz,f.nr);
f.b0_amp_sp(1,:,:) = f.B0_amp;
f.b0_amp_sp = qdspline.constructSpline2d(0, 0, f.nz, f.nr, f.delZ, f.delR, f.b0_amp_sp);

f.brsp = zeros(9,f.nz,f.nr);
f.brsp(1,:,:) = f.BR;
f.brsp = qdspline.constructSpline2d(0, 0, f.nz, f.nr, f.delZ, f.delR, f.brsp);

f.bzsp = zeros(9,f.nz,f.nr);
f.bzsp(1,:,:) = f.BZ;
f.bzsp = qdspline.constructSpline2d(0, 0, f.nz, f.nr, f.delZ, f.delR, f.bzsp);

f.btsp = zeros(9,f.nz,f.nr);
f.btsp(1,:,:) = f.BT;
f.btsp = qdspline.constructSpline2d(0, 0, f.nz, f.nr, f.delZ, f.delR, f.btsp);

%% safety guard
if f.xnorm < 0
    error(['Negative normalization values for length! ','f.xnorm is: ',num2str(f.xnorm)])
end

%%
iplot = 1;
if iplot == 1
    figure('name','phase_space',...
        'unit','normalized',...
        'position',[0.0,0.0,1,0.43],... % figure position
        'DefaultAxesFontSize',15,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',2,...
        'DefaultAxesTickLength',[0.02,0.05]);
    subplot(141)
    contourf(f.R2drz,f.Z2drz,f.BR,20);daspect([1 1 1]);colorbar;title('BR/Baxis (vector with sign)');hold on;
    plot(efit.rbbbs/f.xnorm,efit.zbbbs/f.xnorm,'k-','linewidth',2);hold on;
    plot(efit.rlim/f.xnorm,efit.zlim/f.xnorm,'r-','linewidth',2);hold on;
    plot(efit.rmaxis/f.xnorm,efit.zmaxis/f.xnorm,'go','markersize',5,'MarkerFaceColor','g');
    xlabel('R/R0')
    ylabel('Z/R0')
    
    subplot(142)
    contourf(f.R2drz,f.Z2drz,f.BZ,20);daspect([1 1 1]);colorbar;title('BZ/Baxis (vector with sign)');hold on;
    plot(efit.rbbbs/f.xnorm,efit.zbbbs/f.xnorm,'k-','linewidth',2);hold on;
    plot(efit.rlim/f.xnorm,efit.zlim/f.xnorm,'r-','linewidth',2);hold on;
    plot(efit.rmaxis/f.xnorm,efit.zmaxis/f.xnorm,'go','markersize',5,'MarkerFaceColor','g');
    xlabel('R/R0')
    ylabel('Z/R0')
    
    subplot(143)
    contourf(f.R2drz,f.Z2drz,f.BT,20);daspect([1 1 1]);colorbar;title('BT/Baxis (vector with sign)');hold on;
    plot(efit.rbbbs/f.xnorm,efit.zbbbs/f.xnorm,'k-','linewidth',2);hold on;
    plot(efit.rlim/f.xnorm,efit.zlim/f.xnorm,'r-','linewidth',2);hold on;
    plot(efit.rmaxis/f.xnorm,efit.zmaxis/f.xnorm,'go','markersize',5,'MarkerFaceColor','g');
    xlabel('R/R0')
    ylabel('Z/R0')
    
    subplot(144)
    contourf(f.R2drz,f.Z2drz,f.B0_amp,20);daspect([1 1 1]);colorbar;title('B0\_amp/Baxis (amplitude > 0)');hold on;
    plot(efit.rbbbs/f.xnorm,efit.zbbbs/f.xnorm,'k-','linewidth',2);hold on;
    plot(efit.rlim/f.xnorm,efit.zlim/f.xnorm,'r-','linewidth',2);hold on;
    plot(efit.rmaxis/f.xnorm,efit.zmaxis/f.xnorm,'go','markersize',5,'MarkerFaceColor','g');
    xlabel('R/R0')
    ylabel('Z/R0')
    
end
