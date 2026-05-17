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
global efit cocos fielddir
efit=readEFITgfile([datapath,efitfile]);

if efit.simag > efit.sibry
    efit.simag = - efit.simag;
    efit.sibry = - efit.sibry;
    efit.psizr = - efit.psizr;
end

efit.mpsi = size(efit.qpsi,1);
efit.psiw = efit.sibry - efit.simag;
efit.psi1d = linspace(0,efit.psiw,efit.mpsi)';
efit.dpsi = efit.psiw/(efit.mpsi-1);

efit.eqLdp = operator2.gradient1d(3, efit.psi1d);
efit.ffprim = efit.fpol.*(efit.eqLdp*efit.fpol);
efit.pprime = (efit.eqLdp*efit.pres);
run efit_output;

q_tmp=zeros(efit.mpsi-1,1);
for i = 2:efit.mpsi
    q_tmp(i-1) = (efit.qpsi(i-1)+efit.qpsi(i))*0.5;
end
torpsi_tmp = cumsum(q_tmp*efit.dpsi); 
efit.torpsi = [0;torpsi_tmp];

efit.nr = efit.nw;
efit.nz = efit.nh;
efit.R = linspace(efit.rleft,efit.rleft+efit.rdim,efit.nr);
efit.delR = efit.rdim/(efit.nr-1);
efit.Z = linspace(efit.zmid-0.5*efit.zdim,efit.zmid+0.5*efit.zdim,efit.nz);
efit.delZ = efit.zdim/(efit.nz-1);

%%
% R determines the number of columns, Z determines the number of rows, namely, z is the first
% dimension, r is the second dimension
[efit.R2drz,efit.Z2drz] = meshgrid(efit.R,efit.Z);
% efit.psizr: r is the first dimension, z is the second dimension
% efit.psi2drz: z is the first dimension, r is the second dimension
efit.psi2drz = efit.psizr' - efit.simag;

efit.psisp = zeros(9,efit.nz,efit.nr);
efit.psisp(1,:,:) = efit.psi2drz;
efit.psisp = qdspline.constructSpline2d(0, 0, efit.nz, efit.nr, efit.delZ, efit.delR, efit.psisp);
efit.qsp = zeros(3,efit.mpsi);
efit.qsp(1,:) = efit.qpsi;
efit.qsp = qdspline.constructSpline1d(0, efit.mpsi, efit.dpsi, efit.qsp);
efit.fsp = zeros(3,efit.mpsi);
efit.fsp(1,:) = efit.fpol;
efit.fsp = qdspline.constructSpline1d(0, efit.mpsi, efit.dpsi, efit.fsp);

[~,efit.raxis_ind] = min(abs(efit.R - efit.rmaxis));
[~,efit.zaxis_ind] = min(abs(efit.Z - efit.zmaxis));

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Give coordinate directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get sign of Bt and Ip with respect to counter-clockwise(C.C.W.) phi in EFIT
cocos.RpZ_Ip = sign(efit.current);
cocos.RpZ_Bt = sign(efit.bcentr);

% Psi positive is always from axis to edge.
if cocos.RpZ_Bt < 0 && cocos.RpZ_Ip < 0
    
    % Ip -> out, Bt -> out.
    % Theta positive is along counter-clockwise, zeta positive is along
    % clockwise.
    fielddir = 'zero';
    
elseif cocos.RpZ_Bt > 0 && cocos.RpZ_Ip < 0
    
    % Ip -> out, Bt -> in.
    % Theta positive is along counter-clockwise, zeta positive is along
    % clockwise.
    fielddir = 'one';
    
elseif cocos.RpZ_Bt > 0 && cocos.RpZ_Ip > 0
    
    % Ip -> in, Bt -> in.
    % Theta positive is along clockwise, zeta positive is along
    % counter-clockwise.
    fielddir = 'two';
    
elseif cocos.RpZ_Bt < 0 && cocos.RpZ_Ip > 0
    
    % Ip -> in, Bt -> out.
    % Theta positive is along clockwise, zeta positive is along
    % counter-clockwise.
    fielddir = 'three';
    
end
cocos.rhotp_Bt = cocos.RpZ_Bt*cocos.RpZ_Ip;
%% check efit equilibrium
figure('name','EFIT psi');
contourf(efit.R2drz,efit.Z2drz,efit.psi2drz,100,'edgecolor','none');colorbar;
hold on;
contour(efit.R2drz,efit.Z2drz,efit.psi2drz,[efit.psiw efit.psiw],'linecolor','k','linewidth',2);hold on;
daspect([1 1 1]);
plot(efit.rmaxis,efit.zmaxis,'r*');
xlabel('R (m)','fontsize',16)
ylabel('Z (m)','fontsize',16)
title('EFIT poloidal magnetic flux')
