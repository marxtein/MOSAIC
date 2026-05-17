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
function A = find_Pzeta_right(qpart,apart,E,lam0,lam1,lam_grid_num,sign_gpsi)
global eq unit

E_gtc_norm = E*1000/unit.energy_norm; % convert to GTC normalization

A = struct;
A.lam = linspace(lam0,lam1,lam_grid_num);

A.psi_ind = zeros(1,lam_grid_num);
A.psi = zeros(1,lam_grid_num);
A.gpsi = zeros(1,lam_grid_num);
A.b0 = zeros(1,lam_grid_num);
A.rhopara = zeros(1,lam_grid_num);
A.Pzeta_right = zeros(1,lam_grid_num);


mpsi = 1000;
psi_dum = linspace(0,eq.psiw,mpsi);
term1 = zeros(1,mpsi);
term2 = zeros(1,mpsi);

for j = 1:lam_grid_num
    lam_dum = A.lam(j);
    for i = 1:mpsi
        b0_tmp =  qdspline.spline2d(0, psi_dum(i), 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
        db0_tmp =  qdspline.spline2d(1, psi_dum(i), 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
        gpsi_tmp = qdspline.spline1d(0, psi_dum(i), eq.lsp, eq.dpsi, eq.gpsi);
        dgpsi_tmp = qdspline.spline1d(1, psi_dum(i), eq.lsp, eq.dpsi, eq.gpsi);
        
        if strcmp(sign_gpsi,'negative')
            rhopara = - sqrt(2*E_gtc_norm*(1-lam_dum*b0_tmp)/b0_tmp^2*apart/qpart^2);
        elseif strcmp(sign_gpsi,'positive')
            rhopara = sqrt(2*E_gtc_norm*(1-lam_dum*b0_tmp)/b0_tmp^2*apart/qpart^2);
        else
            error('Wrong ''sign_gpsi'' string. Set ''negative'' or ''positive''.')
        end
        
        term1(i) = - 1 + rhopara*dgpsi_tmp;
        term2(i) = - gpsi_tmp/rhopara*(apart/qpart^2/b0_tmp^2)*(-2*E_gtc_norm/b0_tmp + lam_dum*E_gtc_norm)*db0_tmp;
    end
    
    %plot(psi_dum,term1,'r');hold on;
    %plot(psi_dum,term2,'b')
    %xlim([0.01*psi_dum(end),psi_dum(end)]);
    
    [~,A.psi_ind(j)]=min(abs(term1 - term2));
    A.psi(j) = psi_dum(A.psi_ind(j));
    A.gpsi(j) = qdspline.spline1d(0, A.psi(j), eq.lsp, eq.dpsi, eq.gpsi);
    A.b0(j) =  qdspline.spline2d(0, A.psi(j), 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
    A.rhopara(j) = -sqrt(2*E_gtc_norm*(1-lam_dum*A.b0(j))/A.b0(j)^2*apart/qpart^2);
    A.Pzeta_right(j) = (A.gpsi(j)*A.rhopara(j) - A.psi(j))/eq.psiw;
    
end
end