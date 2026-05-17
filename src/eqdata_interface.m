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
function  f = eqdata_interface

global spdata

%% normalization
f.lsp = spdata.lsp;
f.lst = spdata.lst;

% normalization unit
f.bnorm = spdata.baxis;
f.xnorm = spdata.rmaj/100;
f.fx = f.bnorm*f.xnorm*f.xnorm;
f.psiw = spdata.psiw/f.fx;
f.ped = spdata.ped/f.fx;
f.torped = spdata.torped/f.fx;
f.psi = spdata.psi/f.fx;
f.dpsi = (spdata.psi(end)-spdata.psi(1))/f.fx/(spdata.lsp-1);
f.dtheta = 2.0*pi/(f.lst-1);
f.theta = 0.0:f.dtheta:2.0*pi;

%% 1d f array normalization
f.torpsi = sp_norm.normalizeSpline1d(spdata.torpsi, f.fx, f.fx);

f.qpsi = sp_norm.normalizeSpline1d(spdata.qpsi, 1.0, f.fx);

f.gpsi = sp_norm.normalizeSpline1d(spdata.gpsi, f.xnorm*f.bnorm, f.fx);

f.ipsi = sp_norm.normalizeSpline1d(spdata.ipsi, f.xnorm*f.bnorm, f.fx);

f.rpsi = sp_norm.normalizeSpline1d(spdata.rpsi, 1.0, f.fx); % already normalized

f.ppsi_norm = (1.602e-19*1e6);
f.ppsi = sp_norm.normalizeSpline1d(spdata.ppsi, f.ppsi_norm, f.fx); % from Pascal to cm^-3*eV
%% 2d f array normalization
f.bsp = sp_norm.normalizeSpline2d(spdata.bsp, f.bnorm, f.fx);

f.zsp = sp_norm.normalizeSpline2d(spdata.zsp, f.xnorm, f.fx);

f.xsp = sp_norm.normalizeSpline2d(spdata.xsp, f.xnorm, f.fx);

f.nsp = sp_norm.normalizeSpline2d(spdata.nsp, 1, f.fx);

f.gsp = sp_norm.normalizeSpline2d(spdata.gsp, f.xnorm/f.bnorm, f.fx); % not used

f.delsp = sp_norm.normalizeSpline2d(spdata.delsp, 1/f.xnorm, f.fx);

f.jsp_norm = 1.0/(f.xnorm*4*pi*1e-7);
f.jsp = sp_norm.normalizeSpline2d(spdata.jsp, f.jsp_norm, f.fx);

% use quadratic spline to re-calculate the spline coefficient, even though
% the derivatives are already given from spdata file
f.bsp = qdspline.constructSpline2d(0, 0, f.lsp, f.lst, f.dpsi, f.dtheta, f.bsp);
f.zsp = qdspline.constructSpline2d(0, 0, f.lsp, f.lst, f.dpsi, f.dtheta, f.zsp);
f.xsp = qdspline.constructSpline2d(0, 0, f.lsp, f.lst, f.dpsi, f.dtheta, f.xsp);
f.gsp = qdspline.constructSpline2d(0, 0, f.lsp, f.lst, f.dpsi, f.dtheta, f.gsp);

f.qpsi = qdspline.constructSpline1d(0, f.lsp, f.dpsi, f.qpsi);
f.gpsi = qdspline.constructSpline1d(0, f.lsp, f.dpsi, f.gpsi);
f.ipsi = qdspline.constructSpline1d(0, f.lsp, f.dpsi, f.ipsi);
f.rpsi = qdspline.constructSpline1d(0, f.lsp, f.dpsi, f.rpsi);
f.torpsi = qdspline.constructSpline1d(0, f.lsp, f.dpsi, f.torpsi);
f.ppsi = qdspline.constructSpline1d(0, f.lsp, f.dpsi, f.ppsi);
%% plot magnetic field
%    figure('unit','normalized',...
%         'Position',[0.0 0.0 0.3 0.5],...
%         'DefaultAxesFontSize',20,...
%         'DefaultAxesFontWeight','normal',...
%         'DefaultAxesLineWidth',1,...
%         'DefaultAxesTickLength',[0.01,0.01]);
% 
% % 画磁场等高线
% [C,hContour] = contourf(squeeze(f.xsp(1,:,:)*R0),squeeze(f.zsp(1,:,:)*R0),squeeze(f.bsp(1,:,:)),100);
% set(hContour,'linecolor','none');
% hold on;
% 
% step = 5;      % 稀疏度
% lw = 0.4;      % 线条粗细
% size(f.xsp)
% % 画稀疏网格线
% %hMesh = plot(squeeze(f.xsp(1,1:step:end,1:step:end)), squeeze(f.zsp(1,1:step:end,1:step:end)), 'k-', 'LineWidth', lw);
% %plot(squeeze(f.xsp(1,1:step:end,1:step:end))', squeeze(f.zsp(1,1:step:end,1:step:end))', 'k-', 'LineWidth', lw);
% % 取出网格（只取一次，干净）
% x = squeeze(f.xsp(1,:,:))*R0;
% z = squeeze(f.zsp(1,:,:))*R0;
% 
% % 【正确画法】只画网格线，不断线、自动闭合、最稳定
% hMesh =plot(x(1:step-1:end,:)', z(1:step-1:end,:)', 'k-', 'LineWidth', lw);   % 横线
% plot(x(:,1:step:end),  z(:,1:step:end),  'k-', 'LineWidth', lw);   % 竖线
% % 坐标轴设置
% colorbar;
% cb = colorbar;
% % 给色条加标签，就是B0/Ba
% cb.Label.String = '$B_0/B_a$';
% cb.Label.Interpreter = 'latex';
% cb.Label.FontSize = 20;
% cb.Label.FontName = 'Times New Roman';
% %cb.FontWeight = 'bold';          
% cb.LineWidth = 1; 
% daspect([1 1 1]);
% %title('Equilibrium magnetic field B_0/B_a','fontsize',20);
% xlabel('$X$/m','fontsize',25,'Interpreter','latex');
% ylabel('$Z$/m','fontsize',25,'Interpreter','latex');
% 
% % ===================== 【正确图例，绝对不报错】 =====================
% % 只画网格图例，colorbar 本身就代表 B0/Ba
% legend(hMesh, 'Mesh Grid', ...
%     'Interpreter','latex','Location','best','FontSize',23);
% saveas(gca,'../output/B_mesh','epsc');

end