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
% 三条曲线
xA = PLam_2D.Pzeta(1,:);   yA = PLam_2D.lambda(1,:);
xB = PLam_2D.Pzeta(3,:);   yB = 0*PLam_2D.lambda(3,:);
xC = PLam_2D.Pzeta(3,:);   yC = PLam_2D.lambda(3,:);

% 统一插值参数
N = 3000; % 插值密度（越大越平滑）
tA = linspace(0,1,length(xA));
tB = linspace(0,1,length(xB));
tC = linspace(0,1,length(xC));

ti = linspace(0,1,N);

xA_i = interp1(tA, xA, ti, 'pchip');
yA_i = interp1(tA, yA, ti, 'pchip');

xB_i = interp1(tB, xB, ti, 'pchip');
yB_i = interp1(tB, yB, ti, 'pchip');

xC_i = interp1(tC, xC, ti, 'pchip');
yC_i = interp1(tC, yC, ti, 'pchip');

% ---------- 自动求三条曲线之间的交点 ----------
[xAB, yAB] = polyxpoly(xA_i, yA_i, xB_i, yB_i);
[xBC, yBC] = polyxpoly(xB_i, yB_i, xC_i, yC_i);
%xBC=xBC(2);
%yBC=yBC(2);
xBC=xC_i(end);
yBC=yC_i(end);
[xCA, yCA] = polyxpoly(xC_i, yC_i, xA_i, yA_i);

findNearestIndex = @(xCurve,yCurve,xI,yI) ...
    find( (xCurve - xI).^2 + (yCurve - yI).^2 == ...
          min( (xCurve - xI).^2 + (yCurve - yI).^2 ), 1);
% ---------- 找到交点在各曲线上的索引（用于截取区间） ----------
% A ∩ B
iAB_A = findNearestIndex(xA_i, yA_i, xAB, yAB);
iAB_B = findNearestIndex(xB_i, yB_i, xAB, yAB);

% B ∩ C
iBC_B = findNearestIndex(xB_i, yB_i, xBC, yBC);
iBC_C = findNearestIndex(xC_i, yC_i, xBC, yBC);

% C ∩ A
iCA_C = findNearestIndex(xC_i, yC_i, xCA, yCA);
iCA_A = findNearestIndex(xA_i, yA_i, xCA, yCA);

% ---------- 按闭合方向拼接区域（逆时针） ----------

% A 段：从 CA → AB
xA_seg = xA_i(min(iCA_A,iAB_A):max(iCA_A,iAB_A));
yA_seg = yA_i(min(iCA_A,iAB_A):max(iCA_A,iAB_A));

% B 段：从 AB → BC
xB_seg = xB_i(min(iAB_B,iBC_B):max(iAB_B,iBC_B));
yB_seg = yB_i(min(iAB_B,iBC_B):max(iAB_B,iBC_B));

% C 段：从 BC → CA
xC_seg = xC_i(min(iBC_C,iCA_C):max(iBC_C,iCA_C));
yC_seg = yC_i((min(iBC_C,iCA_C):max(iBC_C,iCA_C)));

% ---------- 合并三个部分形成闭合区域 ----------
x_fill = [xA_seg, xB_seg, xC_seg];
y_fill = [yA_seg, yB_seg, yC_seg];


