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
% Three curves
xA = PLam_2D.Pzeta(1,:);   yA = PLam_2D.lambda(1,:);
xB = PLam_2D.Pzeta(4,:);   yB = PLam_2D.lambda(4,:);
xC = PLam_2D.Pzeta(5,:);   yC = PLam_2D.lambda(5,:);

% Unified interpolation parameters
N = 2000; % Interpolation density (higher = smoother)
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

% ---------- Automatically find intersection points between three curves ----------
[xAB, yAB] = polyxpoly(xA_i, yA_i, xB_i, yB_i);
[xBC, yBC] = polyxpoly(xB_i, yB_i, xC_i, yC_i);
xBC=xC(end);
yBC=yC(end);
[xCA, yCA] = polyxpoly(xC_i, yC_i, xA_i, yA_i);
xCA=xCA(1);
yCA=yCA(1);
findNearestIndex = @(xCurve,yCurve,xI,yI) ...
    find( (xCurve - xI).^2 + (yCurve - yI).^2 == ...
          min( (xCurve - xI).^2 + (yCurve - yI).^2 ), 1);
% ---------- Find the indices of intersection points on each curve (for interval truncation) ----------
% A ∩ B
iAB_A = findNearestIndex(xA_i, yA_i, xAB, yAB);
iAB_B = findNearestIndex(xB_i, yB_i, xAB, yAB);

% B ∩ C
iBC_B = findNearestIndex(xB_i, yB_i, xBC, yBC);
iBC_C = findNearestIndex(xC_i, yC_i, xBC, yBC);

% C ∩ A
iCA_C = findNearestIndex(xC_i, yC_i, xCA, yCA);
iCA_A = findNearestIndex(xA_i, yA_i, xCA, yCA);

% ---------- Stitch regions along closed direction (counterclockwise) ----------

% Segment A: CA → AB
xA_seg = xA_i(min(iCA_A,iAB_A):max(iCA_A,iAB_A));
yA_seg = yA_i(min(iCA_A,iAB_A):max(iCA_A,iAB_A));

% Segment B: AB → BC
xB_seg = xB_i(min(iAB_B,iBC_B):max(iAB_B,iBC_B));
yB_seg = yB_i(min(iAB_B,iBC_B):max(iAB_B,iBC_B));

% Segment C: BC → CA
xC_seg = xC_i(min(iBC_C,iCA_C):max(iBC_C,iCA_C));
yC_seg = yC_i((min(iBC_C,iCA_C):max(iBC_C,iCA_C)));

% ---------- Merge three parts into a closed region ----------
x_fill = [xA_seg, xB_seg, xC_seg];
y_fill = [yA_seg, yB_seg, yC_seg];








% % ---------- Automatically find intersection points between three curves ----------
% [xAB,yAB] = polyxpoly(xA,yA,xB,yB);   % A ∩ B
% [xBC,yBC] = polyxpoly(xB,yB,xC,yC);   % B ∩ C
% xBC=xC(end);
% yBC=yC(end);
% [xCA,yCA] = polyxpoly(xC,yC,xA,yA);   % C ∩ A
% xCA=xCA(1);
% yCA=yCA(1);
% findNearestIndex = @(xCurve,yCurve,xI,yI) ...
%     find( (xCurve - xI).^2 + (yCurve - yI).^2 == ...
%           min( (xCurve - xI).^2 + (yCurve - yI).^2 ), 1);
% % ---------- Find the indices of intersection points on each curve (for interval truncation) ----------
% % A ∩ B
% iAB_A = findNearestIndex(xA, yA, xAB, yAB);
% iAB_B = findNearestIndex(xB, yB, xAB, yAB);
% 
% % B ∩ C
% iBC_B = findNearestIndex(xB, yB, xBC, yBC);
% iBC_C = findNearestIndex(xC, yC, xBC, yBC);
% 
% % C ∩ A
% iCA_C = findNearestIndex(xC, yC, xCA, yCA);
% iCA_A = findNearestIndex(xA, yA, xCA, yCA);
% 
% % ---------- Stitch regions along closed direction (counterclockwise) ----------
% 
% % Segment A: CA → AB
% xA_seg = xA(min(iCA_A,iAB_A):max(iCA_A,iAB_A));
% yA_seg = yA(min(iCA_A,iAB_A):max(iCA_A,iAB_A));
% 
% % Segment B: AB → BC
% xB_seg = xB(min(iAB_B,iBC_B):max(iAB_B,iBC_B));
% yB_seg = yB(min(iAB_B,iBC_B):max(iAB_B,iBC_B));
% 
% % Segment C: BC → CA
% xC_seg = xC(min(iBC_C,iCA_C):max(iBC_C,iCA_C));
% yC_seg = yC((min(iBC_C,iCA_C):max(iBC_C,iCA_C)));
% 
% % ---------- Merge three parts into a closed region ----------
% x_fill = [xA_seg, xB_seg, xC_seg];
% y_fill = [yA_seg, yB_seg, yC_seg];


