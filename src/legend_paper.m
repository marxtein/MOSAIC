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
%---------------------------------------------------------------------------------------------图例
if strcmp(ps_option,'MG')
    posx = min(xlim) + 0.03*range(xlim);
    posy = min(ylim) + 0.79*range(ylim);
elseif strcmp(ps_option,'GB')
    posx = min(xlim) + 0.03*range(xlim);
    posy = min(ylim) + 0.79*range(ylim);
elseif strcmp(ps_option,'PLam')
    posx = min(xlim) + 0.03*range(xlim);
    posy = min(ylim) + 0.79*range(ylim);
end



    w = 0.45*range(xlim);
    h = 0.19*range(ylim);
    fill([posx posx+w posx+w posx],[posy-h/2 posy-h/2 posy+h posy+h],'w',...
        'EdgeColor','k','LineWidth',1.0);


    x0 = posx + 0.08*range(xlim);
    y0 = posy + 0.14*range(ylim);
    a = 0.03 * range(xlim);
    b = 0.03 * range(ylim);
    N = 200; theta = linspace(0,2*pi,N);
    cmap = colormap('jet'); nC = size(cmap,1);
    for k = 1:nC
        r = 1 - (k-1)/(nC-1);
        fill(x0 + a*r*cos(theta), y0-0.015 + b*r*sin(theta), cmap(k,:), ...
            'EdgeColor','none');
    end
    plot(x0 + a*cos(theta), y0-0.015 + b*sin(theta), 'k', 'LineWidth', 1.2);


   

    text(x0 + 0.056*range(xlim), y0-0.02, '$\delta \hat{f}$(GTC)', ...
        'FontSize',23,'FontName','Times New Roman', ...
        'FontWeight','normal','Interpreter','latex', ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
    text(x0 - 0.025*range(xlim), y0-h/1.9, '$\cdot\cdot \ p \in\mathcal{Z}$ ', ...
        'FontSize',23,'FontName','Times New Roman', ...
        'FontWeight','normal','Interpreter','latex', ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
    text(x0 - 0.02*range(xlim), y0-h/1, '- -$\langle q \rangle$ = 3', ...
        'FontSize',23,'FontName','Times New Roman', ...
        'FontWeight','normal','Interpreter','latex', ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
    %--------------------------------------------------------------------------------------------------
