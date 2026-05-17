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
% 起点 (x0, y0)
x0 = 0;
y0 = 30;

% 斜率 m
m = 25.5556;

% 定义线的范围
x = linspace(-0.1, 0.2,100);
y = y0 + m * (x - x0);

% 绘制线
plot(x, y, 'b-', 'LineWidth', 2);
hold on;

% 在线旁边添加序号
text(x(end)+0.08, y(end), '(1)', 'FontSize', 20, 'Color', 'b', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');


x0 = 0;
y0 = 45.09;

% 斜率 m
m = 25.5556;

% 定义线的范围
x = linspace(0.1, 0.4,100);
y = y0 + m * (x - x0);

% 绘制线
plot(x, y, 'b-', 'LineWidth', 2);
hold on;

% 在线旁边添加序号
text(x(end)+0.1, y(end), '(2)', 'FontSize', 20, 'Color', 'b', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');


x0 = 0;
y0 = 51.93;

% 斜率 m
m = 25.5556;

% 定义线的范围
x = linspace(0.1, 0.4,100);
y = y0 + m * (x - x0);

% 绘制线
plot(x, y, 'b-', 'LineWidth', 2);
hold on;

% 在线旁边添加序号
text(x(end)+0.1, y(end), '(3)', 'FontSize', 20, 'Color', 'b', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

x0 = 0;
y0 = 60.22;

% 斜率 m
m = 25.5556;

% 定义线的范围
x = linspace(0.1, 0.4,100);
y = y0 + m * (x - x0);

% 绘制线
plot(x, y, 'b-', 'LineWidth', 2);
hold on;

% 在线旁边添加序号
text(x(end)+0.1, y(end), '(4)', 'FontSize', 20, 'Color', 'b', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');


