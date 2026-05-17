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
% Extract and organize data: m values and corresponding single particle (No.24) push execution times
m_values = [0, 1, 3, 5, 7];
push_times = [0.0236, 0.0732, 0.114, 0.198, 0.224];

% Create a new figure window
figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.45,0.59],... % figure position
        'DefaultAxesFontSize',28,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',1,...
        'DefaultAxesTickLength',[0.005,0.005]);
% Plot the bar chart
b = bar(m_values, push_times);

% Beautify the bar chart color (using a soft blue)
b.FaceColor = [0.2 0.6 0.8]; 
b.EdgeColor = 'none'; % Remove edge lines for a more modern look

% Set X-axis ticks to ensure only the provided m values are displayed
xticks(m_values);

% Add axis labels and title
xlabel('Number of m in the Perturbation Field', 'FontSize', 30, 'FontWeight', 'bold');
ylabel('Simulation Time (s)', 'FontSize', 30, 'FontWeight', 'bold');%title('Comparison of Single Particle Push Execution Speed for Different m Values', 'FontSize', 14, 'FontWeight', 'bold');

% Add specific value labels on top of each bar (applicable for MATLAB R2019b and later)
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(push_times);
text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize',22, 'FontWeight', 'bold');

% Adjust Y-axis limits to leave space for the top labels (increase height by 15%)
ylim([0, max(push_times) * 1.15]);

% Turn on Y-axis grid lines for easier height comparison
grid on;
ax = gca;
ax.YGrid = 'on';
ax.XGrid = 'off'; % No grid needed for the X-axis
ax.GridLineStyle = '--';
ax.GridAlpha = 0.4;
   %saveas(gca,'../output/CoM_2','png');
    print(gcf, '-depsc', '-r600', ['../output/','time_compare','.eps']);