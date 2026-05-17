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
 figure('name','region',...
        'unit','normalized',...
        'position',[0.0,0.0,0.25,0.4],... % figure position
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',1.3,...
        'DefaultAxesTickLength',[0.008,0.008]);

load(['../input/','data3d_' ps_option '.mat']);

run plot_GTC.m

passing_option = 'counter-passing';
filename = [ps_option,'_',passing_option, '.mat'];
save_path = fullfile('../save/', filename);
load(save_path);
char_freq1=char_freq;
run plot_MAS.m

passing_option = 'co-passing'; 
filename = [ps_option,'_',passing_option, '.mat'];
save_path = fullfile('../save/', filename);
load(save_path);

run plot_MAS.m
run legend_paper.m
if strcmp(ps_option,'MG')
    run plot_resonan_no.m
end
saveas(gca,['../output/' 'co&counter' ps_option],'epsc');
