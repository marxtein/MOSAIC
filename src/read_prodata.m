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
function A=read_prodata(path,file)
%read in profile.dat
filename = [path,file];
data = importdata(filename,' ',1);

A.psi = data.data(:,1);
A.x = data.data(:,2);
A.r = data.data(:,3);
A.R = data.data(:,4);
A.Rr = data.data(:,5);
A.Te = data.data(:,6);
A.ne = data.data(:,7);
A.Ti = data.data(:,8);
A.Zeff = data.data(:,9);
A.omega_tor = data.data(:,10);
A.Er = data.data(:,11);
A.ni = data.data(:,12);
A.nimp = data.data(:,13);
A.nf = data.data(:,14);
A.Tf = data.data(:,15);
