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
function A=read_spdata(path,file,option)


%% read in A.dat
filename = [path,file];
fid = fopen(filename,'r');
tline = fgetl(fid);
tline = fgetl(fid);
np = str2num(tline);
A.lsp = np(1);
A.lst = np(2);
A.lemax = np(3);
A.lrmax = np(4);
tline = fgetl(fid);
np = str2num(tline);
A.psiw = np(1);
A.ped = np(2);

spdim_2d = 9; % spline dimension for 2d array
spdim_1d = 3; % spline dimension for 1d array
num = A.lst*spdim_2d*4+spdim_1d*6;

for i = 1:A.lsp
    data1 = fscanf(fid,'%f',num);
    nindex_tmp = 0;
    for s = 1:spdim_2d
        for j = 1:A.lst
            sj = (s-1)*A.lst+j;
            A.bsp(s,i,j) = data1(sj);  % magnetic field
        end
    end
    nindex_tmp = nindex_tmp+sj;
    for s = 1:spdim_2d
        for j = 1:A.lst
            sj = (s-1)*A.lst+j;
            A.xsp(s,i,j) = data1(sj+nindex_tmp);  % X position
        end
    end
    nindex_tmp = nindex_tmp+sj;
    for s = 1:spdim_2d
        for j = 1:A.lst
            sj = (s-1)*A.lst+j;
            A.zsp(s,i,j) = data1(sj+nindex_tmp);  % Z position
        end
    end
    nindex_tmp = nindex_tmp+sj;
    for s = 1:spdim_2d
        for j = 1:A.lst
            sj = (s-1)*A.lst+j;
            A.gsp(s,i,j) = data1(sj+nindex_tmp);  % jacobian
        end
    end
    nindex_tmp = nindex_tmp+sj;
    
    A.qpsi(1:spdim_1d,i) = data1(1+nindex_tmp:spdim_1d+nindex_tmp);   % safety facotor
    nindex_tmp = nindex_tmp+spdim_1d;
    
    A.gpsi(1:spdim_1d,i) = data1(1+nindex_tmp:spdim_1d+nindex_tmp);   % current
    nindex_tmp = nindex_tmp+spdim_1d;
    
    A.ipsi(1:spdim_1d,i) = data1(1+nindex_tmp:spdim_1d+nindex_tmp);   % current
    nindex_tmp = nindex_tmp+spdim_1d;
    
    A.ppsi(1:spdim_1d,i) = data1(1+nindex_tmp:spdim_1d+nindex_tmp);   % pressure
    nindex_tmp = nindex_tmp+spdim_1d;
    
    A.rpsi(1:spdim_1d,i) = data1(1+nindex_tmp:spdim_1d+nindex_tmp);   % radius at outer-midplane
    nindex_tmp = nindex_tmp+spdim_1d;
    
    A.torpsi(1:spdim_1d,i) = data1(1+nindex_tmp:spdim_1d+nindex_tmp); % toroidal psi
    nindex_tmp = nindex_tmp+spdim_1d;
end

data2 = fscanf(fid,'%f',7);
A.krip = data2(1);
A.nrip = data2(2);
A.rmaj = data2(3);
A.d0 = data2(4);
A.brip = data2(5);
A.wrip = data2(6);
A.xrip = data2(7);
if strcmp(option,'mapping')
% read nu, phi = zeta_B + nu
num3 = A.lst*spdim_2d;
for i = 1:A.lsp
    data3 = fscanf(fid,'%f',num3);
    for s = 1:spdim_2d
        for j = 1:A.lst
            sj = (s-1)*A.lst+j;
            A.nsp(s,i,j) = data3(sj);
        end
    end
end

% read poloidal flux psi
num4 = A.lsp;
data4 = fscanf(fid,'%f',num4);
for i = 1:A.lsp
    A.psi(i) = data4(i);
end

% read magnetic field information at axis and separatrix
data5 = fscanf(fid,'%f',2);
A.torped = data5(1);
A.baxis = data5(2);

num6 = A.lst*spdim_2d;
for i = 1:A.lsp
    data6 = fscanf(fid,'%f',num6);
    for s = 1:spdim_2d
        for j = 1:A.lst
            sj = (s-1)*A.lst+j;
            A.delsp(s,i,j) = data6(sj);
        end
    end
end

num7 = A.lst*spdim_2d;
for i = 1:A.lsp
    data7 = fscanf(fid,'%f',num7);
    for s = 1:spdim_2d
        for j = 1:A.lst
            sj = (s-1)*A.lst+j;
            A.jsp(s,i,j) = data7(sj);
        end
    end
end
elseif strcmp(option,'orbit')

A.psi = linspace(0,A.psiw,A.lsp);
A.torped = A.torpsi(1,A.lsp);
A.baxis = A.bsp(1,1,1);
A.rmaj = A.xsp(1,1,1)*100;
A.nsp = zeros(9,A.lsp,A.lst);
A.delsp = zeros(9,A.lsp,A.lst);
A.jsp = zeros(9,A.lsp,A.lst);
end
% process the data structure with poloidal periodicity
A.lst = A.lst+1; % orbit code convention is (2pi/lst), here we consider the offset for GTC convention (2pi/(lst-1))
A.bsp(1:9,1:A.lsp,A.lst) = A.bsp(1:9,1:A.lsp,1);
A.xsp(1:9,1:A.lsp,A.lst) = A.xsp(1:9,1:A.lsp,1);
A.zsp(1:9,1:A.lsp,A.lst) = A.zsp(1:9,1:A.lsp,1);
A.gsp(1:9,1:A.lsp,A.lst) = A.gsp(1:9,1:A.lsp,1);
A.nsp(1:9,1:A.lsp,A.lst) = A.nsp(1:9,1:A.lsp,1);
A.delsp(1:9,1:A.lsp,A.lst) = A.delsp(1:9,1:A.lsp,1);
A.jsp(1:9,1:A.lsp,A.lst) = A.jsp(1:9,1:A.lsp,1);

fclose(fid);

% %% plot magnetic field
% close all;
   figure('unit','normalized',...
        'Position',[0.0 0.0 0.3 0.5],...
        'DefaultAxesFontSize',20,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',1,...
        'DefaultAxesTickLength',[0.01,0.01]);
[C,h]=contourf(squeeze(A.xsp(1,:,:)),squeeze(A.zsp(1,:,:)),squeeze(A.bsp(1,:,:)),100);
set(h,'linecolor','none');
hold on;
step = 4;      % 稀疏度
lw = 0.4; 
x = squeeze(A.xsp(1,:,:));
z = squeeze(A.zsp(1,:,:));

% 【正确画法】只画网格线，不断线、自动闭合、最稳定
hMesh =plot(x(1:step-1:end,:)', z(1:step-1:end,:)', 'k-', 'LineWidth', lw); hold on  % 横线
plot(x(:,1:step:end),  z(:,1:step:end),  'k-', 'LineWidth', lw);   % 竖线
%plot(squeeze(A.xsp(1,:,:)),squeeze(A.zsp(1,:,:)),'k.-');
%hold on;
%plot(squeeze(A.xsp(1,:,:))',squeeze(A.zsp(1,:,:))','k.-');
colorbar;
colorbar;
cb = colorbar;
% 给色条加标签，就是B0/Ba
cb.Label.String = '$B_0$(T)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 20;
cb.Label.FontName = 'Times New Roman';
daspect([1 1 1]);
%title('Equilibrium magnetic field B_0(T)','fontsize',20);
xlabel('$X$/m','fontsize',25,'Interpreter','latex');
ylabel('$Z$/m','fontsize',25,'Interpreter','latex');

% ===================== 【正确图例，绝对不报错】 =====================
% 只画网格图例，colorbar 本身就代表 B0/Ba
legend(hMesh, 'Mesh Grid', ...
    'Interpreter','latex','Location','best','FontSize',23);
saveas(gca,'../output/B_mesh','epsc');
