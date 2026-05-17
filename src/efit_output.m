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
figure('name','efit_1d_data_check',...
    'unit','normalized',...
    'Position',[0.1 0.1 1 1],...
    'DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','normal',...
    'DefaultAxesLineWidth',2,...
    'DefaultAxesTickLength',[0.02,0.05]);
subplot(221)
plot(efit.psi1d,efit.pres,'r.-');
xlabel('$\psi (Weber/rad)$','interpreter','latex');
ylabel('$P(Pa)$','interpreter','latex');
set(gca,'fontsize',20);
grid on;

subplot(222)
plot(efit.psi1d,efit.fpol,'r.-');
xlabel('$\psi (Weber/rad)$','interpreter','latex');
ylabel('$F$','interpreter','latex');
set(gca,'fontsize',20);
grid on;

subplot(223)
plot(efit.psi1d,efit.pprime_raw,'r.-');
hold on;
plot(efit.psi1d,efit.pprime,'b');
legend('pprime','dpres/dp');
xlabel('$\psi (Weber/rad)$','interpreter','latex');
ylabel('$\partial P/\partial\psi$','interpreter','latex');
title('consistency check: pprime and pres');
set(gca,'fontsize',20);
grid on;

subplot(224)
plot(efit.psi1d,efit.ffprim_raw,'r.-');
hold on;
plot(efit.psi1d,efit.ffprim,'b');
legend('ffprim','F*dF/dp');
xlabel('$\psi (Weber/rad)$','interpreter','latex');
ylabel('$F\partial F/\partial\psi$','interpreter','latex');
title('consistency check: ffprim and fpol');
set(gca,'fontsize',20);
grid on;
if ~exist('./output','dir');
    mkdir('output');
end
saveas(gca,'./output/efit_1d','png')

figure('name','efit_q_check',...
    'unit','normalized',...
    'DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','normal',...
    'DefaultAxesLineWidth',2,...
    'DefaultAxesTickLength',[0.02,0.05]);
plot(efit.psi1d,efit.qpsi,'r.-');
xlabel('$\psi (Weber/rad)$','interpreter','latex');
ylabel('$q$','interpreter','latex');
set(gca,'fontsize',20);
set(gca,'ytick',0:1:20);
grid on;
if ~exist('./output','dir');
    mkdir('output');
end
saveas(gca,'./output/q','png')