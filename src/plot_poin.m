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
        'position',[0.0,0.0,0.3,0.45],... % figure position
        'DefaultAxesFontSize',23,...
        'DefaultAxesFontWeight','normal',...
        'DefaultAxesLineWidth',1,...
        'DefaultAxesTickLength',[0.02,0.015]);
%for m=25
for m=1:size(Poini.psi,2)-num_ones
    poin_size=5;
x = Poini.E(:,m);
%x = Poini.p_zeta(:,m);
%y = (Poini.theta(:,m));
y=mod(Poini.theta(:,m)+pi, 2*pi) - pi;
x(find(x==0))=nan;
y(find(y==0))=nan;
validIdx = ~isnan(x) & ~isnan(y);
x= x(validIdx);
y= y(validIdx);

scatter(x, y/pi,poin_size, 'filled');hold on
%scatter(x, y,poin_size, 'filled');hold on

end
axis tight
xlim([PoinE_down,PoinE_up]);
%xlim([-0.1,0.13]);
ylim([-1,1]);
%xlim([26,35]);
box on
xlabel('$E$(keV)','fontsize',33,'Interpreter','latex')
ylabel('$\theta/\pi$','fontsize',33,'Interpreter','latex')
saveas(gca,['../output/Poin_E_theta' num2str(PoinE(1))],'epsc');

% figure('name','region',...
%         'unit','normalized',...
%         'position',[0.0,0.0,0.3,0.45],... % figure position
%         'DefaultAxesFontSize',18,...
%         'DefaultAxesFontWeight','normal',...
%         'DefaultAxesLineWidth',1,...
%         'DefaultAxesTickLength',[0.02,0.015]);
% %for m=25
% for m=1:size(Poini.psi,2)-num_ones
%     poin_size=5;
% %x = Poini.E(:,m);
% x = Poini.p_zeta(:,m);
% %y = (Poini.theta(:,m));
% y=mod(Poini.theta(:,m)+pi, 2*pi) - pi;
% x(find(x==0))=nan;
% y(find(y==0))=nan;
% validIdx = ~isnan(x) & ~isnan(y);
% x= x(validIdx);
% y= y(validIdx);
% 
% scatter(x, y,poin_size, 'filled');hold on
% %scatter(x, y,poin_size, 'filled');hold on
% 
% end
% axis tight
% xlim([PoinPzeta_down,PoinPzeta_up]);
% %xlim([-0.1,0.1]);
% 
% %xlim([26,35]);
% box on
% xlabel('P_\zeta/\psi_w','fontsize',25)
% ylabel('\theta','fontsize',25)
% saveas(gca,'../output/Poin_Pzeta_theta','epsc');

