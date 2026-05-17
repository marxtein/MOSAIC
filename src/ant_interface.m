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
function  f = ant_interface(file,Amp_coeff,mmin,mmax)

f=load(file);

if mmin<f.m1
    mmin = f.m1;
end
if mmax>f.m2
    mmax = f.m2;
end
f.m_modes = linspace(mmin,mmax,mmax-mmin+1);
f.num_modes = length(f.m_modes);

disp(['ant.m_modes = ',num2str(f.m_modes)]);
disp(['ant.num_modes = ',num2str(f.num_modes)]);
disp(['ant.Amp_coeff = ',num2str(Amp_coeff)]);

%%
ndim = f.m2 - f.m1 + 1;
ilow = min(abs(mmin),abs(mmax))+2;%1;
ihigh = max(abs(mmin),abs(mmax))+2;%ndim;

close all;

    A1=[1,1,1.2,1]*1.1;
    
figure('name','raw',...
    'unit','normalized',...
    'position',[0.0,0.0,1,1]./A1,... % figure position
    'DefaultAxesFontSize',20,...
    'DefaultAxesFontWeight','normal',...
    'DefaultAxesLineWidth',1,...
    'DefaultAxesTickLength',[0.01,0.01]);
subplot(331)

for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,real(f.phim_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
    end
xlim([f.rho_tor(1)+0.1 f.rho_tor(end)-0.1]);

set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$Re(\delta\phi_m)/(B_aR_0V_{Ap}/c)$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;

subplot(332)
for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,real(f.aparam_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$Re(\delta A_{||m})/(B_aR_0)$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;

subplot(333)
for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,real(f.Brm_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$Re(\delta B_{rm}/B_a)$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;

subplot(334)
for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,imag(f.phim_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$Im(\delta\phi_m)/(B_aR_0V_{Ap}/c)$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;

subplot(335)
for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,imag(f.aparam_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$Im(\delta A_{||m})/(B_aR_0)$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;

subplot(336)
for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,imag(f.Brm_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$Im(\delta B_{rm}/B_a)$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;

subplot(337)
for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,abs(f.phim_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$|\delta\phi_m|/(B_aR_0V_{Ap}/c)$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;

subplot(338)
for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,abs(f.aparam_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$|\delta A_{||m}|/(B_aR_0)$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;

subplot(339)
for i=ihigh:-1:ilow
    h(i)=plot(f.rho_tor,abs(f.Brm_x(:,i)),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(i+f.m1-1))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',20);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast','fontsize',20)
xlabel('$\sqrt{\psi_T/\psi_{Tw}}$','Interpreter','latex','fontsize',25);
ylabel('$|\delta B_{rm}/B_a|$','Interpreter','latex','fontsize',25);
%title('raw data');
grid on;
 saveas(gca,'../output/ant','epsc');
%%
f.Ldp = operator2.gradient1d(3, f.psi);

f.phi_real=cell(1,f.num_modes);
f.phi_imag=cell(1,f.num_modes);
f.apara_real=cell(1,f.num_modes);
f.apara_imag=cell(1,f.num_modes);
f.dptdp_real=cell(1,f.num_modes);
f.dptdp_imag=cell(1,f.num_modes);
f.dapdp_real=cell(1,f.num_modes);
f.dapdp_imag=cell(1,f.num_modes);

for i = f.num_modes:-1:1
    f.dptdp_real{i} = Amp_coeff*(f.Ldp*real(f.phim_x(:,f.m_modes(i)-f.m1+1)))';
    f.dptdp_imag{i} = Amp_coeff*(f.Ldp*imag(f.phim_x(:,f.m_modes(i)-f.m1+1)))';
    f.dapdp_real{i} = Amp_coeff*(f.Ldp*real(f.aparam_x(:,f.m_modes(i)-f.m1+1)))';
    f.dapdp_imag{i} = Amp_coeff*(f.Ldp*imag(f.aparam_x(:,f.m_modes(i)-f.m1+1)))';
    
    f.phi_real{i} = Amp_coeff*real(f.phim_x(:,f.m_modes(i)-f.m1+1))';
    f.phi_imag{i} = Amp_coeff*imag(f.phim_x(:,f.m_modes(i)-f.m1+1))';
    f.apara_real{i} = Amp_coeff*real(f.aparam_x(:,f.m_modes(i)-f.m1+1))';
    f.apara_imag{i} = Amp_coeff*imag(f.aparam_x(:,f.m_modes(i)-f.m1+1))';
end

%%
A1=[1,1,1.2,1]*1.1;
figure('name','gradient',...
    'unit','normalized',...
    'position',[0.0,0.0,1,1]./A1,... % figure position
    'DefaultAxesFontSize',20,...
    'DefaultAxesFontWeight','normal',...
    'DefaultAxesLineWidth',1,...
    'DefaultAxesTickLength',[0.01,0.01]);

subplot(241)
for i=f.num_modes:-1:1
    h2(i)=plot(f.rho_tor,f.dptdp_real{i},'linewidth',2);hold on;
    lgd2(i)={['m=',num2str(abs(f.m_modes(i)))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',16);
hlgd=legend(h2,lgd2,'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast')
xlabel('$\sqrt{\psi_T}$','Interpreter','latex','fontsize',25);
ylabel('$Re(\partial \phi_m/\partial\psi)$','Interpreter','latex','fontsize',25);
title('antenna data');
grid on;

subplot(242)
for i=f.num_modes:-1:1
    h2(i)=plot(f.rho_tor,f.dapdp_real{i},'linewidth',2);hold on;
    lgd2(i)={['m=',num2str(abs(f.m_modes(i)))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',16);
hlgd=legend(h2,lgd2,'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast')
xlabel('$\sqrt{\psi_T}$','Interpreter','latex','fontsize',25);
ylabel('$Re(\partial\delta A_{||m}/\partial\psi)$','Interpreter','latex','fontsize',25);
title('antenna data');
grid on;

subplot(243)
for i=f.num_modes:-1:1
    h2(i)=plot(f.rho_tor,f.dptdp_imag{i},'linewidth',2);hold on;
    lgd2(i)={['m=',num2str(abs(f.m_modes(i)))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',16);
hlgd=legend(h2,lgd2,'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast')
xlabel('$\sqrt{\psi_T}$','Interpreter','latex','fontsize',25);
ylabel('$Im(\partial \phi_m/\partial\psi)$','Interpreter','latex','fontsize',25);
title('antenna data');
grid on;

subplot(244)
for i=f.num_modes:-1:1
    h2(i)=plot(f.rho_tor,f.dapdp_imag{i},'linewidth',2);hold on;
    lgd2(i)={['m=',num2str(abs(f.m_modes(i)))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',16);
hlgd=legend(h2,lgd2,'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast')
xlabel('$\sqrt{\psi_T}$','Interpreter','latex','fontsize',25);
ylabel('$Im(\partial\delta A_{||m}/\partial\psi)$','Interpreter','latex','fontsize',25);
title('antenna data');
grid on;

subplot(245)
for i=ihigh:-1:ilow
    h(i)=plot(0.5*(f.rho_tor(1:end-1) + f.rho_tor(2:end)),diff(real(f.phim_x(:,i))')./diff(f.psi),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(abs(i+f.m1-1)))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',16);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast')
xlabel('$\sqrt{\psi_T}$','Interpreter','latex','fontsize',25);
ylabel('$Re(\partial \phi_m/\partial\psi)$','Interpreter','latex','fontsize',25);
title('raw data');
grid on;

subplot(246)
for i=ihigh:-1:ilow
    h(i)=plot(0.5*(f.rho_tor(1:end-1) + f.rho_tor(2:end)),diff(real(f.aparam_x(:,i))')./diff(f.psi),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(abs(i+f.m1-1)))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',16);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast')
xlabel('$\sqrt{\psi_T}$','Interpreter','latex','fontsize',25);
ylabel('$Re(\partial\delta A_{||m}/\partial\psi)$','Interpreter','latex','fontsize',25);
title('raw data');
grid on;

subplot(247)
for i=ihigh:-1:ilow
    h(i)=plot(0.5*(f.rho_tor(1:end-1) + f.rho_tor(2:end)),diff(imag(f.phim_x(:,i))')./diff(f.psi),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(abs(i+f.m1-1)))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',16);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast')
xlabel('$\sqrt{\psi_T}$','Interpreter','latex','fontsize',25);
ylabel('$Im(\partial \phi_m/\partial\psi)$','Interpreter','latex','fontsize',25);
title('raw data');
grid on;

subplot(248)
for i=ihigh:-1:ilow
    h(i)=plot(0.5*(f.rho_tor(1:end-1) + f.rho_tor(2:end)),diff(imag(f.aparam_x(:,i))')./diff(f.psi),'linewidth',2);hold on;
    lgd(i)={['m=',num2str(abs(abs(i+f.m1-1)))]};
    hold on;
end
xlim([f.rho_tor(1) f.rho_tor(end)]);
set(gca,'fontsize',16);
hlgd=legend(h(ihigh:-1:ilow),lgd(ihigh:-1:ilow),'location','northwest'); legend('boxoff');
set(hlgd,'location','northeast')
xlabel('$\sqrt{\psi_T}$','Interpreter','latex','fontsize',25);
ylabel('$Im(\partial\delta A_{||m}/\partial\psi)$','Interpreter','latex','fontsize',25);
title('raw data');
grid on;

end