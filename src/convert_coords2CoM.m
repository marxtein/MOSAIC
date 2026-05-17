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
%% calculate Pz, lambda from (R, Z, E and pitch=v||/v)
clc;clear;clear global; close all;
global particle eq spdata plasma ps_option passing_option

tokamak = 'DIIID'; %%
datapath = '../input/';
spfile = 'spdata_rsae.dat'; %spdata_rsae.dat: (g,q<0);
particle.charge = 1.0; % particle charge normalized by elementary charge
particle.mass = 2.0; % particle mass normalized by proton mass
%----
E0 = 20; % keV, lower bound of energy range
E1 = 100;

E_grid = 581; % energy slice (fine grid for interpolation)
Pzeta_grid = 100; % Pzeta grid number

%% read in the equilibrium data
spdata = read_spdata(datapath,spfile,'mapping');
eq = eqdata_interface;
% clear spfile pfile datapath;
plasma.R0 = spdata.rmaj;  %(cm)
plasma.b0 = spdata.baxis*1e4;  %(gauss)    magnetic axis value

%% physics unit
run physics_unit;

%% create phase space boundary
run test_ps_cls;

%% check the poloidal flux at the wall compare with IDS data
% Please note that IDS use poloidal flux of SI units,
% representing the poloidal flux in SI units.
% In HAGIS and White-Chance, psi_p (ψ_p) is defined as
% the flux divided by 2π. Therefore, there is a difference
% of 2π compared to HAGIS's cpsurf.
cpsurf_IDS=eq.psiw*(plasma.b0/10000*(plasma.R0/100)^2)*2*pi;
disp(['psi_p_wall = ',num2str(cpsurf_IDS), '[Wb]']);


%% first calculate psi from R,Z
Rgrid=squeeze(spdata.xsp(1,:,:));
Zgrid=squeeze(spdata.zsp(1,:,:));

psigrid=spdata.psi;
thetagrid=linspace(0,2*pi,spdata.lst );
[psi2d,theta2d]=meshgrid(psigrid,thetagrid);
psi2d=psi2d';theta2d=theta2d';
%%
figure
plot(Rgrid,Zgrid)
hold on
plot(Rgrid(:,1),Zgrid(:,1),'r')
plot(Rgrid(:,end),Zgrid(:,end),'r*')
%get Fpsi, Ftheta of RZ
Fpsi   = scatteredInterpolant(Rgrid(:), Zgrid(:), psi2d(:), 'linear', 'nearest');
Ftheta = scatteredInterpolant(Rgrid(:), Zgrid(:), theta2d(:), 'linear', 'nearest');
%% to test
R=1.8;%[m]
Z=0.0;%[m]
pitch=1.0;
En=40000;%[eV]
en_gtc_unit = En/unit.energy_norm;
psi_val   = Fpsi(R,Z); % code unit, 0 to eq.psiw
theta_val = Ftheta(R,Z);
g = qdspline.spline1d(0, psi_val, eq.lsp, eq.dpsi, eq.gpsi);
b = qdspline.spline2d(0, psi_val, theta_val, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
cmratio = particle.charge/particle.mass;
upara = pitch * sqrt(2*en_gtc_unit/particle.mass);
rho_para = upara/b/cmratio;
mu= (en_gtc_unit -  0.5*particle.mass*upara^2)/b;
pzeta = g*rho_para-psi_val;
%------
lambda = mu/en_gtc_unit;
pzeta_norm = pzeta/eq.psiw;
disp(['lambda=',num2str(lambda),'; pzeta_norm=',num2str(pzeta_norm)])
%%