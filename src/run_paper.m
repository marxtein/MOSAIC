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
%clc;clear;clear global; close all;
global particle eq ant spdata plasma ps_option passing_option ant_option

%tokamak = 'DIIID';
%ps_option = 'MG'; % PE, PLam, MG, GB, PLam_traj,Poincare
%passing_option = 'co-passing'; % co-passing, counter-passing.
%ant_option = 'off';% on, off.
%% create Pzeta-energy-lambda phase space (used for PE and PLam/PLam_traj cases)
E0 = 2; % keV, lower bound of energy range
E1 = 60; %  keV, upper bound of energy range
E_grid = 581; % energy slice (fine grid for interpolation)
Pzeta_grid = 100; % Pzeta grid number

np = 40; % number of particle per row (column number is fixed as 100 (PE, MG and GB) or 200 (PLam), one can change in ps_cls.m)
mstep = 5000; % number of time step
tstep = 10; % time step normalized by omega_cp^-1 (omega_cp is proton cyclotron frequency)
n = 4;
omega0 = 82.03*1000*2*pi;
if ~exist('../output','dir');
    mkdir('../output');
end
if ~exist('../save','dir');
    mkdir('../save');
end
%% choose 2D phase space coordinates
if strcmp(ps_option,'PE') % Pzeta-energy at fixed lambda

    lam_in = 0.4; % mu*Ba/E, pitch angle
    E_grid_num = 100; % energy dimension grid

elseif strcmp(ps_option,'PLam') | strcmp(ps_option,'PLam_traj') | strcmp(ps_option,'PLam_traj2') % Pzeta-lambda at fixed energy

    energy_in = 60; % keV, total energy
    if energy_in > E1 | energy_in < E0
        error('The energy is out of phase space range, which should be between E0 and E1.')
    end

elseif strcmp(ps_option,'MG')  | strcmp(ps_option,'Poincare') % Pzeta-energy at fixed mu

    Eperp_in = 20; % unit is keV, mu*Ba (Ba is magnetic axis value)
    Pzeta0 = -1.7; % unit is eq.psiw
    Pzeta1 = 0.5; % unit is eq.psiw
    Pzeta_grid_num = 200;
    E_max = 90; % unit is keV
    E_grid_num = 100;
    psi_diag_norm = 0.215; % unit is eq.psiw
    if strcmp(ps_option,'Poincare')  % Pzeta-energy at fixed mu
        %Eperp_in = 20; % unit is keV, mu*Ba (Ba is magnetic axis value)

        np = 200; % number of particle
        mstep = 80000; % number of time step
        PoinE_down=PoinE(1);%keV
        PoinE_up=PoinE(2);%keV
        E_prime_input=PoinE(3);%keV,E_prime=E-omega/n*Pzeta
        PoinPzeta_down=-0.1; % unit is eq.psiw
        PoinPzeta_up=0.1; % unit is eq.psiw
        psi_in_norm=0;%unit is eq.psiw
        psi_out_norm=0.8;%unit is eq.psiw
        amp_mod=1/20;
        sign_pa=-1;
        omega_n=82030;%kHz
    end
elseif strcmp(ps_option,'GB') % Pzeta-lambda at fixed mu

    Eperp_in = 20; % unit is keV, mu*Ba (Ba is magnetic axis value)
    Pzeta0 = -1.7; % unit is eq.psiw
    Pzeta1 = 0.8; % unit is eq.psiw
    Pzeta_grid_num = 100;
    lam_min = 0.0;
    lam_grid_num = 100;

end

if strcmp(ps_option,'PLam_traj')
    % sample np number particles uniformly in [lam0 lam1 Pzeta_norm0 Pzeta_norm1]
    traj.lam0 = 0.2;
    traj.lam1 = 1.2;
    traj.Pzeta_norm0 = -1.5;
    traj.Pzeta_norm1 = 0;

    np = 25;
end

if strcmp(ps_option,'PLam_traj2')

    %traj.lam0 = 0.913;
    %traj.lam1 = 0.913;
    traj.lam0 = 0;
    traj.lam1 = 1;
    traj.Pzeta_norm0 = -0.292;
    traj.Pzeta_norm1 = -0.292;

    np = 1000;
end

%% set particle charge/mass and equilibrium path
if strcmp(tokamak,'DIIID')

    particle.charge = 1.0; % particle charge normalized by elementary charge
    particle.mass = 2.0; % particle mass normalized by proton mass
    datapath = '../input/';
    spfile = 'spdata_rsae.dat'; %spdata_rsae.dat: (g,q<0);

elseif strcmp(tokamak,'EAST2')

    particle.charge = -1.0; % particle charge normalized by elementary charge
    particle.mass = 1.0/1837.0; % particle mass normalized by proton mass
    datapath = '../input/';
    spfile = 'spdata_112786_05950.dat';

elseif strcmp(tokamak,'ITER')

    particle.charge = 2.0; % particle charge normalized by elementary charge
    particle.mass = 4.0; % particle mass normalized by proton mass
    datapath = '../input/';
    spfile = 'spdata_eqdsk16HR.dat';

elseif strcmp(tokamak,'Analytical equilibrium')
    plasma.R0 = 83.5  %(cm)
    plasma.b0 = 20125  %(gauss)    magnetic axis value
    particle.charge = 1.0; % particle charge normalized by elementary charge
    particle.mass = 1.0; % particle mass normalized by proton mass
    eq.qcoef=[0.82;1.1;1.0]; % q=q(1) + q(2)*psi_norm +q(3)*psi_norm^2
    eq.psiw = 0.0375;

end

%% read in the equilibrium data
if strcmp(tokamak,'Analytical equilibrium')
    run construct_analytic_eq.m
    if strcmp(ant_option,'on')
        ant = ant_interface('RSAE_ant.mat',1,-13,-11);
    end
else
    spdata = read_spdata(datapath,spfile,'mapping');
    eq = eqdata_interface;
    if strcmp(ant_option,'on')
        ant = ant_interface('RSAE_ant.mat',1,-13,-11);
        ant.omega=0.0027;%82.03 kHz
    end
    clear spfile pfile datapath;
    plasma.R0 = spdata.rmaj;  %(cm)
    plasma.b0 = spdata.baxis*1e4;  %(gauss)    magnetic axis value
end
%% physics unit
run physics_unit;

%% create phase space boundary
run test_ps_cls;

%% push particle and calculate frequency
moduleName = 'Particle Motion Module';
fprintf('Starting module: %s...\n', moduleName);
tic;
run test_push;
elapsedTime = toc;
fprintf('=> Module [%s] execution time: %.4f seconds\n\n', moduleName, elapsedTime);

%% output plot
if strcmp(ps_option,'PLam')
    run plot_PLam_fix_E;
    run plot_line;
elseif strcmp(ps_option,'PLam_traj')
    iplot=1;
    run plot_PLam_fix_E_traj;
    iplot=2;
    run plot_PLam_fix_E_traj;
    if strcmp(ant_option,'off')
    run mosaic_par_benchmark;
    end
elseif strcmp(ps_option,'PLam_traj2')
    save('single_E200_Pzeta0292_refactor.mat','char_freq','PLam_2D')
    run plot_multi_EE_resonance3;

elseif strcmp(ps_option,'PE')
    run plot_PE_fix_lambda;

elseif strcmp(ps_option,'MG') | strcmp(ps_option,'GB')
    %run plot_MG_GB;
    run plot_line;
end


