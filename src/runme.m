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
%% control parameters
ps_option = 'MG'; % PE, PLam, MG, GB, PLam_traj,Poincare
passing_option = 'co-passing'; % co-passing, counter-passing.
ant_option = 'off';% on, off.

%% 
%1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test particle orbit 
%with/without antenna
clc;clear;clear global; close all;
tokamak = 'DIIID';
ps_option = 'PLam_traj'; % PLam_traj
passing_option = 'co-passing';
ant_option = 'off';% on, off.!!!!!!!!!!!!!!!!!!!!!change this!!!!!!!!!!!!!!!!!

run run_paper.m;
%% 
%2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resonance line calculation 
%without antenna
clc;clear;clear global; close all;
tokamak = 'DIIID';
ps_option = 'PLam'; %PLam, MG, GB!!!!!!!!!!!!!!!!!!!!!change this!!!!!!!!!!!!!!!!!
ant_option = 'off';% on, off.
passing_option = 'co-passing';% co-passing, counter-passing.!!!!!!!!!!!!!!!!!!!!!change this!!!!!!!!!!!!!!!!!

run run_paper.m;
filename = [ps_option,'_',passing_option, '.mat'];
save_path = fullfile('../save/', filename);
save(save_path,'char_freq');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resonance line calculation and comparison with the perturbation distribution of GTC
%without antenna
clc;clear;clear global; close all;
tokamak = 'DIIID';
ps_option = 'GB'; %PLam, MG, GB!!!!!!!!!!!!!!!!!!!!!change this!!!!!!!!!!!!!!!!!
ant_option = 'off';% on, off.


passing_option = 'co-passing';
run run_paper.m;
filename = [ps_option,'_',passing_option, '.mat'];
save_path = fullfile('../save/', filename);
save(save_path,'char_freq');

passing_option = 'counter-passing'; 
run run_paper.m;
filename = [ps_option,'_',passing_option, '.mat'];
save_path = fullfile('../save/', filename);
save(save_path,'char_freq');

run plot_MAS_GTC.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%% 
%4.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kinetic Poincare plot
tokamak = 'DIIID';
ps_option = 'Poincare'; %PLam, MG, GB,Poincare
passing_option = 'co-passing';%don't matter 
ant_option = 'on';% on, off.

Eperp_in = 20; % unit is keV, mu*Ba (Ba is magnetic axis value)

PoinE=[26 34.5 30];%1.[26 34.5 30];%2.[46 53 45.09];%3.[54 61 51.93];%4.[64 70 60.22];%%[E_dowm E_up E_prime=E-omega/n*Pzeta]unit is keV!!!!!!!!!!!!!!!change this!!!!!!!!!!!
theta_rand = 1;   % 1 → include a small random phase perturbation in the initial particle distribution to
                  %     produce smoother and more visually pleasing Poincaré plots;;
                  % 0 → uniform initialization without any stochastic component.
run run_paper.m;
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%