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

global particle PE_2D char_freq MG PLam_2D GB ps_option passing_option

qpart = particle.charge;
apart = particle.mass;

np_buffer = max(1,floor(0.01*np)); %avoiding touch the orbit boundary

ind_low = 1; % start index for phase-space scan, dont change
if strcmp(ps_option,'PE')

    ind_high = length(PE_2D.E_uni_grid);
    E_max = PE_2D.E_uni_grid(end);

    if strcmp(passing_option,'co-passing')
        %plot(PE_2D.Pzeta_bound_right(1,:),PE_2D.E_uni_grid,'linewidth',3);hold on;
        %plot(PE_2D.Pzeta_bound_right(3,:),PE_2D.E_uni_grid,'linewidth',3);hold on;
    elseif strcmp(passing_option,'counter-passing')
        %plot(PE_2D.Pzeta_bound_left(2,:),PE_2D.E_uni_grid,'linewidth',3);hold on;
        %plot(PE_2D.Pzeta_bound_left(3,:),PE_2D.E_uni_grid,'linewidth',3);hold on;
    else
        error('Wrong passing_option.')
    end
    xlabel('Pzeta (\psi_w)');
    ylabel('E (keV)');

elseif strcmp(ps_option,'MG') | strcmp(ps_option,'Poincare')

    ind_high = length(MG.E_uni_grid);
    E_max = MG.E_uni_grid(end);

    % plot(MG.Pzeta_bound0,MG.E_uni_grid,'linewidth',3);hold on;
    % plot(MG.Pzeta_bound4,MG.E_uni_grid,'linewidth',3);hold on;
    %plot(MG.Pzeta_bound5,MG.E_uni_grid1,'linewidth',3);hold on;
    xlabel('Pzeta (\psi_w)');
    ylabel('E (keV)');

elseif strcmp(ps_option,'GB')

    if strcmp(passing_option,'co-passing')
        ind_high = length(GB.lambda_uni_grid);
        %plot(GB.Pzeta_bound0,GB.lambda_uni_grid,'linewidth',3);hold on;
        %plot(GB.Pzeta_bound1,GB.lambda_uni_grid,'linewidth',3);hold on;
    elseif strcmp(passing_option,'counter-passing')
        ind_high = length(GB.lambda_uni_grid_counter_passing);
        %plot(GB.Pzeta_bound0_counter_passing,GB.lambda_uni_grid_counter_passing,'linewidth',3);hold on;
        %plot(GB.Pzeta_bound1_counter_passing,GB.lambda_uni_grid_counter_passing,'linewidth',3);hold on;
    else
        error('Wrong passing_option.')
    end

    xlabel('Pzeta (\psi_w)');
    ylabel('\lambda');

elseif strcmp(ps_option,'PLam')

    if strcmp(passing_option,'co-passing')
        ind_high = length(PLam_2D.lam_uni_grid); % co-passing
    elseif strcmp(passing_option,'counter-passing')
        ind_high = length(PLam_2D.lam_uni_grid_counter_passing); % counter-passing
    else
        error('Wrong passing_option')
    end
    E_max = energy_in;

elseif strcmp(ps_option,'PLam_traj') | strcmp(ps_option,'PLam_traj2')

    ind_high = 1; % end index for phase space scan, dont change for this case.
    E_max = energy_in; %keV

end

%%
char_freq = struct;
char_freq.omega_phi = [];
char_freq.omega_d = [];
char_freq.omega_b = [];
char_freq.q_avrg = [];
char_freq.Pzeta_norm = [];
char_freq.E = [];
char_freq.lambda = [];
char_freq.exist = [];
char_freq.confine = [];
char_freq.type = [];
char_freq.psi_init = [];
char_freq.theta_init = [];
char_freq.rhopara_init = [];
if strcmp(ps_option,'PE')|strcmp(ps_option,'MG') |strcmp(ps_option,'GB')|strcmp(ps_option,'PLam')
    h=waitbar(0,'please wait');
end
for i = ind_low:ind_high
    t0=tic;
    if strcmp(ps_option,'PE')|strcmp(ps_option,'MG') |strcmp(ps_option,'GB')|strcmp(ps_option,'PLam')
    str=['Processing...',num2str((i - ind_low)/(ind_high - ind_low + 1)*100),'%'];
    waitbar((i - ind_low)/(ind_high - ind_low + 1),h,str);
    end
    if strcmp(ps_option,'PE')

        if strcmp(passing_option,'co-passing')
            pzeta_tmp = linspace(PE_2D.Pzeta_bound_right(1,i),PE_2D.Pzeta_bound_right(3,i),np+2*np_buffer);
        elseif strcmp(passing_option,'counter-passing')
            pzeta_tmp = linspace(PE_2D.Pzeta_bound_left(2,i),PE_2D.Pzeta_bound_left(3,i),np+2*np_buffer);
        else
            error('Wrong passing_option')
        end
        lambda_tmp = ones(1,np+2*np_buffer)*PE_2D.lambda0;
        E_tmp = ones(1,np+2*np_buffer)*PE_2D.E_uni_grid(i);
        E_ratio = PE_2D.E_uni_grid(i)/E_max;

    elseif strcmp(ps_option,'MG') | strcmp(ps_option,'Poincare')
        if strcmp(passing_option,'co-passing')
            pzeta_tmp = linspace(MG.Pzeta_bound0(i),MG.Pzeta_bound4(i),np+2*np_buffer);
        elseif strcmp(passing_option,'counter-passing')
            pzeta_tmp = linspace(MG.Pzeta_bound2(i),MG.Pzeta_bound3(i),np+2*np_buffer);
        else
            error('Wrong passing_option')
        end
        lambda_tmp = ones(1,np+2*np_buffer)*(MG.E_perp_axis/MG.E_uni_grid(i));

        E_tmp = ones(1,np+2*np_buffer)*MG.E_uni_grid(i);
        E_ratio = MG.E_uni_grid(i)/E_max;

    elseif strcmp(ps_option,'GB')

        if strcmp(passing_option,'co-passing')
            pzeta_tmp = linspace(GB.Pzeta_bound0(i),GB.Pzeta_bound1(i),np+2*np_buffer);
            lambda_tmp = ones(1,np+2*np_buffer)*GB.lambda_uni_grid(i);
            E_tmp = GB.E_perp_axis./lambda_tmp;
            E_ratio = GB.lambda_uni_grid(1)/GB.lambda_uni_grid(i);
        elseif strcmp(passing_option,'counter-passing')
            pzeta_tmp = linspace(GB.Pzeta_bound0_counter_passing(i),GB.Pzeta_bound1_counter_passing(i),np+2*np_buffer);
            lambda_tmp = ones(1,np+2*np_buffer)*GB.lambda_uni_grid_counter_passing(i);
            E_tmp = GB.E_perp_axis./lambda_tmp;
            E_ratio = GB.lambda_uni_grid_counter_passing(1)/GB.lambda_uni_grid_counter_passing(i);
        else
            error('Wrong passing_option')
        end

    elseif strcmp(ps_option,'PLam')

        if strcmp(passing_option,'co-passing')

            if qpart > 0
                pzeta_tmp = linspace(PLam_2D.Pzeta_bound_left1(i),PLam_2D.Pzeta_bound_right(i),np+2*np_buffer);
            elseif qpart < 0 % For negative charge electron, the locations of co- and counter-passing particles exchange
                pzeta_tmp = linspace(PLam_2D.Pzeta_bound_left1(i),PLam_2D.Pzeta_bound_right_co_EE(i),np+2*np_buffer);
            else
                error('''qpart'' should be a number either positive or negative.');
            end
            lambda_tmp = ones(1,np+2*np_buffer)*max(0,PLam_2D.lam_uni_grid(i));

        elseif strcmp(passing_option,'counter-passing')

            if qpart > 0
                pzeta_tmp = linspace(PLam_2D.Pzeta_bound_counter_passing_left(i),PLam_2D.Pzeta_bound_counter_passing_right(i),np+2*np_buffer);
                lambda_tmp = ones(1,np+2*np_buffer)*PLam_2D.lam_uni_grid_counter_passing(i);
            elseif qpart < 0 % For negative charge electron, the locations of co- and counter-passing particles exchange
                pzeta_tmp = linspace(PLam_2D.Pzeta_bound_left1(i),PLam_2D.Pzeta_bound_right(i),np+2*np_buffer);
                lambda_tmp = ones(1,np+2*np_buffer)*max(0,PLam_2D.lam_uni_grid(i));
            else
                error('''qpart'' should be a number either positive or negative.');
            end

        else
            error('Wrong passing_option')
        end
        E_tmp = ones(1,np+2*np_buffer)*E_max;
        E_ratio = 1; % affect actual time step in push, don't change

    elseif strcmp(ps_option,'PLam_traj')

        P_num = 5;
        l_num = floor(np/5);

        lam = linspace(traj.lam0,traj.lam1,l_num);
        Pzeta_norm = linspace(traj.Pzeta_norm0,traj.Pzeta_norm1,P_num);
        [Pzeta_2D,lam_2D] = meshgrid(Pzeta_norm,lam);

        pzeta_tmp = reshape(Pzeta_2D,l_num*P_num,1);
        lambda_tmp = reshape(lam_2D,l_num*P_num,1);
        E_tmp = E_max*ones(l_num*P_num,1);

        if length(pzeta_tmp) < np
            pzeta_tmp(end+1:np) = pzeta_tmp(end);
            lambda_tmp(end+1:np) = lambda_tmp(end);
            E_tmp(end+1:np) = E_tmp(end);
        end

        E_ratio = 1;
        np_buffer = 0;

    elseif strcmp(ps_option,'PLam_traj2')

        lambda_tmp = linspace(traj.lam0,traj.lam1,np)';
        pzeta_tmp = traj.Pzeta_norm0*ones(np,1);
        E_tmp = E_max*ones(np,1);

        E_ratio = 1;
        np_buffer = 0;

    end

    orb_plot_opt = 'off';

    eq_curr = 1; % 0 for without and 1 for with curvature drift.
    eq_gradB = 1; % 0 for without and 1 for with grad-B drift.
    if strcmp(ant_option,'on') && strcmp(ps_option,'Poincare')
        moduleName = 'Poincare_perturb';
        fprintf('Starting module: %s...\n', moduleName);
        tic;
        run poincare_perturb;
        elapsedTime = toc;
        fprintf('=> Module [%s] execution time: %.4f seconds\n\n', moduleName, elapsedTime);
        run plot_poin.m
        break
    else
        if strcmp(ant_option,'on')
            A = cal_char_freq_perturb(np,np_buffer,mstep,tstep/sqrt(E_ratio),qpart,apart,eq_curr,eq_gradB,E_tmp,pzeta_tmp,lambda_tmp,orb_plot_opt);
        else
            A = cal_char_freq(np,np_buffer,mstep,tstep/sqrt(E_ratio),qpart,apart,eq_curr,eq_gradB,E_tmp,pzeta_tmp,lambda_tmp,orb_plot_opt);
        end
        char_freq.omega_phi = [char_freq.omega_phi;A.omega_phi];
        char_freq.omega_d = [char_freq.omega_d;A.omega_d];
        char_freq.omega_b = [char_freq.omega_b;A.omega_b];
        char_freq.q_avrg = [char_freq.q_avrg;A.q_avrg];
        char_freq.Pzeta_norm = [char_freq.Pzeta_norm;A.Pzeta_norm];
        char_freq.E = [char_freq.E;A.E];
        char_freq.lambda = [char_freq.lambda;A.lambda];
        char_freq.exist =  [char_freq.exist;A.exist];
        char_freq.confine = [char_freq.confine;A.confine];
        char_freq.type = [char_freq.type;A.type];
        if isfield(A,'part_pos')
            char_freq.psi_init = [char_freq.psi_init;squeeze(A.part_pos(1,1,:))];
            char_freq.theta_init = [char_freq.theta_init;squeeze(A.part_pos(1,2,:))];
            char_freq.rhopara_init = [char_freq.rhopara_init;squeeze(A.part_pos(1,4,:))];
        end
        disp(['Push No.' num2str(i),' particle group(','np=',num2str(np),'): ',num2str(toc(t0)),' (second)']);
    end
end
if strcmp(ps_option,'PE')|strcmp(ps_option,'MG') |strcmp(ps_option,'GB')|strcmp(ps_option,'PLam')
delete(h);
end
