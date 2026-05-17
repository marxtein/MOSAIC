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
%Contains:
classdef ps_cls
    methods(Static)
        %%
        function A = surf_Pzeta_E_lambda(qpart,apart,E0_keV,E1_keV,E_grid,Pzeta_grid,by)
            global eq unit
            
            A = struct;
            A.charge = qpart; % unit is elementary charge (e)
            A.mass = apart; % unit is proton mass (mp)
            A.E0 = E0_keV;
            A.E1 = E1_keV;
            A.E_grid = E_grid;
            A.Pzeta_grid = Pzeta_grid;
            
            energy_1d_gtc_norm = linspace(A.E0,A.E1,A.E_grid)*1000/unit.energy_norm; % convert to GTC unit
            
            
            bmax = qdspline.spline2d(0, eq.psiw, pi, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            bmin = qdspline.spline2d(0, eq.psiw, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            baxis = eq.bsp(1,1,1);
            gpsi_wall = abs(eq.gpsi(1,end));
            gpsi_axis = abs(eq.gpsi(1,1));
            
            A.boundary = 7; % number of boundary line
            A.Pzeta = zeros(A.boundary,A.E_grid,A.Pzeta_grid);
            A.lambda = zeros(A.boundary,A.E_grid,A.Pzeta_grid);
            
            b_high_1d = zeros(A.E_grid,A.Pzeta_grid);
            b_low_1d = zeros(A.E_grid,A.Pzeta_grid);
            for i = 1:A.E_grid % formulas in the loop are normalized using GTC unit
                %
                A.Pzeta(1,i,:) = linspace(- eq.psiw - sqrt(2*energy_1d_gtc_norm(i)*apart/qpart^2)*gpsi_wall/bmin,...
                    - eq.psiw + sqrt(2*energy_1d_gtc_norm(i)*apart/qpart^2)*gpsi_wall/bmin,A.Pzeta_grid);
                
                A.lambda(1,i,:) = 1/bmin - 1/(2*energy_1d_gtc_norm(i))*(qpart^2/apart)*(A.Pzeta(1,i,:) + eq.psiw).^2*bmin/gpsi_wall^2;
                %
                A.Pzeta(2,i,:) = linspace(- eq.psiw - sqrt(2*energy_1d_gtc_norm(i)*apart/qpart^2)*gpsi_wall/bmax,...
                    - eq.psiw + sqrt(2*energy_1d_gtc_norm(i)*apart/qpart^2)*gpsi_wall/bmax,A.Pzeta_grid);
                
                A.lambda(2,i,:) = 1/bmax - 1/(2*energy_1d_gtc_norm(i))*(qpart^2/apart)*(A.Pzeta(2,i,:) + eq.psiw).^2*bmax/gpsi_wall^2;
                %
                A.Pzeta(3,i,:) = linspace(- sqrt(2*energy_1d_gtc_norm(i)*apart/qpart^2)*gpsi_axis/baxis,...
                    sqrt(2*energy_1d_gtc_norm(i)*apart/qpart^2)*gpsi_axis/baxis,A.Pzeta_grid);
                
                A.lambda(3,i,:) = 1/baxis - 1/(2*energy_1d_gtc_norm(i))*(qpart^2/apart)*A.Pzeta(3,i,:).^2*baxis/gpsi_axis^2;
                %
                A.Pzeta(4,i,:) = linspace(- eq.psiw,0,A.Pzeta_grid);
                
                for j = 1:A.Pzeta_grid
                    b_high_1d(i,j) =  qdspline.spline2d(0, - A.Pzeta(4,i,j), pi, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                end
                A.lambda(4,i,:) = 1./b_high_1d(i,:);
                %
                A.Pzeta(5,i,:) = linspace(- eq.psiw,0,A.Pzeta_grid);
                
                for j = 1:A.Pzeta_grid
                    b_low_1d(i,j) =  qdspline.spline2d(0, - A.Pzeta(5,i,j), 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                end
                A.lambda(5,i,:) = 1./b_low_1d(i,:);
                %
                A.Pzeta(6,i,:) = linspace(- eq.psiw,- eq.psiw,A.Pzeta_grid);
                
                A.lambda(6,i,:) = linspace(max(A.lambda(2,i,:)),max(A.lambda(1,i,:)),A.Pzeta_grid);

            end
            % transform the units of A.energy and A.Pzeta back to (keV, psiw) being consistent with inputs
            A.energy = repmat(energy_1d_gtc_norm,[A.boundary 1 A.Pzeta_grid])*(unit.energy_norm/1000);% unit is keV
            %size(A.energy)
            
            A.Pzeta = A.Pzeta/eq.psiw; % unit is psiw
            %by=1;
            if by==1
            for i = 1:A.E_grid

                if eq.gpsi(1,1)>0
                    sign_gpsi = 'positive';
                else
                    sign_gpsi = 'negative';
                end
                %A.energy(7,i,2)
                if A.energy(7,i,1)==0
                    A.Pzeta(7,i,:)=0;
                    A.lambda(6,i,:) =0;
                else
%size(A.energy(7,i,1))
                    PLam_2D_new = ps_cls.cro_sec_Pzeta_lambda(A.energy(7,i,1),A,sign_gpsi);

                    %A.Pzeta(7,i,:)=PLam_2D_new.Pzeta_bound_right;
                    A.Pzeta(7,i,:) = interp1( linspace(0,1,length(PLam_2D_new.Pzeta_bound_right)), ...
                          PLam_2D_new.Pzeta_bound_right, ...
                          linspace(0,1,size(A.Pzeta,3)) );
                    %A.lambda(7,i,:) =PLam_2D_new.lam_uni_grid;
                    A.lambda(7,i,:) = interp1( linspace(0,1,length(PLam_2D_new.lam_uni_grid)), ...
                          PLam_2D_new.lam_uni_grid, ...
                          linspace(0,1,size(A.Pzeta,3)) );
                    % i
                end

                %aaaaa=85468

            end
            end
            
        end
        %%
        function A = cro_sec_Pzeta_E(lam_in,E_grid_num,ps_3D)
            
            A = struct;
            A.lambda0 = lam_in;
            A.energy = zeros(ps_3D.boundary,ps_3D.E_grid);
            A.Pzeta = zeros(ps_3D.boundary,ps_3D.E_grid);
            A.inter_sec = cell(1,ps_3D.boundary);
            A.inter_sec(:) = {'Yes'}; %'Yes'or'No', use 'Yes' for initialization.
            
            for j = 1:ps_3D.boundary
                %For each E_grid, the same pzeta index doesn't correspond to the same Pzeta value between different boundary lines.
                %Thus, use E_grid for 2D loop in the following.
                if A.lambda0 > max(ps_3D.lambda(j,1,:)) | A.lambda0 < min(ps_3D.lambda(j,1,:))
                    warning(['No.',num2str(j),' boundary line doesnt intersect with lambda =',num2str(lam_in)]);
                    A.inter_sec{j} = 'No';
                    continue;
                end
                for i = 1:ps_3D.E_grid
                    [~,ind]=min(abs(ps_3D.lambda(j,i,:) - A.lambda0));
                    A.energy(j,i) = ps_3D.energy(j,i,ind);
                    A.Pzeta(j,i) = ps_3D.Pzeta(j,i,ind);
                end
            end
            
            A.E0 = ps_3D.E0;
            A.E1 = ps_3D.E1;
            A.E_uni_grid = linspace(A.E0,A.E1,E_grid_num);
            
            A.Pzeta_bound_right = zeros(ps_3D.boundary,E_grid_num);
            A.Pzeta_bound_left = zeros(ps_3D.boundary,E_grid_num);
            Pzeta_pass_axis = [-1 -1 0]; %Don't change! (Pzeta axis position of passing boundaries)
            
            for j = 1:3
                if strcmp(A.inter_sec{j},'Yes')
                    ind_right = find(A.Pzeta(j,:) > Pzeta_pass_axis(j));
                    Pzeta_right = A.Pzeta(j,ind_right);
                    E_right = A.energy(j,ind_right);
                    A.Pzeta_bound_right(j,:) = interp1(E_right, Pzeta_right, A.E_uni_grid, 'spline');
                    
                    ind_left = find(A.Pzeta(j,:) < Pzeta_pass_axis(j));
                    Pzeta_left = A.Pzeta(j,ind_left);
                    E_left = A.energy(j,ind_left);
                    A.Pzeta_bound_left(j,:) = interp1(E_left, Pzeta_left, A.E_uni_grid, 'spline');
                end
            end
            for j = 4:6
                if strcmp(A.inter_sec{j},'Yes')
                    A.Pzeta_bound_right(j,:) = interp1(A.energy(j,:), A.Pzeta(j,:), A.E_uni_grid, 'spline');
                    A.Pzeta_bound_left(j,:) = A.Pzeta_bound_right(j,:);
                end
            end
            
            %{
            ind_tmp1 = find(A.Pzeta(1,:) > -1);
            Pzeta_tmp1 = A.Pzeta(1,ind_tmp1);
            E_tmp1 = A.energy(1,ind_tmp1);
            %plot(Pzeta_tmp1,E_tmp1,'linewidth',3);hold on;
            
            ind_tmp2 = find(A.Pzeta(3,:) > 0);
            Pzeta_tmp2 = A.Pzeta(3,ind_tmp2);
            E_tmp2 = A.energy(3,ind_tmp2);
            %plot(Pzeta_tmp2,E_tmp2,'linewidth',3);hold on;
            
            A.E0 = ps_3D.E0;
            A.E1 = ps_3D.E1;
            A.E_uni_grid = linspace(A.E0,A.E1,E_grid_num);
            A.Pzeta_bound0 = interp1(E_tmp1, Pzeta_tmp1, A.E_uni_grid, 'spline');
            A.Pzeta_bound1 = interp1(E_tmp2, Pzeta_tmp2, A.E_uni_grid, 'spline');
            %}
        end
        %%
        function A = cro_sec_Pzeta_lambda(energy_in,ps_3D,sign_gpsi)
            A = struct;
            A.energy_in = energy_in; %keV
            
            A.Pzeta = zeros(ps_3D.boundary,ps_3D.Pzeta_grid);
            A.lambda = zeros(ps_3D.boundary,ps_3D.Pzeta_grid);
            energy_1d = squeeze(ps_3D.energy(1,:,1));
            [~,E_index] = min(abs(energy_1d - A.energy_in));
            for i = 1:ps_3D.boundary
                A.Pzeta(i,:) = squeeze(ps_3D.Pzeta(i,E_index,:));
                A.lambda(i,:) = squeeze(ps_3D.lambda(i,E_index,:));
            end
            
            A.energy_out = energy_1d(E_index);
            
            
            lam_grid_num = 200;
            Atmp = ps_cls.find_Pzeta_right(ps_3D.charge,ps_3D.mass,A.energy_out,min(A.lambda(1,:)),0.99*max(A.lambda(1,:)),lam_grid_num,sign_gpsi);
            A.lam_uni_grid = Atmp.lam;
            A.Pzeta_bound_left1 = interp1(A.lambda(1,find(A.Pzeta(1,:) > -1)), A.Pzeta(1,find(A.Pzeta(1,:) > -1)), A.lam_uni_grid, 'spline');
            A.Pzeta_bound_left2 = interp1(A.lambda(1,find(A.Pzeta(1,:) < -1)), A.Pzeta(1,find(A.Pzeta(1,:) < -1)), A.lam_uni_grid, 'spline');
            A.Pzeta_bound_right = max(Atmp.Pzeta_right,A.Pzeta_bound_left1);
            
            Pzeta_bound_right_co_EE_raw = [A.Pzeta(3,find(A.Pzeta(3,:) < 0)) flip(A.Pzeta(5,:))];
            lambda_bound_right_co_EE_raw = [A.lambda(3,find(A.Pzeta(3,:) < 0)) flip(A.lambda(5,:))];
            %A.Pzeta_bound_right_co_EE = interp1(lambda_bound_right_co_EE_raw, Pzeta_bound_right_co_EE_raw, A.lam_uni_grid, 'spline');
            A.Pzeta_bound_right_co_EE = interp1(lambda_bound_right_co_EE_raw, Pzeta_bound_right_co_EE_raw, A.lam_uni_grid, 'pchip');
            %A.Pzeta(7,:) = linspace(min(A.Pzeta_bound_left1),max(A.Pzeta_bound_left1),ps_3D.Pzeta_grid);
            %A.lambda(7,:) = interp1(A.Pzeta_bound_left1, A.lam_uni_grid, A.Pzeta(7,:), 'spline');
            %A.Pzeta(8,:) = linspace(min(A.Pzeta_bound_right),max(A.Pzeta_bound_right),ps_3D.Pzeta_grid);
            %A.lambda(8,:) = interp1(A.Pzeta_bound_right, A.lam_uni_grid, A.Pzeta(8,:), 'spline');
            
            Atmp2 = ps_cls.new_find_counter_passing_boundary(A);
            A.lam_uni_grid_counter_passing = Atmp2.lam_uni_grid_zone1;
            A.Pzeta_bound_counter_passing_left = Atmp2.Pzeta_bound0_zone1;
            A.Pzeta_bound_counter_passing_right = Atmp2.Pzeta_bound1_zone1;
        end
        %%
        function A = new_find_counter_passing_boundary(PLam_2D) % new algorithm to determine counter-passing boundaries, more robust in high energy regime.
            % TCP for trapped-counter passing; MACP for magnetic axis co-passing
            Pzeta_TCP_raw = [PLam_2D.Pzeta(2,find(PLam_2D.Pzeta(2,:) < -1)) PLam_2D.Pzeta(4,:)];
            lambda_TCP_raw = [PLam_2D.lambda(2,find(PLam_2D.Pzeta(2,:) < -1)) PLam_2D.lambda(4,:)];
            lambda_uni_grid = linspace(min(lambda_TCP_raw),max(lambda_TCP_raw),length(lambda_TCP_raw));
            %Pzeta_TCP_uni_grid = interp1(lambda_TCP_raw, Pzeta_TCP_raw, lambda_uni_grid, 'spline');
            Pzeta_TCP_uni_grid = interp1(lambda_TCP_raw, Pzeta_TCP_raw, lambda_uni_grid, 'pchip');
            
            Pzeta_MACP_raw = PLam_2D.Pzeta(3,find(PLam_2D.Pzeta(3,:) < 0));
            lambda_MACP_raw = PLam_2D.lambda(3,find(PLam_2D.Pzeta(3,:) < 0));
            %Pzeta_MACP_uni_grid = interp1(lambda_MACP_raw, Pzeta_MACP_raw, lambda_uni_grid, 'spline');
            Pzeta_MACP_uni_grid = interp1(lambda_MACP_raw, Pzeta_MACP_raw, lambda_uni_grid, 'pchip');
            
            lambda_intersection = [];
            ind_intersection = [];
            for i = 2:length(lambda_uni_grid)
                Delta0 = Pzeta_TCP_uni_grid(i-1) - Pzeta_MACP_uni_grid(i-1);
                Delta1 = Pzeta_TCP_uni_grid(i) - Pzeta_MACP_uni_grid(i);
                if Delta0*Delta1 < 0
                    % always choose lower lambda grid to guarantee Pzeta_left<Pzeta_right
                    lambda_intersection  = [lambda_intersection , lambda_uni_grid(i-1)];
                    ind_intersection = [ind_intersection, i-1];
                end
            end
            
            if isempty(lambda_intersection)
                warning('Cannot find the boundaries for confined counter-passing particles.');
                A.lam_uni_grid_zone1 = lambda_uni_grid;
                A.Pzeta_bound0_zone1 = Pzeta_TCP_uni_grid;
                A.Pzeta_bound1_zone1 = Pzeta_MACP_uni_grid;
            else
                lam0 = min(lambda_uni_grid);
                lam1 = lambda_intersection(end);
                A.lam_uni_grid_zone1 = linspace(lam0,lam1,length(PLam_2D.lam_uni_grid));
                %A.Pzeta_bound0_zone1 = interp1(lambda_uni_grid(1:ind_intersection(end)), Pzeta_TCP_uni_grid(1:ind_intersection(end)), A.lam_uni_grid_zone1, 'spline');
                %A.Pzeta_bound1_zone1 = interp1(lambda_uni_grid(1:ind_intersection(end)), Pzeta_MACP_uni_grid(1:ind_intersection(end)), A.lam_uni_grid_zone1, 'spline');
                A.Pzeta_bound0_zone1 = interp1(lambda_uni_grid(1:ind_intersection(end)), Pzeta_TCP_uni_grid(1:ind_intersection(end)), A.lam_uni_grid_zone1, 'cubic');
                A.Pzeta_bound1_zone1 = interp1(lambda_uni_grid(1:ind_intersection(end)), Pzeta_MACP_uni_grid(1:ind_intersection(end)), A.lam_uni_grid_zone1, 'cubic');
            end
        end
        %%
        function A = find_counter_passing_boundary(PLam_2D) % old algorithm from phase_space2.m
            
            ind_counp = find(PLam_2D.Pzeta(2,:) < -1);
            ind_cp_neg = find(PLam_2D.Pzeta(3,:) < 0);
            Pzeta_cp_neg = PLam_2D.Pzeta(3,ind_cp_neg);
            lambda_cp_neg = PLam_2D.lambda(3,ind_cp_neg);
            
            Pzeta_trap_low = Pzeta_cp_neg;
            lambda_trap_low = interp1(PLam_2D.Pzeta(4,:),PLam_2D.lambda(4,:),Pzeta_trap_low, 'spline');
            [~,ind_inter]=min(abs(lambda_trap_low - lambda_cp_neg));
            lambda_inter = lambda_trap_low(ind_inter);
            Pzeta_inter = Pzeta_trap_low(ind_inter);
            
            Pzeta_bound0_A = PLam_2D.Pzeta(2,ind_counp);
            lam_bound0_A = PLam_2D.lambda(2,ind_counp);
            Pzeta_bound0_B = PLam_2D.Pzeta(4,find(PLam_2D.lambda(4,:)<lambda_inter));
            lam_bound0_B = PLam_2D.lambda(4,find(PLam_2D.lambda(4,:)<lambda_inter));
            lam_bound0_tmp = [lam_bound0_A lam_bound0_B];
            Pzeta_bound0_tmp = [Pzeta_bound0_A Pzeta_bound0_B];
            
            A.lam_uni_grid_zone1 = linspace(min(lam_bound0_tmp),max(lam_bound0_tmp),length(PLam_2D.lam_uni_grid));
            A.Pzeta_bound0_zone1 = interp1(lam_bound0_tmp,Pzeta_bound0_tmp,A.lam_uni_grid_zone1, 'spline');
            A.Pzeta_bound1_zone1 = interp1(lambda_cp_neg, Pzeta_cp_neg, A.lam_uni_grid_zone1, 'spline');
            
            ind = A.lam_uni_grid_zone1(find(A.Pzeta_bound0_zone1>A.Pzeta_bound1_zone1));
            if ~isempty(ind)
                A.lam_uni_grid_zone1(ind) = [];
                A.Pzeta_bound0_zone1(ind) = [];
                A.Pzeta_bound1_zone1(ind) = [];
            end
        end
        %%
        function A = find_Pzeta_right(qpart,apart,energy_in,lam0,lam1,lam_grid_num,sign_gpsi)
            global eq unit
            
            energy_gtc_norm = energy_in*1000/unit.energy_norm; % convert to GTC normalization
            
            A = struct;
            A.lam = linspace(lam0,lam1,lam_grid_num);
            
            A.psi_ind = zeros(1,lam_grid_num);
            A.psi = zeros(1,lam_grid_num);
            A.gpsi = zeros(1,lam_grid_num);
            A.b0 = zeros(1,lam_grid_num);
            A.rhopara = zeros(1,lam_grid_num);
            A.Pzeta_right = zeros(1,lam_grid_num);
            
            
            mpsi = 1000;
            psi_dum = linspace(0,eq.psiw,mpsi);
            term1 = zeros(1,mpsi);
            term2 = zeros(1,mpsi);
            
            for j = 1:lam_grid_num
                lam_dum = A.lam(j);
                for i = 1:mpsi
                    b0_tmp =  qdspline.spline2d(0, psi_dum(i), 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                    db0_tmp =  qdspline.spline2d(1, psi_dum(i), 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                    gpsi_tmp = qdspline.spline1d(0, psi_dum(i), eq.lsp, eq.dpsi, eq.gpsi);
                    dgpsi_tmp = qdspline.spline1d(1, psi_dum(i), eq.lsp, eq.dpsi, eq.gpsi);
                    
                    parallel_energy_term = max(1-lam_dum*b0_tmp,0);
                    if strcmp(sign_gpsi,'negative')
                        rhopara = - sqrt(2*energy_gtc_norm*(parallel_energy_term)/b0_tmp^2*apart/qpart^2);
                    elseif strcmp(sign_gpsi,'positive')
                        rhopara = sqrt(2*energy_gtc_norm*(parallel_energy_term)/b0_tmp^2*apart/qpart^2);
                    else
                        error('Wrong ''sign_gpsi'' string. Set ''negative'' or ''positive''.')
                    end
                    
                    term1(i) = - 1 + rhopara*dgpsi_tmp;
                    term2(i) = - gpsi_tmp/rhopara*(apart/qpart^2/b0_tmp^2)*(-2*energy_gtc_norm/b0_tmp + lam_dum*energy_gtc_norm)*db0_tmp;
                end
                
                %plot(psi_dum,term1,'r');hold on;
                %plot(psi_dum,term2,'b')
                %xlim([0.01*psi_dum(end),psi_dum(end)]);
                
                [~,A.psi_ind(j)]=min(abs(term1 - term2));
                A.psi(j) = psi_dum(A.psi_ind(j));
                A.gpsi(j) = qdspline.spline1d(0, A.psi(j), eq.lsp, eq.dpsi, eq.gpsi);
                A.b0(j) =  qdspline.spline2d(0, A.psi(j), 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                if strcmp(sign_gpsi,'negative')
                    A.rhopara(j) = -sqrt(2*energy_gtc_norm*(1-lam_dum*A.b0(j))/A.b0(j)^2*apart/qpart^2);
                elseif strcmp(sign_gpsi,'positive')
                    A.rhopara(j) = sqrt(2*energy_gtc_norm*(1-lam_dum*A.b0(j))/A.b0(j)^2*apart/qpart^2);
                else
                    error('Wrong ''sign_gpsi'' string. Set ''negative'' or ''positive''.')
                end
                A.Pzeta_right(j) = (A.gpsi(j)*A.rhopara(j) - A.psi(j))/eq.psiw;
                
            end
        end
    end
end
