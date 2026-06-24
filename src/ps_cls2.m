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
classdef ps_cls2
    methods(Static)
        %%
        function A =surf_Pzeta_E_mu(qpart,apart,E0_keV,E1_keV,E_grid,Pzeta0,Pzeta1,Pzeta_grid,E_max,E_grid_num,psi_diag_norm)
            global eq unit
            A = struct;
            A.charge = qpart; % unit is elementary charge (e)
            A.mass = apart; % unit is proton mass (mp)
            A.E0 = E0_keV;
            A.E1 = E1_keV;
            A.E_grid = E_grid;
            A.Pzeta_grid = Pzeta_grid;

            energy_1d_gtc_norm = linspace(A.E0,A.E1,A.E_grid)*1000/unit.energy_norm; % convert to GTC unit
            A.boundary = 7; % number of boundary line
            A.Pzeta = zeros(A.boundary,A.E_grid,A.Pzeta_grid);
            A.mu = zeros(A.boundary,A.E_grid,A.Pzeta_grid);
            A.energy = repmat(energy_1d_gtc_norm,[A.boundary 1 A.Pzeta_grid])*(unit.energy_norm/1000);% unit is keV
            size(A.energy)

            A.Pzeta = A.Pzeta/eq.psiw; % unit is psiw
            for i = 1:A.E_grid

                if eq.gpsi(1,1)>0
                    sign_gpsi = 'positive';
                else
                    sign_gpsi = 'negative';
                end
                %A.energy(7,i,2)
                if A.energy(7,i,1)==0
                    A.Pzeta(7,i,:)=0;
                    A.energy(6,i,:) =0;
                else

                    MG = ps_cls2.Pzeta_E_fixed_mu(A.charge,A.mass,A.energy(7,i,1),Pzeta0,Pzeta1,Pzeta_grid,E_max,Pzeta_grid,psi_diag_norm);

                    A.Pzeta(1,i,:)=MG.pzeta_norm;
                    A.Pzeta(2,i,:)=MG.pzeta_norm;
                    A.Pzeta(3,i,:)=MG.pzeta_norm;
                    A.Pzeta(4,i,:)=MG.pzeta_norm;
                    A.Pzeta(5,i,:)=MG.pzeta_norm;
                    A.Pzeta(6,i,:)=MG.pzeta_norm;
                    A.Pzeta(7,i,:)=MG.Pzeta_bound5;
                    A.mu(1:6,i,:)=MG.E_pzeta(1:6,:);
                    A.mu(7,i,:)=MG.E_uni_grid1;

                end

                %aaaaa=85468

            end
        end

        %%
        function A = Pzeta_E_fixed_mu(qpart,apart,Eperp_in,Pzeta0,Pzeta1,Pzeta_grid_num,E_max,E_grid_num,psi_diag_norm)
            global unit eq
            A= struct;

            A.E_perp_axis = Eperp_in; % keV
            mu_norm = A.E_perp_axis*1000/unit.energy_norm;
            A.pzeta_norm0 = Pzeta0;
            A.pzeta_norm1 = Pzeta1;
            A.pzeta_num = Pzeta_grid_num;
            A.pzeta_norm = linspace(A.pzeta_norm0,A.pzeta_norm1,A.pzeta_num);

            bmax = qdspline.spline2d(0, eq.psiw, pi, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            bmin = qdspline.spline2d(0, eq.psiw, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            baxis = eq.bsp(1,1,1);
            gpsi_wall = abs(eq.gpsi(1,end));
            gpsi_axis = abs(eq.gpsi(1,1));

            A.E_pzeta = zeros(6,A.pzeta_num);
            A.E_pzeta(1,:) = (qpart^2/apart)*(A.pzeta_norm*eq.psiw + eq.psiw).^2/(2*gpsi_wall^2)*bmin^2 + mu_norm*bmin;
            A.E_pzeta(2,:) = (qpart^2/apart)*(A.pzeta_norm*eq.psiw + eq.psiw).^2/(2*gpsi_wall^2)*bmax^2 + mu_norm*bmax;
            A.E_pzeta(3,:) = (qpart^2/apart)*(A.pzeta_norm*eq.psiw).^2/(2*gpsi_axis^2)*baxis^2 + mu_norm*baxis;

            for i = 1:A.pzeta_num
                crit_pzeta = (A.pzeta_norm(i) + 1)*A.pzeta_norm(i);
                if crit_pzeta <= 0
                    b_high_tmp = qdspline.spline2d(0, - A.pzeta_norm(i)*eq.psiw, pi, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                    A.E_pzeta(4,i) = mu_norm*b_high_tmp;
                    b_low_tmp = qdspline.spline2d(0, - A.pzeta_norm(i)*eq.psiw, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                    A.E_pzeta(5,i) = mu_norm*b_low_tmp;
                else
                    A.E_pzeta(4,i) = nan;
                    A.E_pzeta(5,i) = nan;
                end
            end

            b_diag = qdspline.spline2d(0, psi_diag_norm*eq.psiw, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            g_diag = qdspline.spline1d(0, psi_diag_norm*eq.psiw, eq.lsp, eq.dpsi, eq.gpsi);
            A.E_pzeta(6,:) = (qpart^2/apart)*(A.pzeta_norm*eq.psiw + psi_diag_norm*eq.psiw).^2/(2*g_diag^2)*b_diag^2 + mu_norm*b_diag;

            E_convert = unit.energy_norm/1000; %keV
            A.E_pzeta = A.E_pzeta*E_convert; % convert the normalized A.E_pzeta value to SI unit keV

            ind_tmp1 = find(A.pzeta_norm > -1);
            Pzeta_tmp1 = A.pzeta_norm(ind_tmp1);
            E_tmp1 = A.E_pzeta(1,ind_tmp1);

            ind_tmp2 = find(A.pzeta_norm > 0);
            Pzeta_tmp2 = A.pzeta_norm(ind_tmp2);
            E_tmp2 = A.E_pzeta(3,ind_tmp2);

            ind_tmp3 = find(A.pzeta_norm < -1);
            ind_tmp31 = find(A.pzeta_norm > -1 & A.pzeta_norm < -0.1);
            ind_tmp32=[ind_tmp3 ind_tmp31];
            Pzeta_tmp3 = A.pzeta_norm(ind_tmp3);
            Pzeta_tmp31 = A.pzeta_norm(ind_tmp31);
            Pzeta_tmp32 = [Pzeta_tmp3 Pzeta_tmp31];
            E_tmp3 = A.E_pzeta(2,ind_tmp3);
            E_tmp31 = A.E_pzeta(4,ind_tmp31);
            E_tmp32=[E_tmp3 E_tmp31];

            ind_tmp4 = find(A.pzeta_norm < 0);
            Pzeta_tmp4 = A.pzeta_norm(ind_tmp4);
            E_tmp4 = A.E_pzeta(3,ind_tmp4);

            %A.E_uni_grid = linspace(ceil(min(A.E_pzeta(3,:))),E_max,E_grid_num);
            A.E_uni_grid = linspace(ceil(min(A.E_pzeta(1,:))),E_max,E_grid_num);
            A.Pzeta_bound0 = interp1(E_tmp1, Pzeta_tmp1, A.E_uni_grid, 'spline');
            A.Pzeta_bound1 = interp1(E_tmp2, Pzeta_tmp2, A.E_uni_grid, 'spline');
            A.Pzeta_bound2 = interp1(E_tmp32, Pzeta_tmp32, A.E_uni_grid, 'spline');
            A.Pzeta_bound3 = interp1(E_tmp4, Pzeta_tmp4, A.E_uni_grid, 'spline');
            diff_y = A.Pzeta_bound2 -A.Pzeta_bound3;                        % 两条曲线的差

            sign_change = diff_y(1:end-1) .* diff_y(2:end) < 0;  % 判断符号是否变化
            cross_idx = find(sign_change);            % 交点区间的下标
            %size(cross_idx)
            if isempty(cross_idx)==0
            A.Pzeta_bound2(1:cross_idx(end))=nan;
            A.Pzeta_bound3(1:cross_idx(end))=nan;
            end
            if eq.gpsi(1,1)>0
                sign_gpsi = 'positive';
            else
                sign_gpsi = 'negative';
            end

            Atmp = ps_cls2.find_Pzeta_right_fixed_mu(qpart,apart,Eperp_in,0.99*Eperp_in/min(A.E_pzeta(1,:)),0.99*Eperp_in/max(A.E_uni_grid),E_grid_num,sign_gpsi);
            A.E_uni_grid1 = Eperp_in./Atmp.lam;
            A.Pzeta_bound5=Atmp.Pzeta_right;
            A.Pzeta_bound4 = interp1(A.E_uni_grid1, Atmp.Pzeta_right, A.E_uni_grid, 'spline');
        end

        %%
        function A =surf_Pzeta_lambda_mu(qpart,apart,E0_keV,E1_keV,E_grid,Pzeta0,Pzeta1,Pzeta_grid,lam_min,lam_grid_num)
            global eq unit
            A = struct;
            A.charge = qpart; % unit is elementary charge (e)
            A.mass = apart; % unit is proton mass (mp)
            A.E0 = E0_keV;
            A.E1 = E1_keV;
            A.E_grid = E_grid;
            A.Pzeta_grid = Pzeta_grid;

            energy_1d_gtc_norm = linspace(A.E0,A.E1,A.E_grid)*1000/unit.energy_norm; % convert to GTC unit
            A.boundary = 7; % number of boundary line
            A.Pzeta = zeros(A.boundary,A.E_grid,A.Pzeta_grid);
            A.lambda = zeros(A.boundary,A.E_grid,A.Pzeta_grid);
            A.energy = repmat(energy_1d_gtc_norm,[A.boundary 1 A.Pzeta_grid])*(unit.energy_norm/1000);% unit is keV
            size(A.energy)

            A.Pzeta = A.Pzeta/eq.psiw; % unit is psiw
            for i = 1:A.E_grid

                if eq.gpsi(1,1)>0
                    sign_gpsi = 'positive';
                else
                    sign_gpsi = 'negative';
                end
                %A.energy(7,i,2)
                if A.energy(7,i,1)==0
                    A.Pzeta(7,i,:)=0;
                    A.energy(6,i,:) =0;
                else

                    %GB = ps_cls2.Pzeta_lambda_fixed_mu(particle.charge,particle.mass,Eperp_in,Pzeta0,Pzeta1,Pzeta_grid_num,lam_min,lam_grid_num,sign_gpsi);


                    GB = ps_cls2.Pzeta_lambda_fixed_mu(A.charge,A.mass,A.energy(7,i,1),Pzeta0,Pzeta1,Pzeta_grid,lam_min,lam_grid_num,sign_gpsi);

                    A.Pzeta(1,i,:)=GB.pzeta_norm;
                    A.Pzeta(2,i,:)=GB.pzeta_norm;
                    A.Pzeta(3,i,:)=GB.pzeta_norm;
                    A.Pzeta(4,i,:)=GB.pzeta_norm;
                    A.Pzeta(5,i,:)=GB.pzeta_norm;
                    A.Pzeta(6,i,:)=GB.Pzeta_bound0;
                    A.Pzeta(7,i,:)=GB.Pzeta_bound1;
                    A.lambda(1:5,i,:)=GB.lambda_pzeta(1:5,:);
                    A.lambda(6,i,:)=GB.lambda_uni_grid;
                    A.lambda(7,i,:)=GB.lambda_uni_grid;

                end

                %aaaaa=85468

            end
        end

        %%
        function A = Pzeta_lambda_fixed_mu(qpart,apart,Eperp_in,Pzeta0,Pzeta1,Pzeta_grid_num,lambda_min,lambda_grid_num,sign_gpsi)
            global unit eq

            A= struct;
            A.E_perp_axis = Eperp_in; % keV
            mu_norm = A.E_perp_axis*1000/unit.energy_norm;
            A.pzeta_norm0 = Pzeta0;
            A.pzeta_norm1 = Pzeta1;
            A.pzeta_num = Pzeta_grid_num;
            A.pzeta_norm = linspace(A.pzeta_norm0,A.pzeta_norm1,A.pzeta_num);

            bmax = qdspline.spline2d(0, eq.psiw, pi, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            bmin = qdspline.spline2d(0, eq.psiw, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            baxis = eq.bsp(1,1,1);
            gpsi_wall = abs(eq.gpsi(1,end));
            gpsi_axis = abs(eq.gpsi(1,1));

            A.lambda_pzeta = zeros(5,A.pzeta_num);

            tmp = 0.5*(qpart^2/apart)*(A.pzeta_norm*eq.psiw + eq.psiw).^2*bmin^2/(gpsi_wall^2) + mu_norm*bmin;
            A.lambda_pzeta(1,:) = mu_norm./tmp;
            tmp = 0.5*(qpart^2/apart)*(A.pzeta_norm*eq.psiw + eq.psiw).^2*bmax^2/(gpsi_wall^2) + mu_norm*bmax;
            A.lambda_pzeta(2,:) = mu_norm./tmp;
            tmp = 0.5*(qpart^2/apart)*(A.pzeta_norm*eq.psiw).^2*baxis^2/(gpsi_axis^2) + mu_norm*baxis;
            A.lambda_pzeta(3,:) = mu_norm./tmp;


            for i = 1:A.pzeta_num
                crit_pzeta = (A.pzeta_norm(i) + 1)*A.pzeta_norm(i);
                if crit_pzeta <= 0
                    b_high_tmp = qdspline.spline2d(0, - A.pzeta_norm(i)*eq.psiw, pi, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                    A.lambda_pzeta(4,i) = 1./b_high_tmp;
                    b_low_tmp = qdspline.spline2d(0, - A.pzeta_norm(i)*eq.psiw, 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                    A.lambda_pzeta(5,i) = 1./b_low_tmp;
                else
                    A.lambda_pzeta(4,i) = nan;
                    A.lambda_pzeta(5,i) = nan;
                end
            end

            %{
            Pzeta_tmp0 = A.pzeta_norm(find(A.pzeta_norm > -1));
            lambda_tmp0 = A.lambda_pzeta(1,find(A.pzeta_norm > -1));
            a = A.lambda_pzeta(5,find(~isnan(A.lambda_pzeta(5,:))));
            ax = A.pzeta_norm(find(~isnan(A.lambda_pzeta(5,:))));
            b = A.lambda_pzeta(3,find(A.pzeta_norm > 0));
            bx = A.pzeta_norm(find(A.pzeta_norm > 0));
            lambda_tmp1 = [ax bx];
            Pzeta_tmp1 = [a b];
            A.lambda_uni_grid = linspace(lambda_min,max(A.lambda_pzeta(1,:)),lambda_grid_num);
            A.Pzeta_bound0 = interp1(lambda_tmp0, Pzeta_tmp0, A.lambda_uni_grid, 'spline');
            A.Pzeta_bound1 = interp1(lambda_tmp1, Pzeta_tmp1, A.lambda_uni_grid, 'spline');
            %}

            lambda_min = max(lambda_min,A.lambda_pzeta(3,end));
            Pzeta_tmp0 = A.pzeta_norm(find(A.pzeta_norm > -1));
            lambda_tmp0 = A.lambda_pzeta(1,find(A.pzeta_norm > -1));
            Atmp = ps_cls2.find_Pzeta_right_fixed_mu(qpart,apart,Eperp_in,lambda_min,0.99*max(A.lambda_pzeta(1,:)),lambda_grid_num,sign_gpsi);

            A.lambda_uni_grid = Atmp.lam;
            A.Pzeta_bound0 = interp1(lambda_tmp0, Pzeta_tmp0, A.lambda_uni_grid, 'spline');
            A.Pzeta_bound1 = Atmp.Pzeta_right;

            %% counter passing boundary
            pzeta_left = [A.pzeta_norm(find(A.pzeta_norm < -1)) A.pzeta_norm(find(~isnan(A.lambda_pzeta(4,:))))];
            lambda_left = [A.lambda_pzeta(2,find(A.pzeta_norm < -1)) A.lambda_pzeta(4,find(~isnan(A.lambda_pzeta(4,:))))];
            %plot(pzeta_left,lambda_left,'o');hold on;

            pzeta_right = A.pzeta_norm(find(A.pzeta_norm < 0));
            lambda_right = A.lambda_pzeta(3,find(A.pzeta_norm < 0));
            %plot(pzeta_right,lambda_right,'o');hold on;

            lam1_count_p = min(lambda_left(end),lambda_right(end));
            lam0_count_p = max(lambda_left(1), lambda_right(1));
            lambda_uni_grid_tmp = linspace(lam0_count_p,lam1_count_p,lambda_grid_num);

            pzeta_bound0_tmp = interp1(lambda_left, pzeta_left, lambda_uni_grid_tmp, 'spline');%cubic
            pzeta_bound1_tmp = interp1(lambda_right, pzeta_right, lambda_uni_grid_tmp, 'spline');%cubic

            lambda_uni_grid_tmp_cut = lambda_uni_grid_tmp(find(pzeta_bound0_tmp<pzeta_bound1_tmp));
            A.lambda_uni_grid_counter_passing = linspace(lambda_uni_grid_tmp_cut(1),lambda_uni_grid_tmp_cut(end),lambda_grid_num);

            A.Pzeta_bound0_counter_passing = interp1(lambda_uni_grid_tmp, pzeta_bound0_tmp, A.lambda_uni_grid_counter_passing, 'spline');%cubic
            A.Pzeta_bound1_counter_passing = interp1(lambda_uni_grid_tmp, pzeta_bound1_tmp, A.lambda_uni_grid_counter_passing, 'spline');%cubic

            %plot(A.Pzeta_bound0_counter_passing,A.lambda_uni_grid_counter_passing,'*');hold on;
            %plot(A.Pzeta_bound1_counter_passing,A.lambda_uni_grid_counter_passing,'*');hold on;
      

        end

        %%
        function A = find_Pzeta_right_fixed_mu(qpart,apart,Eperp_in,lam0,lam1,lam_grid_num,sign_gpsi)
            global eq unit

            A = struct;

            A.E_perp_axis = Eperp_in; % keV
            mu_norm = A.E_perp_axis*1000/unit.energy_norm;
            A.lam = linspace(lam0,lam1,lam_grid_num);
            energy_gtc_norm = mu_norm./A.lam; % convert to GTC normalization

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
                        rhopara = - sqrt(2*energy_gtc_norm(j)*(parallel_energy_term)/b0_tmp^2*apart/qpart^2);
                    elseif strcmp(sign_gpsi,'positive')
                        rhopara = sqrt(2*energy_gtc_norm(j)*(parallel_energy_term)/b0_tmp^2*apart/qpart^2);
                    else
                        error('Wrong ''sign_gpsi'' string. Set ''negative'' or ''positive''.')
                    end

                    term1(i) = - 1 + rhopara*dgpsi_tmp;
                    term2(i) = - gpsi_tmp/rhopara*(apart/qpart^2/b0_tmp^2)*(-2*energy_gtc_norm(j)/b0_tmp + lam_dum*energy_gtc_norm(j))*db0_tmp;
                end

                %plot(psi_dum,term1,'r');hold on;
                %plot(psi_dum,term2,'b')
                %xlim([0.01*psi_dum(end),psi_dum(end)]);

                [~,A.psi_ind(j)]=min(abs(term1 - term2));
                A.psi(j) = psi_dum(A.psi_ind(j));
                A.gpsi(j) = qdspline.spline1d(0, A.psi(j), eq.lsp, eq.dpsi, eq.gpsi);
                A.b0(j) =  qdspline.spline2d(0, A.psi(j), 0, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
                if strcmp(sign_gpsi,'negative')
                    A.rhopara(j) = -sqrt(2*energy_gtc_norm(j)*(1-lam_dum*A.b0(j))/A.b0(j)^2*apart/qpart^2);
                elseif strcmp(sign_gpsi,'positive')
                    A.rhopara(j) = sqrt(2*energy_gtc_norm(j)*(1-lam_dum*A.b0(j))/A.b0(j)^2*apart/qpart^2);
                else
                    error('Wrong ''sign_gpsi'' string. Set ''negative'' or ''positive''.')
                end
                A.Pzeta_right(j) = (A.gpsi(j)*A.rhopara(j) - A.psi(j))/eq.psiw;

            end

        end

    end
end
