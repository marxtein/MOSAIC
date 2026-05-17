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
function A = cal_char_freq(np,np_buffer,mstep,tstep,qpart,apart,eq_curr,eq_gradB,E_input,pzeta_input,lambda_input,orb_plot_opt)
global eq unit passing_option ps_option

%% load
%1 zpart(1,m): \psi
%2 zpart(2,m): \theta
%3 zpart(3,m): \zeta
%4 zpart(4,m): \rho_{||}
%5 zpart(5,m): \zeta without 2*pi periodicity
%6 zpart(6,m): sqrt(\mu)
%7 zpart(7,m): energy
%8 zpart(8,m): pzeta
%9 zpart(9,m): pitch
%10 zpart(10,m): \theta without 2*pi periodicity
%11 zpart(11,m): q

cmratio = qpart/apart;
c_inv = 1/qpart;

A = struct;
A.E = E_input(np_buffer+1:end-np_buffer);
A.lambda = lambda_input(np_buffer+1:end-np_buffer);
A.Pzeta_norm = pzeta_input(np_buffer+1:end-np_buffer);
A.exist = cell(1,np); %'yes','no': whether the solution of particle coordinate exists
A.type = cell(1,np); %'co-passing','counter-passing','trapped','potato','stagnation','inexistence'
A.confine = cell(1,np);
A.confine(:) = {'confined'}; %'confined','lost'
init = struct;
%figure;
for m = 1:np
    % Coordinate transformation accuracy. 5000 for ion; 50000 for electron



    if qpart > 0
        mpsi = 5000;
    elseif qpart < 0
        mpsi = 50000;
    else
        error('''qpart'' should be a number either positive or negative.');
    end
    psi_dum = linspace(0,eq.psiw,mpsi);
    rho_para2_A = zeros(1,mpsi);
    rho_para2_B = zeros(1,mpsi);
    E_gtc_unit = A.E(m)*1000/unit.energy_norm;
    Pzeta_gtc_unit = A.Pzeta_norm(m)*eq.psiw;
    lambda = A.lambda(m);
    
    psi_intersection = [];
    ind_intersection = [];
    theta_intersection = [];
    
    th_birth = 0; % only should be 0 or pi.
    n = 1; % only should be 1.
    while isempty(psi_intersection)
        if n > 2
            break;
        end
        
        for i = 1:mpsi
            b0_tmp =  qdspline.spline2d(0, psi_dum(i), th_birth, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            gpsi_tmp = qdspline.spline1d(0, psi_dum(i), eq.lsp, eq.dpsi, eq.gpsi);
            rho_para2_A(i) = (Pzeta_gtc_unit + psi_dum(i))^2/gpsi_tmp^2;
            rho_para2_B(i) = 2*(apart/qpart^2)/b0_tmp^2*(1 - lambda*b0_tmp)*E_gtc_unit;
        end
        
        for i = 2:mpsi
            Delta0 = rho_para2_A(i-1) - rho_para2_B(i-1);
            Delta1 = rho_para2_A(i) - rho_para2_B(i);
            if Delta0*Delta1 < 0
                %{
                if abs(Delta0) > abs(Delta1)
                    psi_intersection  = [psi_intersection , psi_dum(i)];
                    ind_intersection = [ind_intersection, i];
                else
                    psi_intersection  = [psi_intersection , psi_dum(i-1)];
                    ind_intersection = [ind_intersection, i-1];
                end
                %}
                delp_dum = psi_dum(i) - psi_dum(i-1);
                x_dum = abs(Delta0)/(abs(Delta0) + abs(Delta1))*delp_dum;
                psi_intersection  = [psi_intersection , psi_dum(i-1)+x_dum];
                
                theta_intersection = [theta_intersection, th_birth];
            end
        end
        n = n + 1;
        th_birth = th_birth + pi;
    end
    clear th_birth;
    
    if isempty(psi_intersection)
        A.exist{m} = 'no';
        A.confine{m} = 'non-defined';
        A.type{m} = 'inexistence';
        continue;
    else
        A.exist{m} = 'yes';
    end
    
    if qpart > 0
        % note that the maximal length of *_intersection array is two, so choose 1 or 2.
        if strcmp(passing_option,'co-passing')
            ind_tmp = min(2,length(psi_intersection));
        elseif strcmp(passing_option,'counter-passing')
            ind_tmp = min(1,length(psi_intersection));
        else
            error('Wrong passing_option.')
        end
    elseif qpart < 0 % For negative charge electron, the locations of co- and counter-passing particles exchange
        % note that the maximal length of *_intersection array is two, so choose 1 or 2.
        if strcmp(passing_option,'co-passing')
            ind_tmp = min(1,length(psi_intersection));
        elseif strcmp(passing_option,'counter-passing')
            ind_tmp = min(2,length(psi_intersection));
        else
            error('Wrong passing_option.')
        end
    else
        error('''qpart'' should be a number either positive or negative.')
    end
    %init.psi_ind(m) = ind_intersection(ind_tmp);
    init.psi(m) = psi_intersection(ind_tmp);
    init.theta(m) = theta_intersection(ind_tmp);
    clear ind_tmp;
    
    init.zeta(m) = 0;
    init.gpsi(m) = qdspline.spline1d(0, init.psi(m), eq.lsp, eq.dpsi, eq.gpsi);
    init.rho_para(m) = (Pzeta_gtc_unit + init.psi(m))/init.gpsi(m); % guarantee correct sign
    init.mu(m) = lambda*E_gtc_unit;
    
    %plot(psi_dum,rho_para2_A,'r');hold on;
    %plot(psi_dum,rho_para2_B,'b');hold on;
    %plot(init.psi(m),rho_para2_A(init.psi_ind(m)),'o','linewidth',1,'markersize',10,'color','k');hold on;
    %xlabel('psi (GTC\_unit)');
    %ylabel('\rho_{||} (GTC\_unit)');
    %%
    zpart(1,m) = init.psi(m);
    zpart(2,m) = init.theta(m);
    zpart(3,m) = init.zeta(m);
    zpart(4,m) = init.rho_para(m);
    zpart(6,m) = sqrt(init.mu(m));
    
    zpart(5,m) = zpart(3,m);
    zpart(10,m) = zpart(2,m);
end

% push time loop
var_num = 12;
part_pos=zeros(mstep,var_num,np);

%h=waitbar(0,'please wait');

for m = 1:np
    %t0=tic;
    %str=['Processing...',num2str(((m-1)*mstep)/(np*mstep)*100),'%'];
    %%waitbar(((m-1)*mstep)/(np*mstep),h,str);
    if strcmp(A.exist{m},'no')
        continue;
    end
    
    for istep = 1:mstep
        
        if zpart(1,m) > eq.psiw
            A.confine{m} = 'lost';
            break;
        end
        
        for irk = 1:2
            
            if irk == 1
                dtime = 0.5*tstep;
                zpart0(1:6,m) = zpart(1:6,m);
                zpart0(10,m) = zpart(10,m);
            else
                dtime = tstep;
            end
            
            pdum = zpart(1,m);
            tdum = zpart(2,m);
            zdum = zpart(3,m);
            
            mu = zpart(6,m)*zpart(6,m);
            
            g = qdspline.spline1d(0, pdum, eq.lsp, eq.dpsi, eq.gpsi);
            ri = qdspline.spline1d(0, pdum, eq.lsp, eq.dpsi, eq.ipsi);
            gp = qdspline.spline1d(1, pdum, eq.lsp, eq.dpsi, eq.gpsi);
            rip = qdspline.spline1d(1, pdum, eq.lsp, eq.dpsi, eq.ipsi);
            q = qdspline.spline1d(0, pdum, eq.lsp, eq.dpsi, eq.qpsi);
            deni = 1/(g*q + ri);
            
            b = qdspline.spline2d(0, pdum, tdum, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            dbdp = qdspline.spline2d(1, pdum, tdum, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            dbdt = qdspline.spline2d(2, pdum, tdum, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            upara = zpart(4,m)*b*cmratio;
            dedb = zpart(4,m)*zpart(4,m)*b*cmratio + c_inv*zpart(6,m)*zpart(6,m);
            
            energy = mu*b + 0.5*apart*upara^2;
            pzeta = g*zpart(4,m)-pdum;
            pitch = mu/energy;
            
            
            rdot = ((gp*zpart(4,m)*eq_curr-1)*dedb*dbdt)*deni;
            
            pdot = eq_gradB*(-dedb*g*dbdt)*deni;
            
            tdot = (upara*b*(1-zpart(4,m)*gp*eq_curr) + g*dedb*dbdp*eq_gradB)*deni;
            
            zdot = (upara*b*(q+zpart(4,m)*rip*eq_curr) - ri*dedb*dbdp*eq_gradB)*deni;
            
            zpart(1,m) = zpart0(1,m) + pdot*dtime;
            zpart(2,m) = zpart0(2,m) + tdot*dtime;
            zpart(3,m) = zpart0(3,m) + zdot*dtime;
            zpart(4,m) = zpart0(4,m) + rdot*dtime;
            zpart(5,m) = zpart0(5,m) + zdot*dtime;
            zpart(10,m) = zpart0(10,m) + tdot*dtime;
            
            zpart(2,m) = mod(zpart(2,m),2*pi);
            zpart(3,m) = mod(zpart(3,m),2*pi);
            
            zpart(7,m) = energy*unit.energy_norm/1000;
            zpart(8,m) = pzeta/eq.psiw;
            zpart(9,m) = pitch;
            zpart(11,m) = q;
            if irk == 2
                Rtmp = qdspline.spline2d(0, zpart(1,m), zpart(2,m), eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.xsp);
                Ztmp = qdspline.spline2d(0, zpart(1,m), zpart(2,m), eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.zsp);
                part_pos(istep,:,m) = [zpart(1,m) zpart(2,m) zpart(3,m) zpart(4,m) zpart(5,m) zpart(7,m) zpart(8,m) zpart(9,m) zpart(10,m) zpart(11,m) Rtmp Ztmp];
            end
            
        end
    end
    
    
   % disp(['No.' num2str(m),' particle push: ',num2str(toc(t0)),' (second)']);
end
%delete(h);
if strcmp(orb_plot_opt,'on')
    figure;
    for m = 1:np
        if strcmp(A.exist{m},'no')
            continue;
        end
        plot(part_pos(:,end-1,m),part_pos(:,end,m),'r.');hold on;
        [~,ind] = min(abs(init.psi(m) - eq.psi));
        plot(squeeze(eq.xsp(1,ind,:)),squeeze(eq.zsp(1,ind,:)),'k--','linewidth',2);hold on;
    end
    plot(squeeze(eq.xsp(1,end,:)),squeeze(eq.zsp(1,end,:)),'k-','linewidth',2);hold on;
    plot(eq.xsp(1,1,1),eq.zsp(1,1,1),'*','linewidth',1,'markersize',10);hold on;
    daspect([1 1 1]);
    grid on;
end

A.omega_b = zeros(1,np);
A.q_avrg = zeros(1,np);
A.omega_phi = zeros(1,np);
A.omega_d = zeros(1,np);

for m = 1:np
    if strcmp(A.exist{m},'no') | strcmp(A.confine{m},'lost')
        continue;
    end
    pol_cir_num =  abs((part_pos(end,9,m)-part_pos(1,9,m))/(2*pi)); % The number of poloidal orbit
    
    if all(sign(part_pos(:,4,m))==sign(part_pos(1,4,m))) & pol_cir_num > 1 % co-passing, counter-passing
        
        sign_theta = sign(part_pos(end,9,m) - part_pos(1,9,m));
        [~,step_one_taub] = min(abs(part_pos(:,9,m) - part_pos(1,9,m) - sign_theta*2*pi));
        A.omega_b(m) = (part_pos(step_one_taub,9,m)-part_pos(1,9,m))/(step_one_taub*tstep/unit.omegacp); % Poloidal freq, the unit is (rad/s)
        A.q_avrg(m) = mean(part_pos(1:step_one_taub,10,m));
        dzeta_one_taub = (part_pos(step_one_taub,5,m) - part_pos(1,5,m)); % Toroidal angle increment after one poloidal time period
        tau_b = abs(2*pi/A.omega_b(m)); % Time period of one poloidal orbit, the unit is (second)
        A.omega_phi(m) = dzeta_one_taub/tau_b; % toroidal frequency
        A.omega_d(m) = (dzeta_one_taub - (tau_b*A.omega_b(m))*A.q_avrg(m))/tau_b; % precession frequency
        if A.q_avrg(m) < 0 % BT direction is opposite to Ip
            if sign(part_pos(1,4,m)/qpart) > 0 % positive BT corresponds to negative Ip
                A.type{m} = 'counter-passing';
            elseif sign(part_pos(1,4,m)/qpart) < 0 % negative BT corresponds to positive Ip
                A.type{m} = 'co-passing';
            end
        elseif A.q_avrg(m) > 0
            if sign(part_pos(1,4,m)/qpart) > 0 % positive BT corresponds to positive Ip
                A.type{m} = 'co-passing';
            elseif sign(part_pos(1,4,m)/qpart) < 0 % negative BT corresponds to negative Ip
                A.type{m} = 'counter-passing';
            end
        end
        
    else % trapped, potato, stagnation
        [pks,locs] = findpeaks(part_pos(:,5,m));
        if isempty(locs)
            [~,locs] = findpeaks(part_pos(:,9,m));
            pks = part_pos(locs,5,m);
        end
        A.omega_d(m) = (pks(end)-pks(1))/((locs(end)-locs(1))*tstep/unit.omegacp);
        A.omega_phi(m) = A.omega_d(m);
        tau_b = (locs(end) - locs(1))/(length(locs) - 1)*tstep/unit.omegacp;
        A.q_avrg(m) = mean(part_pos(locs(1):locs(end),10,m));
        A.omega_b(m) = 2*pi/tau_b;
        if all(sign(part_pos(:,4,m))==sign(part_pos(1,4,m)))
            A.type{m} = 'stagnation';
        else
            if pol_cir_num < 1
                A.type{m} = 'trapped';
            else
                A.type{m} = 'potato';
            end
        end
    end
    
end

%% Only for the usage of testing particle orbit
if strcmp(ps_option,'PLam_traj')
    A.part_pos = part_pos;
    A.init = init;
    A.E_motion = A.part_pos(:,6,:);
    A.Pzeta_norm = A.part_pos(:,7,:);
    A.pitch = A.part_pos(:,8,:);
    A.R = A.part_pos(:,end-1,:);
    A.Z = A.part_pos(:,end,:);
end

end
