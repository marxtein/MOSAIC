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
%iplot=3;    %1:lambda-Pzeta fix E,2:E-Pzeta fix lambda
if strcmp(ps_option,'MG')
    iplot = 3;
elseif strcmp(ps_option,'GB')
    iplot = 4;
elseif strcmp(ps_option,'PLam')
    iplot = 1;
end
tstep=4;
ii=6;%[1,2,3,4,5,6,7]
%ii=6;
iii=3;

allname={'Particle Flux';'Energy Flux';'Angmomentum Flux';'delta E';'delta u_{para}';'delta f';'delta angmomentum'};
parkind={'ti';'te';'fi';'fe'};
Pzeta_in=0.5;



draw_i=3;%change 1:i;3:fi

emax_inv=0.25;egrids=101;%change
lam_inv=0.75;lambdagrids=100;%change
Pzeta_inv=0.65;Pzetagrids=100;%change
mub_inv=0.25;mubgrids=100;
% Set figure properties
meshi0=1.86317834;%1.2998%1.86317834;%gtc.out first meshi
meshf0=23.3529551;%23.786174%17.2%23.786174%23.3529551;%20.783799%23.3529551;%gtc.out first meshf


%%Main plotting
Pzeta=(Pzeta_in +1)*(((draw_i == 1) +(draw_i == 3))*Pzeta_inv*Pzetagrids);
Pzetagrid=ceil(Pzeta);

Pzetaaxis=(0:Pzetagrids-1)/(((draw_i == 1) +(draw_i == 3))*Pzeta_inv*Pzetagrids)-1;%/Pzeta_inv;
Eaxis=(0:egrids-1)/(((draw_i == 1) /meshi0+(draw_i == 3) /meshf0)*emax_inv*egrids);
lambdaaxis=(0:lambdagrids-1)/(((draw_i == 1) +(draw_i == 3))*lam_inv*lambdagrids);
clear Bplot
%figure;
switch iplot
    case 1
        tstep=2;
        clear x y z
        %E=E_in*((draw_i == 1) /1.86317834+(draw_i == 3) /23.2529351)*emax_inv*egrids;
         E_in=energy_in;
         E=E_in*((draw_i == 1) /meshi0+(draw_i == 3) /meshf0)*emax_inv*egrids;
        Egrid=ceil(E)+1;
        x=Pzetaaxis;
        y=lambdaaxis;
         Bplot=((draw_i == 1)*Data3d.data3di +(draw_i == 3)*Data3d.data3df);
        %Bplot=(draw_i == 1)*data3di
        z(:,:)=sum(Bplot(tstep,Egrid-2:Egrid+2,ii,:,:),2);

        thre = max(max(abs(z)))/60; % Set a threshold, e.g. 0.1, to determine which values should be ignored
        z(abs(z) < thre) = NaN; % Set values close to 0 to NaN
        [C, h] = contourf(x, y,z, 100);hold on;colormap('jet');set(h, 'linecolor', 'none');
        cb=colorbar('northoutside');
        cb.FontSize = 15;
        cb.Label.FontName = 'Times New Roman';
        cb.FontWeight = 'bold';
        cb.LineWidth = 1;
        % filename = [char(allname(ii)),'_', char(parkind(draw_i)),'E','_',num2str(energy_in),char(passingp_option),char(ps_option), '.jpg'];full_path = fullfile(path, filename);saveas(gcf, full_path);
        caxis([min(min(z)) max(max(z))])
    case 2

        clear x y z
        lambda_in=lam_in;
        lambda=lambda_in*((draw_i == 1) +(draw_i == 3))*lam_inv*lambdagrids;
        Lambdagrid=ceil(lambda);
        x=Pzetaaxis;
        y=Eaxis;
         Bplot=((draw_i == 1)*Data3d.data3di +(draw_i == 3)*Data3d.data3df);
        z(:,:)=Bplot(tstep,:,ii,Lambdagrid,:);

        thre = max(max(abs(z)))/100; % Set a threshold, e.g. 0.1, to determine which values should be ignored  
        z(abs(z) < thre) = NaN; % Set values close to 0 to NaN  
     
        [C, h] = contourf(x, y,z, 100);hold on;colormap('jet');set(h, 'linecolor', 'none');
        cb=colorbar('northoutside');
        cb.FontSize = 15;
        cb.Label.FontName = 'Times New Roman';
        cb.FontWeight = 'bold';
        cb.LineWidth = 1;

    case 3
        clear x y z
        mub_in=Eperp_in;
        mub=mub_in*((draw_i == 1) /meshi0+(draw_i == 3) /meshf0)*mub_inv*mubgrids;
        mubgrid=ceil(mub);
        x=Pzetaaxis;
        y=Eaxis;
        Bplot=((draw_i == 1)*Data3d.data3di +(draw_i == 3)*Data3d.data3df);
        z(:,:)=Bplot(tstep,:,ii,mubgrid,:);
        thre = max(max(abs(z)))/25; % Set a threshold, e.g. 0.1, to determine which values should be ignored  
        z(abs(z) < thre) = NaN; % Set values close to 0 to NaN 
     
        
        
        [C, h1] = contourf(x, y,z, 100);hold on;colormap('jet');set(h1, 'linecolor', 'none');hold on
        cb=colorbar('northoutside');
        cb.FontSize = 15;
        cb.Label.FontName = 'Times New Roman';
        cb.FontWeight = 'bold';
        cb.LineWidth = 1;
        %title([char(allname(ii)),'_ ',char(parkind(draw_i))],'fontsize',30);
        %caxis([min(min(z)) max(max(z))])
             
        %saveas(gca,'./pzeta_E_fix_mu','epsc');
        %axis tight
        caxis([min(min(z)) max(max(z))])
   
    case 4
        clear x y z
        lambdaaxis1=(-1:lambdagrids-1)/(((draw_i == 1) +(draw_i == 3))*lam_inv*lambdagrids);
        mub_in=Eperp_in;
        mub=mub_in*((draw_i == 1) /1.86317834+(draw_i == 3) /23.2529351)*mub_inv*mubgrids;
        mubgrid=ceil(mub);
        x=Pzetaaxis;
        y=lambdaaxis1;
        Bplot=((draw_i == 1)*Data3d.data3di +(draw_i == 3)*Data3d.data3df);
        z(:,:)=Bplot(tstep,:,ii,mubgrid,:);
thre = max(max(abs(z)))/50; % Set a threshold, e.g. 0.1, to determine which values should be ignored  
        z(abs(z) < thre) = NaN; % Set values close to 0 to NaN  
        [C, h] = contourf(x, y,z, 100);hold on;colormap('jet');set(h, 'linecolor', 'none');
        cb=colorbar('northoutside');
        cb.FontSize = 15;
        cb.Label.FontName = 'Times New Roman';
        cb.FontWeight = 'bold';
        cb.LineWidth = 1;
        %title([char(allname(ii)),'_ ',char(parkind(draw_i))],'fontsize',30);xlabel('E/keV','fontsize',15);ylabel('\lambda','fontsize',15);
        

         caxis([min(min(z)) max(max(z))])
        
end
%cb = colorbar;
