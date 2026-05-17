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
pdum=zpart(1,m);
tdum=zpart(2,m);
zdum=zpart(3,m);
mu=zpart(6,m)*zpart(6,m);
g=qdspline.spline1d(0,pdum,eq.lsp,eq.dpsi,eq.gpsi);
ri=qdspline.spline1d(0,pdum,eq.lsp,eq.dpsi,eq.ipsi);
q=qdspline.spline1d(0,pdum,eq.lsp,eq.dpsi,eq.qpsi);
b=qdspline.spline2d(0,pdum,tdum,eq.lsp,eq.lst,eq.dpsi,eq.dtheta,eq.bsp);
upara=zpart(4,m)*b*cmratio;
b_inv=1.0/b;
energy=mu*b+0.5*apart*upara^2;
pzeta=g*zpart(4,m)-pdum;
pitch=mu/energy;
energy_prime=energy-(ant.omega/ant.n)*qpart*pzeta;
theta_prime=(ant.n*zdum-ant.omega*tstep*(istep));
int_theta=floor(theta_prime/2/pi);
mod_theta=mod(theta_prime,2*pi);
if(abs(int_theta-int_theta0)>0.01)
    if(mod_theta<mod_theta0)
        large=zpart(:,m);
        small=zpart0(:,m);
    else
        large=zpart0(:,m);
        small=zpart(:,m);
    end
    psii=large(1)-(large(3)-(2*pi))*(small(1)-large(1))/(small(3)+(2*pi)-large(3));
    u_psii=large(4)-(large(3)-(2*pi))*(small(4)-large(4))/(small(3)+(2*pi)-large(3));
    u_thetai=large(6)-(large(3)-(2*pi))*(small(6)-large(6))/(small(3)+(2*pi)-large(3));
    p_zetai=large(8)-(large(3)-(2*pi))*(small(8)-large(8))/(small(3)+(2*pi)-large(3));
    Ei=large(12)-(large(3)-(2*pi))*(small(12)-large(12))/(small(3)+(2*pi)-large(3));

    
    if(abs(zpart(2,m)-zpart0(2,m))<pi)
        thetai=large(2)-(large(3)-(2*pi))*(small(2)-large(2))/(small(3)+(2*pi)-large(3));
    end
    if(abs(zpart(2,m)-zpart0(2,m))>pi)
        if(small(2)>large(2))then
            thetai=small(2)-small(3)*(small(2)-large(2)-(2*pi))/(small(3)+(2*pi)-large(3));
        elseif(small(2)<large(2))
            thetai=small(2)+(2*pi)-small(3)*(small(2)+(2*pi)-large(2))/(small(3)+(2*pi)-large(3));
        end
    end
    thetai=mod(thetai,2*pi);

row_now = size(Poini.psi, 1) + 1;
Poini.psi(row_now, m)      = psii;
Poini.u_para(row_now, m)   = u_psii;
Poini.theta(row_now, m)    = thetai;
Poini.p_zeta(row_now, m)   = p_zetai;
Poini.E(row_now, m)        = Ei;
% Poini.psi(:,m)=[Poini.psi(:,m) ,psii];
% Poini.u_para(:,m)=[Poini.u_para(:,m) ,u_psii];
% Poini.theta(:,m)=[Poini.theta(:,m) ,thetai];
% Poini.p_zeta(:,m)=[Poini.p_zeta(:,m) ,p_zetai];
% Poini.E(:,m)=[Poini.E(:,m) ,Ei];
%Poini.u_para(:,m)=[Poini.u_para(:,m) u_psii];
%Poini.u_para(:,m)=[Poini.u_para(:,m) u_psii];
end

