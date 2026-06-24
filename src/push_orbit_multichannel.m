function B = push_orbit_multichannel(zpart0, mstep, tstep, qpart, apart, ...
                                     eq_curr, eq_gradB, ant_multi, amp_mod)
% PUSH_ORBIT_MULTICHANNEL  Multi-channel perturbed guiding-center pusher.
%
% Streamlined RK2 (midpoint) integrator lifted from poincare_perturb.m,
% with the antenna inner-m loop replaced by a call to eval_perturb_field
% that sums over ant_multi.channels(k). Records E, P_zeta and per-channel
% E'_k = E - (omega_k/n_k)*P_zeta every step.
%
% Inputs
%   zpart0  : initial zpart matrix (14 x np), as in poincare_perturb. Required
%             slots used here: 1=psi, 2=theta, 3=zeta, 4=rho_||, 6=sqrt(mu).
%   mstep   : number of time steps.
%   tstep   : time step in omegacp^-1.
%   qpart   : particle charge (normalized).
%   apart   : particle mass (normalized).
%   eq_curr, eq_gradB : standard MOSAIC physics switches (use 1, 1).
%   ant_multi : multi-channel antenna struct (carries .channels(k)).
%   amp_mod  : amplitude modulator (e.g., 1/20 like poincare_perturb).
%
% Output struct B:
%   B.t       : (mstep x 1)  time in omegacp^-1
%   B.psi     : (mstep x np)
%   B.theta   : (mstep x np)
%   B.zeta    : (mstep x np)
%   B.rho_par : (mstep x np)
%   B.E       : (mstep x np)  total energy (MOSAIC norm; keV via *unit.energy_norm/1000)
%   B.Pzeta   : (mstep x np)  P_zeta (norm by eq.psiw if you want)
%   B.Eprime  : (mstep x np x nc)  E'_k per channel = E - (omega_k/n_k)*qpart*P_zeta
%   B.channels: copy of ant_multi.channels (for reference)
%   B.confine : cell array per-particle 'inside' or 'lost'.

global eq

nc       = numel(ant_multi.channels);
np       = size(zpart0, 2);
cmratio  = qpart/apart;
c_inv    = 1/qpart;

B.t       = (1:mstep)' * tstep;
B.psi     = nan(mstep, np);
B.theta   = nan(mstep, np);
B.zeta    = nan(mstep, np);
B.rho_par = nan(mstep, np);
B.E       = nan(mstep, np);
B.Pzeta   = nan(mstep, np);
B.Eprime  = nan(mstep, np, nc);
B.channels = ant_multi.channels;
B.confine = repmat({'inside'}, 1, np);

omega_over_n = zeros(1, nc);
for k = 1:nc
    omega_over_n(k) = ant_multi.channels(k).omega / ant_multi.channels(k).n;
end

zpart  = zpart0;
zpart0_buf = zpart0;

for m = 1:np
    for istep = 1:mstep
        if zpart(1, m) > eq.psiw
            B.confine{m} = 'lost';
            break;
        end

        for irk = 1:2
            if irk == 1
                dtime = 0.5*tstep;
                zpart0_buf(:, m) = zpart(:, m);
            else
                dtime = tstep;
            end

            pdum = zpart(1, m);
            tdum = zpart(2, m);
            zdum = zpart(3, m);
            mu   = zpart(6, m) * zpart(6, m);

            g    = qdspline.spline1d(0, pdum, eq.lsp, eq.dpsi, eq.gpsi);
            ri   = qdspline.spline1d(0, pdum, eq.lsp, eq.dpsi, eq.ipsi);
            gp   = qdspline.spline1d(1, pdum, eq.lsp, eq.dpsi, eq.gpsi);
            rip  = qdspline.spline1d(1, pdum, eq.lsp, eq.dpsi, eq.ipsi);
            q    = qdspline.spline1d(0, pdum, eq.lsp, eq.dpsi, eq.qpsi);
            deni = 1/(g*q + ri);

            b    = qdspline.spline2d(0, pdum, tdum, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            dbdp = qdspline.spline2d(1, pdum, tdum, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
            dbdt = qdspline.spline2d(2, pdum, tdum, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);

            upara  = zpart(4, m) * b * cmratio;
            dedb   = zpart(4, m)*zpart(4, m)*b*cmratio + c_inv*mu;
            b_inv  = 1.0/b;
            energy = mu*b + 0.5*apart*upara^2;
            pzeta  = g*zpart(4, m) - pdum;

            t_now = tstep * (istep + 0.5*irk);
            F = eval_perturb_field(ant_multi.channels, pdum, tdum, zdum, t_now, amp_mod);

            apara    = F.apara;
            dptdp    = F.dptdp;  dptdt    = F.dptdt;  dptdz    = F.dptdz;
            dapdp    = F.dapdp;  dapdt    = F.dapdt;  dapdz    = F.dapdz;
            paparapt = F.paparapt;

            lam     = apara*b_inv;
            dlamdp  = (dapdp - lam*dbdp)*b_inv;
            dlamdt  = (dapdt - lam*dbdt)*b_inv;
            dlamdz  = dapdz*b_inv;
            plampt  = paparapt*b_inv;

            rdot = ((gp*(zpart(4, m)*eq_curr + lam) + g*dlamdp - 1)*(dedb*dbdt + dptdt))*deni ...
                 - ((q + (zpart(4, m) + lam)*rip + ri*dlamdp)*(dptdz))*deni ...
                 + ((ri*dlamdz - g*dlamdt)*(dedb*dbdp + dptdp))*deni - plampt;

            pdot = eq_gradB*(-dedb*g*dbdt)*deni ...
                 + (ri*dptdz - g*dptdt)*deni ...
                 + upara*b*(g*dlamdt - ri*dlamdz)*deni;

            tdot = (upara*b*(1 - (zpart(4, m)*eq_curr + lam)*gp - g*dlamdp) ...
                  + g*(dedb*dbdp*eq_gradB + dptdp))*deni;

            zdot = (upara*b*(q + (zpart(4, m)*eq_curr + lam)*rip + ri*dlamdp) ...
                  - ri*(dedb*dbdp*eq_gradB + dptdp))*deni;

            zpart(1, m) = zpart0_buf(1, m) + pdot*dtime;
            zpart(2, m) = zpart0_buf(2, m) + tdot*dtime;
            zpart(3, m) = zpart0_buf(3, m) + zdot*dtime;
            zpart(4, m) = zpart0_buf(4, m) + rdot*dtime;

            zpart(2, m) = mod(zpart(2, m), 2*pi);
            zpart(3, m) = mod(zpart(3, m), 2*pi);
        end

        % Diagnose at end of step
        pdum = zpart(1, m); tdum = zpart(2, m); zdum = zpart(3, m);
        mu   = zpart(6, m)*zpart(6, m);
        g    = qdspline.spline1d(0, pdum, eq.lsp, eq.dpsi, eq.gpsi);
        b    = qdspline.spline2d(0, pdum, tdum, eq.lsp, eq.lst, eq.dpsi, eq.dtheta, eq.bsp);
        upara  = zpart(4, m)*b*cmratio;
        energy = mu*b + 0.5*apart*upara^2;
        pzeta  = g*zpart(4, m) - pdum;

        B.psi(istep, m)     = pdum;
        B.theta(istep, m)   = tdum;
        B.zeta(istep, m)    = zdum;
        B.rho_par(istep, m) = zpart(4, m);
        B.E(istep, m)       = energy;
        B.Pzeta(istep, m)   = pzeta;
        for k = 1:nc
            B.Eprime(istep, m, k) = energy - omega_over_n(k)*qpart*pzeta;
        end
    end
end
end
