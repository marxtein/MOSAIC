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
function A = readEFITgfile(gfile)

% open file for reading
findex = fopen(gfile,'r');
if findex < 0 %error
    error(['Could not open the gfile ', gfile])
end

% read data

%     header, pest parameter, # horizontal R grid points, # vertical z grid points
%      READ(kueql,1000)(etitl(i),i=1,5),date,ipestg,nx,nz
      A.case = fscanf(findex,'%48c',1);
      A.idum = fscanf(findex,'%d',1);
      A.nw = fscanf(findex,'%d',1);
      A.nh = fscanf(findex,'%d',1);

%     horz size (m) of computational box, vert size (m) of computational box, R of B location, R of boundary, middle z value
%      READ(kueql,1010) xdim,zdim,rcnt,redge,zmid
      A.rdim = fscanf(findex, '%e', 1);
      A.zdim = fscanf(findex, '%e', 1);
      A.rcentr = fscanf(findex, '%e', 1);
      A.rleft = fscanf(findex, '%e', 1);
      A.zmid = fscanf(findex, '%e', 1);

%     major radius (at magnetic axis), z at mag axis, polflux at mag axis, polflux at limiter, B at specified R location (see above)
%      READ(kueql,1010) xma,zma,psimax,psilim,btor
      A.rmaxis = fscanf(findex, '%e', 1);
      A.zmaxis = fscanf(findex, '%e', 1);
      A.simag = fscanf(findex, '%e', 1);
      A.sibry = fscanf(findex, '%e', 1);
      A.bcentr = fscanf(findex, '%e', 1);

%     current in Amps, 
%      READ(kueql,1010) totcur,psimx(1),psimx(2),xax(1),xax(2)
      A.current = fscanf(findex, '%e', 1);
      tmp = fscanf(findex, '%e', 2);
      tmp = fscanf(findex, '%e', 2);

%      READ(kueql,1010) zax(1),zax(2),psisep,xsep,zsep
      tmp = fscanf(findex, '%e', 2);
      tmp = fscanf(findex, '%e', 1);
      tmp = fscanf(findex, '%e', 1);
      A.xdum = fscanf(findex, '%e', 1);

      %poloidal current function... R*B_toroidal (Wb/m)
%      READ(kueql,1010) (sf(i),     i = 1,npv)
      A.fpol = fscanf(findex, '%e', A.nw);

      %plasma pressure in Joule/m^3
%      READ(kueql,1010) (sp(i),     i = 1,npv)
      A.pres = fscanf(findex, '%e', A.nw);

      %f df/dpsi
%      READ(kueql,1010) (sffp(i),   i = 1,npv)
      A.ffprim_raw = fscanf(findex, '%e', A.nw);

      % dp/dpsi
%      READ(kueql,1010) (spp(i),    i = 1,npv)
      A.pprime_raw = fscanf(findex, '%e', A.nw);

      % polflux on rectangular grid
%      READ(kueql,1010) ((psi(i,j), i = 1,nx), j = 1,nz)
      tmppsi = fscanf(findex, '%e', A.nw*A.nh);
      A.psizr = reshape(tmppsi,A.nw,A.nh);

      % q values on flux surfaces
%      READ(kueql,1010) (efitq(i),  i = 1,npv)
      A.qpsi = fscanf(findex, '%e', A.nw);

      % number of boundary points, and number of limiter points
%      READ(kueql,1020) neftbd,neftlm
      A.nbbbs = fscanf(findex, '%d', 1);
      A.limitr = fscanf(findex, '%d', 1);

      % plot the boundary
%      READ(kueql,1010) (efitbdy(i),i = 1,nefbdy)
      tmpbdy = fscanf(findex, '%e', 2*A.nbbbs);
      tmpbdy = reshape(tmpbdy,2,A.nbbbs);
      A.rbbbs = tmpbdy(1,:);
      A.zbbbs = tmpbdy(2,:);

      % plot the limiter
%      READ(kueql,1010) (efitlim(i),i = 1,neflim)
      tmplim = fscanf(findex, '%e', 2*A.limitr);
      tmplim = reshape(tmplim,2,A.limitr);
      A.rlim = tmplim(1,:);
      A.zlim = tmplim(2,:);

%     check some switches, and characteristic radius of rotation
      A.kvtor = fscanf(findex, '%d', 1);
      A.rvtor = fscanf(findex, '%e', 1);
      A.nmass = fscanf(findex, '%d', 1);

%     rotation related quantities; optional
      if A.kvtor > 0
        A.pressw = fscanf(findex, '%e', A.nw);
        A.pwprim = fscanf(findex, '%e', A.nw);
      end

%     mass related properties; optional
      if A.nmass > 0
        A.dmion = fscanf(findex, '%e', A.nw);
      end

%     toroidal flux
      A.rhovn = fscanf(findex, '%e', A.nw);
      
      fclose(findex);

