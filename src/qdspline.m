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
%Contains: constructSpline1d, constructSpline2d, spline1d, spline2d, matrix2d_spline2d, matrix3d_spline2d
%DONT use iflag=1 for the first point.
classdef qdspline
    methods(Static)
        
        function f = constructSpline1d(iflag, nsp, delx, y)
            %iflag=0: y=y1+y2*x+y3*x^2 for the first point
            %iflag=1: y=y1+y2*sqrt(x)+y3*x for the first point  DONT use this!
            ftmp=y;
            if iflag == 0
                
                ftmp(2,1) = (4.0*ftmp(1,2)-ftmp(1,3)-3.0*ftmp(1,1))/(2.0*delx);  % y=y1+y2*x+y3*x^2
                ftmp(3,1) = (ftmp(1,2)-ftmp(1,1)-ftmp(2,1)*delx)/(delx*delx);
                
                for i=2:nsp-2
                    ipp = min(i+2, nsp);
                    ftmp(2,i) = -ftmp(2,i-1) + 2.0*(ftmp(1,i)-ftmp(1,i-1))/delx;
                    ftmp(1,i+1) = 0.5*delx*ftmp(2,i) + 0.25*ftmp(1,ipp) + 0.75*ftmp(1,i); % smooth
                end
                
                ftmp(2,nsp-1) = -ftmp(2,nsp-2) + 2.0*(ftmp(1,nsp-1)-ftmp(1,nsp-2))/delx;
                ftmp(2,nsp) = -ftmp(2,nsp-1) + 2.0*(ftmp(1,nsp)-ftmp(1,nsp-1))/delx;
                
                for i=2:nsp-1
                    ftmp(3,i) = (ftmp(2,i+1)-ftmp(2,i))/(2.0*delx);
                end
                
            elseif iflag == 1
                ftmp(2,1) = (2.0*ftmp(1,2)-ftmp(1,3)-ftmp(1,1))/((2.0-sqrt(2.0))*sqrt(delx));
                ftmp(3,1) = (ftmp(1,2)-ftmp(1,1)-ftmp(2,1)*sqrt(delx))/delx;
                
                ftmp(2,2) = 0.5*ftmp(2,1)/sqrt(delx) + ftmp(3,1);
                ftmp(3,2) = (ftmp(1,3)-ftmp(1,2)-delx*ftmp(2,2))/(delx*delx);
                
                for i=3:nsp-2
                    ipp=min(i+2,nsp);
                    ftmp(2,i) = -ftmp(2,i-1) + 2.0*(ftmp(1,i)-ftmp(1,i-1))/delx;
                    % smooth
                    f(1,i+1)=0.5*delx*ftmp(2,i) + 0.25*ftmp(1,ipp) + 0.75*ftmp(1,i);
                end
                
                ftmp(2,nsp-1) = -ftmp(2,nsp-2) + 2.0*(ftmp(1,nsp-1)-ftmp(1,nsp-2))/delx;
                ftmp(2,nsp) = -ftmp(2,nsp-1) + 2.0*(ftmp(1,nsp)-ftmp(1,nsp-1))/delx;
                
                for i=3:nsp-1
                    ftmp(3,i) = (ftmp(2,i+1)-ftmp(2,i))/(2.0*delx);
                end
                ftmp(3,nsp) = 0.0;
            end
            f=ftmp;
        end
        
        function f = constructSpline2d(iflagx, iflagy, nx, ny, delx, dely, z)
            
            ftmp = z;
            for j=1:ny
                ddum1(1,:) = ftmp(1,:,j);
                ddum1(2:3,:) =0.0;
                ddum2 = qdspline.constructSpline1d(iflagx, nx, delx, ddum1);
                
                ftmp(1,:,j) = ddum2(1,:);
                ftmp(2,:,j) = ddum2(2,:);
                ftmp(3,:,j) = ddum2(3,:);
            end
            
            for i=1:nx
                for s=1:3
                    ddum3(1,:) = ftmp(s,i,:);
                    ddum3(2:3,:) = 0.0;
                    ddum4 = qdspline.constructSpline1d(iflagy, ny, dely, ddum3);
                    ftmp(s,i,:) = ddum4(1,:);
                    ftmp(s+3,i,:) = ddum4(2,:);
                    ftmp(s+6,i,:) = ddum4(3,:);
                end
            end
            f=ftmp;
        end
        
        
        function f = spline1d(iflag, pdum, nsp, delx, y)
            %iflag=0: value;
            %iflag=1: first derivative
            %iflag=2: second derivative
            i = max(1,min(nsp-1,ceil(pdum/delx)));
            dpx = pdum-delx*(i-1);      % allow pdum to be larger than the maximum equilibrium value
            dp2 = dpx*dpx;
            
            if iflag == 0
                f = y(1,i) + dpx*y(2,i) + dp2*y(3,i);
            elseif iflag == 1
                f = y(2,i) + 2*dpx*y(3,i);
            elseif iflag == 2
                f = 2*y(3,i);
            end
        end
        
        function f = spline2d(deriv, x, y, nx, ny, delx, dely, z)
            %deriv=0: f;
            %deriv=1: \partial f/\partial x;
            %deriv=2: \partial f\partial y;
            %deriv=3: \partial^f/(\partial x \partial y)
            i = max(1,min(nx-1,ceil(x/delx)));
            j = max(1,min(ny-1,ceil(y/dely)));
            
            dx = x-delx*(i-1);% allow pdum to be larger than the maximum equilibrium value
            dy = y-dely*(j-1);% allow pdum to be larger than the maximum equilibrium value
            
            dx2=dx.*dx;
            dy2=dy.*dy;
            
            if deriv==0
                f = z(1,i,j) + z(2,i,j).*dx + z(3,i,j).*dx2...
                    + (z(4,i,j) + z(5,i,j).*dx + z(6,i,j).*dx2).*dy...
                    + (z(7,i,j) + z(8,i,j).*dx + z(9,i,j).*dx2).*dy2;
            elseif deriv==1
                f = z(2,i,j) + 2.0*dx*z(3,i,j)...
                    + (z(5,i,j) + 2.0*dx*z(6,i,j))*dy...
                    + (z(8,i,j) + 2.0*dx*z(9,i,j))*dy2;
            elseif deriv==2
                f = z(4,i,j) + z(5,i,j)*dx + z(6,i,j)*dx2...
                    + 2.0*dy*(z(7,i,j) + z(8,i,j)*dx + z(9,i,j)*dx2);
            elseif deriv==3
                f = z(5,i,j) + 2.0*dx*z(6,i,j)...
                    + 2.0*dy*(z(8,i,j) + 2.0*dx*z(9,i,j));
            end
        end
        
        
        function f = matrix2d_spline2d(deriv, x, y, nx, ny, delx, dely, z)
            %deriv=0: f;
            %deriv=1: \partial f/\partial x;
            %deriv=2: \partial f\partial y;
            %deriv=3: \partial^f/(\partial x \partial y)
            ix = max(1,min(nx-1,ceil(x/delx)));
            jy = max(1,min(ny-1,ceil(y/dely)));
            
            dx = x-delx*(ix-1);
            dy = y-dely*(jy-1);
            
            dx2=dx.*dx;
            dy2=dy.*dy;
            
            index_x=max(size(x));
            index_y=max(size(y));
            
            if deriv==0
                
                for j=1:index_y
                    for i =1:index_x
                        f(i,j) = z(1,ix(i),jy(j)) + z(2,ix(i),jy(j)).*dx(i) + z(3,ix(i),jy(j)).*dx2(i)...
                            + (z(4,ix(i),jy(j)) + z(5,ix(i),jy(j)).*dx(i) + z(6,ix(i),jy(j)).*dx2(i))*dy(j)...
                            + (z(7,ix(i),jy(j)) + z(8,ix(i),jy(j)).*dx(i) + z(9,ix(i),jy(j)).*dx2(i))*dy2(j);
                    end
                end
                
            elseif deriv==1
                
                for j=1:index_y
                    for i =1:index_x
                        f(i,j) = z(2,ix(i),jy(j)) + 2.0*dx(i)*z(3,ix(i),jy(j))...
                            + (z(5,ix(i),jy(j)) + 2.0*dx(i)*z(6,ix(i),jy(j)))*dy(j)...
                            + (z(8,ix(i),jy(j)) + 2.0*dx(i)*z(9,ix(i),jy(j)))*dy2(j);
                    end
                end
                
            elseif deriv==2
                
                for j=1:index_y
                    for i =1:index_x
                        f(i,j) = z(4,ix(i),jy(j)) + z(5,ix(i),jy(j))*dx(i) + z(6,ix(i),jy(j))*dx2(i)...
                            + 2.0*dy(j)*(z(7,ix(i),jy(j)) + z(8,ix(i),jy(j))*dx(i) + z(9,ix(i),jy(j))*dx2(i));
                    end
                end
                
            elseif deriv==3
                
                for j=1:index_y
                    for i =1:index_x
                        f(i,j) = z(5,ix(i),jy(j)) + 2.0*dx(i)*z(6,ix(i),jy(j))...
                            + 2.0*dy(j)*(z(8,ix(i),jy(j)) + 2.0*dx(i)*z(9,ix(i),jy(j)));
                    end
                end
                
            end
        end      
        
        % pass x,y as matrix, output f as matrix, y(mpsi,mtheta,mtoroidal) is boozer theta for field aligned mesh
        function f = matrix3d_spline2d(deriv, x, y, nx, ny, delx, dely, z)
            %deriv=0: f;
            %deriv=1: \partial f/\partial x;
            %deriv=2: \partial f\partial y;
            %deriv=3: \partial^f/(\partial x \partial y)
            ix = max(1,min(nx-1,ceil(x/delx)));
            jy = max(1,min(ny-1,ceil(y/dely)));
            
            
            dx = x-delx*(ix-1);
            dy = y-dely*(jy-1);
            
            dx2=dx.*dx;
            dy2=dy.*dy;
            
            index=size(y);
            
            if deriv == 0
                for k=1:index(3)
                    for j=1:index(2)
                        for i =1:index(1)
                            f(i,j,k) = z(1,ix(i),jy(i,j,k)) + z(2,ix(i),jy(i,j,k))*dx(i) + z(3,ix(i),jy(i,j,k))*dx2(i)...
                                + (z(4,ix(i),jy(i,j,k)) + z(5,ix(i),jy(i,j,k))*dx(i) + z(6,ix(i),jy(i,j,k))*dx2(i))*dy(i,j,k)...
                                + (z(7,ix(i),jy(i,j,k)) + z(8,ix(i),jy(i,j,k))*dx(i) + z(9,ix(i),jy(i,j,k))*dx2(i))*dy2(i,j,k);
                        end
                    end
                end
                
            elseif deriv == 1
                for k=1:index(3)
                    for j=1:index(2)
                        for i =1:index(1)
                            f(i,j,k) = z(2,ix(i),jy(i,j,k)) + 2.0*dx(i)*z(3,ix(i),jy(i,j,k))...
                                + (z(5,ix(i),jy(i,j,k)) + 2.0*dx(i)*z(6,ix(i),jy(i,j,k)))*dy(i,j,k)...
                                + (z(8,ix(i),jy(i,j,k)) + 2.0*dx(i)*z(9,ix(i),jy(i,j,k)))*dy2(i,j,k);
                        end
                    end
                end
            elseif deriv == 2
                for k=1:index(3)
                    for j=1:index(2)
                        for i =1:index(1)
                            f(i,j,k) = (z(4,ix(i),jy(i,j,k)) + z(5,ix(i),jy(i,j,k))*dx(i) + z(6,ix(i),jy(i,j,k))*dx2(i))...
                                + 2.0*dy(i,j,k)*(z(7,ix(i),jy(i,j,k)) + z(8,ix(i),jy(i,j,k))*dx(i) + z(9,ix(i),jy(i,j,k))*dx2(i));
                        end
                    end
                end
            elseif deriv == 3
                for k=1:index(3)
                    for j=1:index(2)
                        for i =1:index(1)
                            f(i,j,k) = (z(5,ix(i),jy(i,j,k)) + 2.0*dx(i)*z(6,ix(i),jy(i,j,k)))...
                                + 2.0*dy(i,j,k)*(z(8,ix(i),jy(i,j,k)) + 2.0*dx(i)*z(9,ix(i),jy(i,j,k)));
                        end
                    end
                end
            end
            
        end
        
    end
end