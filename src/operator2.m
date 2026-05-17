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
%Contains: gradient1d/d2dx1d/d4dx1d  (1d vector array input,output L square matrix),
%          gradient2d/d2dx2d/d4dx2d  (2 1d arrays input,output Lx(nx*nx), Ly(ny*ny) and Lxbig((nx*ny)*(nx*ny)))
classdef operator2
    methods(Static)
        %%
        function f = gradient1d(iflag, x)
            %iflag=0: Dirichlet boundary condition at x=0 and x=L
            %iflag=1: periodic boundary condition
            %iflag=2: Dirichlet boundary condition at x=0 and Neumann boundary condition at x=L
            %iflag=3: extraplation values at boundary grids for equilibrium quantity gradient
            
            n = max(size(x));
            L = sparse(n,n);
            dx = diff(x);
            if iflag==0
                
                a1 = 1;
                L = L ...
                    + sparse([1],[1],[a1],n,n);
                
                a2 = -(dx(1)-dx(2))/(dx(1)*dx(2));
                a3 = dx(1)/(dx(2)*(dx(1)+dx(2)));
                L = L ...
                    + sparse([2 2],[2 3],[a2 a3],n,n);
                
                for i=3:n-2
                    a(i-1) = -dx(i)/(dx(i-1)*(dx(i-1)+dx(i)));
                    a(i) = -(dx(i-1)-dx(i))/(dx(i-1)*dx(i));
                    a(i+1) = dx(i-1)/(dx(i)*(dx(i)+dx(i-1)));
                    L = L ...
                        + sparse([i i i],[i-1 i i+1],[a(i-1) a(i) a(i+1)],n,n);
                end
                
                an2 = -dx(n-1)/(dx(n-2)*(dx(n-1)+dx(n-2)));
                an1 = (dx(n-1)-dx(n-2))/(dx(n-1)*dx(n-2));
                
                L = L ...
                    + sparse([n-1 n-1],[n-2 n-1],[an2 an1],n,n);
                
                an = 1;
                L = L ...
                    + sparse([n],[n],[an],n,n);
                
                
            elseif iflag==1
                
                an1 = -dx(1)/(dx(n-1)*(dx(n-1)+dx(1)));
                a1 = -(dx(n-1)-dx(1))/(dx(n-1)*dx(1));
                a2 = dx(n-1)/(dx(1)*(dx(1)+dx(n-1)));
                L = L ...
                    + sparse([1 1 1],[n-1 1 2],[an1 a1 a2],n,n);
                
                for i=2:n-1
                    a(i-1) = -dx(i)/(dx(i-1)*(dx(i-1)+dx(i)));
                    a(i) = -(dx(i-1)-dx(i))/(dx(i-1)*dx(i));
                    a(i+1) = dx(i-1)/(dx(i)*(dx(i)+dx(i-1)));
                    L = L...
                        + sparse([i i i],[i-1 i i+1],[a(i-1) a(i) a(i+1)],n,n);
                end
                
                an1 = -dx(1)/(dx(n-1)*(dx(n-1)+dx(1)));
                an = -(dx(n-1)-dx(1))/(dx(n-1)*dx(1));
                a2 = dx(n-1)/(dx(1)*(dx(1)+dx(n-1)));
                L = L ...
                    + sparse([n n n],[2 n-1 n],[a2 an1 an],n,n);
                
            elseif iflag==2
                %f'(x1)=0(i.e.,f(x1)=f(x2)), f(xn)=0
                a1 = 1;
                L = L ...
                    + sparse([1],[1],[a1],n,n);
                
                a2 = -dx(1)/(dx(2)*(dx(1)+dx(2)));
                a3 = dx(1)/(dx(2)*(dx(1)+dx(2)));
                L = L ...
                    + sparse([2 2],[2 3],[a2 a3],n,n);
                
                for i=3:n-2
                    a(i-1) = -dx(i)/(dx(i-1)*(dx(i-1)+dx(i)));
                    a(i) = -(dx(i-1)-dx(i))/(dx(i-1)*dx(i));
                    a(i+1) = dx(i-1)/(dx(i)*(dx(i)+dx(i-1)));
                    L = L ...
                        + sparse([i i i],[i-1 i i+1],[a(i-1) a(i) a(i+1)],n,n);
                end
                
                an2 = -dx(n-1)/(dx(n-2)*(dx(n-1)+dx(n-2)));
                an1 = (dx(n-1)-dx(n-2))/(dx(n-1)*dx(n-2));
                
                L = L ...
                    + sparse([n-1 n-1],[n-2 n-1],[an2 an1],n,n);
                
                an = 1;
                L = L ...
                    + sparse([n],[n],[an],n,n);
                
            elseif iflag==3
                
                a1 = -(2*dx(1)+dx(2))/(dx(1)*(dx(1)+dx(2)));
                a2 = (dx(1)+dx(2))/(dx(1)*dx(2));
                a3 = -dx(1)/(dx(2)*(dx(1)+dx(2)));
                L = L ...
                    + sparse([1 1 1],[1 2 3],[a1 a2 a3],n,n);
                
                for i=2:n-1
                    a(i-1) = -dx(i)/(dx(i-1)*(dx(i-1)+dx(i)));
                    a(i) = -(dx(i-1)-dx(i))/(dx(i-1)*dx(i));
                    a(i+1) = dx(i-1)/(dx(i)*(dx(i)+dx(i-1)));
                    L = L ...
                        + sparse([i i i],[i-1 i i+1],[a(i-1) a(i) a(i+1)],n,n);
                end
                
                
                an2 = dx(n-1)/(dx(n-2)*(dx(n-1)+dx(n-2)));
                an1 = -(dx(n-1)+dx(n-2))/(dx(n-1)*dx(n-2));
                an  = (2*dx(n-1)+dx(n-2))/(dx(n-1)*(dx(n-1)+dx(n-2)));
                L = L ...
                    + sparse([n n n],[n-2 n-1 n],[an2 an1 an],n,n);
            end
            
            f = L;
        end
        
        function [Lx,Ly,Lxbig] = gradient2d(iflagx, iflagy, x, y)
            % iflagx=0/1/2: x use Dirichlet/periodic/Neumann boundary condition
            % iflagy=0/1/2: y use Dirichlet/periodic/Neumann boundary condition
            % x is the first dimension coordinate
            % y is the second dimension coordinate
            
            nx = max(size(x));
            ny = max(size(y));
            Lx = operator2.gradient1d(iflagx, x);
            Ly = operator2.gradient1d(iflagy, y);
            
            Lx_cell = repmat({Lx},1,ny);
            Lxbig = blkdiag(Lx_cell{:});
        end
        %%
        function f=d2dx1d(iflag, x)
            %iflag=0: Dirichlet boundary condition at both x=0 and x=L
            %iflag=1: periodic boundary condition
            %iflag=2: Neumann boundary condition at x=0 and Dirichlet boundary at x=L
            %iflag=3: extraplation values at boundary grids for equilibrium quantity laplacian
            n = max(size(x));
            L = sparse(n,n);
            dx = diff(x);
            if iflag==0
                %f(x1)=0, f(xn)=0
                a1 = 1;
                L = L ...
                    + sparse([1],[1],[a1],n,n);
                
                a2 = -2/(dx(1)*dx(2));
                a3 = (2/(dx(1)*dx(2)))*(dx(1)/(dx(1)+dx(2)));
                L = L ...
                    + sparse([2 2],[2 3],[a2 a3],n,n);
                
                for i=3:n-2
                    a(i-1) = 2*dx(i)/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                    a(i) = -2/(dx(i)*dx(i-1));
                    a(i+1) = 2*dx(i-1)/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                    L = L ...
                        + sparse([i i i],[i-1 i i+1],[a(i-1) a(i) a(i+1)],n,n);
                end
                
                an2 = (2/(dx(n-1)*dx(n-2)))*(dx(n-1)/(dx(n-1)+dx(n-2))); %a_{n-2}
                an1 = -2/(dx(n-1)*dx(n-2)); %a_{n-1}
                L = L ...
                    + sparse([n-1 n-1],[n-2 n-1],[an2 an1],n,n);
                
                an = 1;
                L = L ...
                    + sparse([n],[n],[an],n,n);
                
            elseif iflag==1
                
                an1 = 2*dx(1)/(dx(1)*dx(n-1)*(dx(1)+dx(n-1)));
                a1 = -2/(dx(1)*dx(n-1));
                a2 = 2*dx(n-1)/(dx(1)*dx(n-1)*(dx(1)+dx(n-1)));
                L = L ...
                    + sparse([1 1 1],[n-1 1 2],[an1 a1 a2],n,n);
                
                for i=2:n-1
                    a(i-1) = 2*dx(i)/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                    a(i) = -2/(dx(i)*dx(i-1));
                    a(i+1) = 2*dx(i-1)/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                    L = L ...
                        + sparse([i i i],[i-1 i i+1],[a(i-1) a(i) a(i+1)],n,n);
                end
                
                an1 = 2*dx(1)/(dx(1)*dx(n-1)*(dx(1)+dx(n-1)));
                an = -2/(dx(1)*dx(n-1));
                a2 = 2*dx(n-1)/(dx(1)*dx(n-1)*(dx(1)+dx(n-1)));
                L = L ...
                    + sparse([n n n],[2 n-1 n],[a2 an1 an],n,n);
            elseif iflag==2
                %f'(x1)=0(i.e.,f(x1)=f(x2)), f(xn)=0
                a1 = 1.0;
                L = L ...
                    + sparse(1,1,a1,n,n);
                
                a2 = -(2/(dx(1)*dx(2)))*(dx(1)/(dx(1)+dx(2)));
                a3 = (2/(dx(1)*dx(2)))*(dx(1)/(dx(1)+dx(2)));
                L = L ...
                    + sparse([2 2],[2 3],[a2 a3],n,n);
                
                for i=3:n-2
                    a(i-1) = 2*dx(i)/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                    a(i) = -2/(dx(i)*dx(i-1));
                    a(i+1) = 2*dx(i-1)/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                    L = L ...
                        + sparse([i i i],[i-1 i i+1],[a(i-1) a(i) a(i+1)],n,n);
                end
                
                an2 = (2/(dx(n-1)*dx(n-2)))*(dx(n-1)/(dx(n-1)+dx(n-2))); %a_{n-2}
                an1 = -2/(dx(n-1)*dx(n-2)); %a_{n-1}
                L = L ...
                    + sparse([n-1 n-1],[n-2 n-1],[an2 an1],n,n);
                
                an = 1;
                L = L ...
                    + sparse([n],[n],[an],n,n);
                
            elseif iflag==3
                a1 = (2/(dx(1)*dx(2)))*(dx(2)/(dx(1)+dx(2)));
                a2 = -(2/(dx(1)*dx(2)));
                a3 = (2/(dx(1)*dx(2)))*(dx(1)/(dx(1)+dx(2)));
                L = L ...
                    + sparse([1 1 1],[1 2 3],[a1 a2 a3],n,n);
                
                for i=2:n-1
                    a(i-1) = 2*dx(i)/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                    a(i) = -2/(dx(i)*dx(i-1));
                    a(i+1) = 2*dx(i-1)/(dx(i)*dx(i-1)*(dx(i)+dx(i-1)));
                    L = L ...
                        + sparse([i i i],[i-1 i i+1],[a(i-1) a(i) a(i+1)],n,n);
                end
                
                an2 = (2/(dx(n-1)*dx(n-2)))*(dx(n-1)/(dx(n-1)+dx(n-2))); %a_{n-2}
                an1 = -2/(dx(n-1)*dx(n-2)); %a_{n-1}
                an = (2/(dx(n-1)*dx(n-2)))*(dx(n-2)/(dx(n-1)+dx(n-2))); %a_{n}
                
                L = L ...
                    + sparse([n n n],[n-2 n-1 n],[an2 an1 an],n,n);
            end
            f = L;
        end
        
        function [L2x,L2y,L2xbig] = d2dx2d(iflagx, iflagy, x, y)
            % iflagx=0/1/2: L2x use Dirichlet/periodic/Neumann boundary condition
            % iflagy=0/1/2: L2y use Dirichlet/periodic/Neumann boundary condition
            % x is the first dimension coordinate
            % y is the second dimension coordinate
            nx = max(size(x));
            ny = max(size(y));
            
            L2x = operator2.d2dx1d(iflagx, x);
            L2y = operator2.d2dx1d(iflagy, y);
            
            L2x_cell = repmat({L2x},1,ny);
            L2xbig = blkdiag(L2x_cell{:});
        end
                %%
        function f=d3dx1d(x)
            % Only Dirichlet boundary condition at both x=0 and x=L
            
            n = max(size(x));
            L = sparse(n,n);
            dx = diff(x);
            
            %f(x1)=0, f(xn)=0
            a11 = 1;
            L = L ...
                + sparse([1],[1],[a11],n,n);
            
            C2tmp(2) = -(dx(1)-dx(2))/(dx(2)*dx(1));
            C2tmp(3) = dx(1)/(dx(2)*(dx(1)+dx(2)));
            
            A2tmp(2) = 2/(dx(2)*dx(1));
            A2tmp(3) = 2/(dx(3)*dx(2));
            
            B2tmp_diag_high = dx(3)/(dx(3)+dx(2));
            B2tmp(3) = dx(1)/(dx(2)+dx(1));
            B2tmp(4) = dx(2)/(dx(3)+dx(2));
            
            a22 =  -C2tmp(2)*A2tmp(2) + C2tmp(3)*A2tmp(3)*B2tmp_diag_high;
            a23 = C2tmp(2)*A2tmp(2)*B2tmp(3) - C2tmp(3)*A2tmp(3);
            a24 = C2tmp(3)*A2tmp(3)*B2tmp(4);           
            L = L ...
                + sparse([2 2 2],[2 3 4],[a22 a23 a24],n,n);
           
            C3tmp(2) = -dx(3)/(dx(2)*(dx(2)+dx(3)));
            C3tmp(3) = -(dx(2)-dx(3))/(dx(3)*dx(2));
            C3tmp(4) = dx(2)/(dx(3)*(dx(2)+dx(3)));
            
            A3tmp(2) = 2/(dx(2)*dx(1));
            A3tmp(3) = 2/(dx(3)*dx(2));
            A3tmp(4) = 2/(dx(4)*dx(3));
            
            B3tmp(2) = dx(3)/(dx(3)+dx(2));
            B3tmp_diag_low = dx(1)/(dx(2)+dx(1));
            B3tmp_diag_high = dx(4)/(dx(4)+dx(3));
            B3tmp(4) = dx(2)/(dx(3)+dx(2));
            B3tmp(5) = dx(3)/(dx(4)+dx(3));
            
            a32 = -C3tmp(2)*A3tmp(2) + C3tmp(3)*A3tmp(3)*B3tmp(2);
            a33 = C3tmp(2)*A3tmp(2)*B3tmp_diag_low - C3tmp(3)*A3tmp(3) + C3tmp(4)*A3tmp(4)*B3tmp_diag_high;
            a34 = C3tmp(3)*A3tmp(3)*B3tmp(4) - C3tmp(4)*A3tmp(4);
            a35 = C3tmp(4)*A3tmp(4)*B3tmp(5);            
            L = L ...
                + sparse([3 3 3 3],[2 3 4 5],[a32 a33 a34 a35],n,n);
            
            for i=4:n-3
                
                Ctmp(i-1) = -dx(i)/(dx(i-1)*(dx(i-1)+dx(i)));
                Ctmp(i) = -(dx(i-1)-dx(i))/(dx(i)*dx(i-1));
                Ctmp(i+1) = dx(i-1)/(dx(i)*(dx(i-1)+dx(i)));
                
                Atmp(i-1) = 2/(dx(i-1)*dx(i-2));
                Atmp(i) = 2/(dx(i)*dx(i-1));
                Atmp(i+1) = 2/(dx(i+1)*dx(i));
                
                Btmp(i-2) = dx(i-1)/(dx(i-1)+dx(i-2));
                Btmp(i-1) = dx(i)/(dx(i)+dx(i-1));
                Btmp_diag_low = dx(i-2)/(dx(i-1)+dx(i-2));
                Btmp_diag_high = dx(i+1)/(dx(i+1)+dx(i));
                Btmp(i+1) = dx(i-1)/(dx(i)+dx(i-1));
                Btmp(i+2) = dx(i)/(dx(i+1)+dx(i));
                
                a(i-2) = Ctmp(i-1)*Atmp(i-1)*Btmp(i-2);
                a(i-1) = -Ctmp(i-1)*Atmp(i-1) + Ctmp(i)*Atmp(i)*Btmp(i-1);
                a(i) = Ctmp(i-1)*Atmp(i-1)*Btmp_diag_low - Ctmp(i)*Atmp(i) + Ctmp(i+1)*Atmp(i+1)*Btmp_diag_high;
                a(i+1) = Ctmp(i)*Atmp(i)*Btmp(i+1) - Ctmp(i+1)*Atmp(i+1);
                a(i+2) = Ctmp(i+1)*Atmp(i+1)*Btmp(i+2);
                
                L = L ...
                    + sparse([i i i i i],[i-2 i-1 i i+1 i+2],[a(i-2) a(i-1) a(i) a(i+1) a(i+2)],n,n);
            end
                       
            i = n-2;
            Ctmp(i-1) = -dx(i)/(dx(i-1)*(dx(i-1)+dx(i)));
            Ctmp(i) = -(dx(i-1)-dx(i))/(dx(i)*dx(i-1));
            Ctmp(i+1) = dx(i-1)/(dx(i)*(dx(i-1)+dx(i)));
                
            Atmp(i-1) = 2/(dx(i-1)*dx(i-2));
            Atmp(i) = 2/(dx(i)*dx(i-1));
            Atmp(i+1) = 2/(dx(i+1)*dx(i));
            
            Btmp(i-2) = dx(i-1)/(dx(i-1)+dx(i-2));
            Btmp(i-1) = dx(i)/(dx(i)+dx(i-1));
            Btmp_diag_low = dx(i-2)/(dx(i-1)+dx(i-2));
            Btmp_diag_high = dx(i+1)/(dx(i+1)+dx(i));
            Btmp(i+1) = dx(i-1)/(dx(i)+dx(i-1));
                       
            a(i-2) = Ctmp(i-1)*Atmp(i-1)*Btmp(i-2);
            a(i-1) = -Ctmp(i-1)*Atmp(i-1) + Ctmp(i)*Atmp(i)*Btmp(i-1);
            a(i) = Ctmp(i-1)*Atmp(i-1)*Btmp_diag_low - Ctmp(i)*Atmp(i) + Ctmp(i+1)*Atmp(i+1)*Btmp_diag_high;
            a(i+1) = Ctmp(i)*Atmp(i)*Btmp(i+1) - Ctmp(i+1)*Atmp(i+1);          
            L = L ...
                + sparse([i i i i],[i-2 i-1 i i+1],[a(i-2) a(i-1) a(i) a(i+1)],n,n);
                        
            i = n-1;
            Ctmp(i-1) = -dx(i)/(dx(i-1)*(dx(i-1)+dx(i)));
            Ctmp(i) = -(dx(i-1)-dx(i))/(dx(i)*dx(i-1));
            
            Atmp(i-1) = 2/(dx(i-1)*dx(i-2));
            Atmp(i) = 2/(dx(i)*dx(i-1));
            
            Btmp(i-2) = dx(i-1)/(dx(i-1)+dx(i-2));
            Btmp(i-1) = dx(i)/(dx(i)+dx(i-1));
            Btmp_diag_low = dx(i-2)/(dx(i-1)+dx(i-2));
            
            a(i-2) = Ctmp(i-1)*Atmp(i-1)*Btmp(i-2);
            a(i-1) = -Ctmp(i-1)*Atmp(i-1) + Ctmp(i)*Atmp(i)*Btmp(i-1);
            a(i) = Ctmp(i-1)*Atmp(i-1)*Btmp_diag_low - Ctmp(i)*Atmp(i);           
            L = L ...
                + sparse([i i i],[i-2 i-1 i],[a(i-2) a(i-1) a(i)],n,n);
                 
            ann = 1;
            L = L ...
                + sparse([n],[n],[ann],n,n);
            
            f = L;
        end
        
        function [L3x,L3y,L3xbig] = d3dx2d(x, y)
            % L3x use Dirichlet boundary condition
            % L3y use Dirichlet boundary condition
            % x is the first dimension coordinate
            % y is the second dimension coordinate
            nx = max(size(x));
            ny = max(size(y));
            
            L3x = operator2.d3dx1d(x);
            L3y = operator2.d3dx1d(y);
            
            L3x_cell = repmat({L3x},1,ny);
            L3xbig = blkdiag(L3x_cell{:});
        end       
        %%
        function f=d4dx1d(x)
            % Only Dirichlet boundary condition at both x=0 and x=L
            
            n = max(size(x));
            L = sparse(n,n);
            dx = diff(x);
            
            %f(x1)=0, f(xn)=0
            a11 = 1;
            L = L ...
                + sparse([1],[1],[a11],n,n);
            
            A2tmp(2) = 2/(dx(2)*dx(1));
            A2tmp(3) = 2/(dx(3)*dx(2));
            
            B2tmp_diag_high = dx(3)/(dx(3)+dx(2));
            B2tmp(3) = dx(1)/(dx(2)+dx(1));
            B2tmp(4) = dx(2)/(dx(3)+dx(2));
            
            a22 =  A2tmp(2)^2 + A2tmp(2)*B2tmp(3)*A2tmp(3)*B2tmp_diag_high;
            a23 = -A2tmp(2)^2*B2tmp(3) - A2tmp(2)*B2tmp(3)*A2tmp(3);
            a24 = A2tmp(2)*B2tmp(3)*A2tmp(3)*B2tmp(4);           
            L = L ...
                + sparse([2 2 2],[2 3 4],[a22 a23 a24],n,n);
            
            A3tmp(2) = 2/(dx(2)*dx(1));
            A3tmp(3) = 2/(dx(3)*dx(2));
            A3tmp(4) = 2/(dx(4)*dx(3));
            
            B3tmp(2) = dx(3)/(dx(3)+dx(2));
            B3tmp_diag_low = dx(1)/(dx(2)+dx(1));
            B3tmp_diag_high = dx(4)/(dx(4)+dx(3));
            B3tmp(4) = dx(2)/(dx(3)+dx(2));
            B3tmp(5) = dx(3)/(dx(4)+dx(3));
            
            a32 = -A3tmp(3)*B3tmp(2)*A3tmp(2) - A3tmp(3)^2*B3tmp(2);
            a33 = A3tmp(3)*B3tmp(2)*A3tmp(2)*B3tmp_diag_low + A3tmp(3)^2 + A3tmp(3)*B3tmp(4)*A3tmp(4)*B3tmp_diag_high;
            a34 = -A3tmp(3)^2*B3tmp(4) - A3tmp(3)*B3tmp(4)*A3tmp(4);
            a35 = A3tmp(3)*B3tmp(4)*A3tmp(4)*B3tmp(5);            
            L = L ...
                + sparse([3 3 3 3],[2 3 4 5],[a32 a33 a34 a35],n,n);
            
            for i=4:n-3
                
                Atmp(i-1) = 2/(dx(i-1)*dx(i-2));
                Atmp(i) = 2/(dx(i)*dx(i-1));
                Atmp(i+1) = 2/(dx(i+1)*dx(i));
                
                Btmp(i-2) = dx(i-1)/(dx(i-1)+dx(i-2));
                Btmp(i-1) = dx(i)/(dx(i)+dx(i-1));
                Btmp_diag_low = dx(i-2)/(dx(i-1)+dx(i-2));
                Btmp_diag_high = dx(i+1)/(dx(i+1)+dx(i));
                Btmp(i+1) = dx(i-1)/(dx(i)+dx(i-1));
                Btmp(i+2) = dx(i)/(dx(i+1)+dx(i));
                
                a(i-2) = Atmp(i)*Btmp(i-1)*Atmp(i-1)*Btmp(i-2);
                a(i-1) = -Atmp(i)*Btmp(i-1)*Atmp(i-1) - Atmp(i)^2*Btmp(i-1);
                a(i) = Atmp(i)*Btmp(i-1)*Atmp(i-1)*Btmp_diag_low + Atmp(i)^2 + Atmp(i)*Btmp(i+1)*Atmp(i+1)*Btmp_diag_high;
                a(i+1) = -Atmp(i)^2*Btmp(i+1) - Atmp(i)*Btmp(i+1)*Atmp(i+1);
                a(i+2) = Atmp(i)*Btmp(i+1)*Atmp(i+1)*Btmp(i+2);
                
                L = L ...
                    + sparse([i i i i i],[i-2 i-1 i i+1 i+2],[a(i-2) a(i-1) a(i) a(i+1) a(i+2)],n,n);
            end
                       
            i = n-2;
            Atmp(i-1) = 2/(dx(i-1)*dx(i-2));
            Atmp(i) = 2/(dx(i)*dx(i-1));
            Atmp(i+1) = 2/(dx(i+1)*dx(i));
            
            Btmp(i-2) = dx(i-1)/(dx(i-1)+dx(i-2));
            Btmp(i-1) = dx(i)/(dx(i)+dx(i-1));
            Btmp_diag_low = dx(i-2)/(dx(i-1)+dx(i-2));
            Btmp_diag_high = dx(i+1)/(dx(i+1)+dx(i));
            Btmp(i+1) = dx(i-1)/(dx(i)+dx(i-1));
                       
            a(i-2) = Atmp(i)*Btmp(i-1)*Atmp(i-1)*Btmp(i-2);
            a(i-1) = -Atmp(i)*Btmp(i-1)*Atmp(i-1) - Atmp(i)^2*Btmp(i-1);
            a(i) = Atmp(i)*Btmp(i-1)*Atmp(i-1)*Btmp_diag_low + Atmp(i)^2 + Atmp(i)*Btmp(i+1)*Atmp(i+1)*Btmp_diag_high;
            a(i+1) = -Atmp(i)^2*Btmp(i+1) - Atmp(i)*Btmp(i+1)*Atmp(i+1);           
            L = L ...
                + sparse([i i i i],[i-2 i-1 i i+1],[a(i-2) a(i-1) a(i) a(i+1)],n,n);
                        
            i = n-1;
            Atmp(i-1) = 2/(dx(i-1)*dx(i-2));
            Atmp(i) = 2/(dx(i)*dx(i-1));
            
            Btmp(i-2) = dx(i-1)/(dx(i-1)+dx(i-2));
            Btmp(i-1) = dx(i)/(dx(i)+dx(i-1));
            Btmp_diag_low = dx(i-2)/(dx(i-1)+dx(i-2));
            
            a(i-2) = Atmp(i)*Btmp(i-1)*Atmp(i-1)*Btmp(i-2);
            a(i-1) = -Atmp(i)*Btmp(i-1)*Atmp(i-1) - Atmp(i)^2*Btmp(i-1);
            a(i) = Atmp(i)*Btmp(i-1)*Atmp(i-1)*Btmp_diag_low + Atmp(i)^2;           
            L = L ...
                + sparse([i i i],[i-2 i-1 i],[a(i-2) a(i-1) a(i)],n,n);
                 
            ann = 1;
            L = L ...
                + sparse([n],[n],[ann],n,n);
            
            f = L;
        end
        
        function [L4x,L4y,L4xbig] = d4dx2d(x, y)
            % L4x use Dirichlet boundary condition
            % L4y use Dirichlet boundary condition
            % x is the first dimension coordinate
            % y is the second dimension coordinate
            nx = max(size(x));
            ny = max(size(y));
            
            L4x = operator2.d4dx1d(x);
            L4y = operator2.d4dx1d(y);
            
            L4x_cell = repmat({L4x},1,ny);
            L4xbig = blkdiag(L4x_cell{:});
        end       
    end
end
