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
classdef sp_norm
    
    methods(Static)
        
        % function used for renormalizing 1d quantity with quadratic form: spline1d(1:3,1:lsp)
        function f = normalizeSpline1d(spline1d,physicalUnit,fluxUnit)
            
            % spline1d is two dimensional
            a = size(spline1d);
            if (ndims(spline1d)~=2 || a(1)~=3)
                error('Input grid must be two dimensional.');
            end
            
            tmp = spline1d/physicalUnit;
            
            tmp(2,:) = tmp(2,:)*fluxUnit;
            
            tmp(3,:) = tmp(3,:)*fluxUnit*fluxUnit;
            
            f=tmp;
            
        end
        
        % function used for renormalizing 2d quantity with quadratic form: spline2d(1:9,1:lsp,1:lst)
        function f=normalizeSpline2d(spline2d,physicalUnit,fluxUnit)
            
            % spline2d is three dimensional
            a = size(spline2d);
            if (ndims(spline2d)~=3 || a(1)~=9)
                error('Input grid must be three dimensional.');
            end
            
            tmp = spline2d/physicalUnit;
            
            tmp(2:3:end,:) = tmp(2:3:end,:)*fluxUnit;
            
            tmp(3:3:end,:) = tmp(3:3:end,:)*fluxUnit*fluxUnit;
            
            f=tmp;
            
        end
        
    end
    
end