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
function P = read_poincare_particles(path, filename)
% READ_POINCARE_PARTICLES
% Read Poincare data, group by particle, and remove energy outliers
%
% Input:
%   path     - Data path (string)
%   filename - File name (e.g. 'poincare.txt')
%
% Output:
%   P        - Structure array grouped by particle (cleaned)

    % === 1. Read parameters (preserve original logic) ===
    %run(fullfile(path, 'read_para_w.m'));

    % === 2. Concatenate path ===
    path_full = fullfile(path);

    % === 3. Read data ===
    data = readtable(fullfile(path_full, filename));

    % === 4. Convert to struct ===
    Poincare = struct( ...
        'psi', data.psi, ...
        'theta', data.theta, ...
        'R', data.R, ...
        'Z', data.Z, ...
        'q', data.q, ...
        'particle_ID', data.particle_ID, ...
        'time_step', data.time_step, ...
        'u_psi', data.zeta, ...
        'zeta', data.u_theta, ...
        'P_zeta', data.u_zeta, ...
        'E', data.E ...
    );

    % === 5. Group by particle ===
    num_particles = max(Poincare.particle_ID);
    fields = fieldnames(Poincare);

    P = struct();

    for i = 1:num_particles
        idx = (Poincare.particle_ID == i);

        for k = 1:length(fields)
            P(i).(fields{k}) = Poincare.(fields{k})(idx);
        end
    end

    % === 6. Remove bad points (energy anomalies)===
    for i = 1:length(P)

        E = P(i).E;
        bad = isoutlier(E, 'median');
        idx = ~bad;

        fields = fieldnames(P(i));

        for k = 1:length(fields)
            P(i).(fields{k}) = P(i).(fields{k})(idx);
        end
    end

end