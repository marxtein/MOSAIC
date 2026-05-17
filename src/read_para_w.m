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

% % read gtc.out parameters

file1=[path,'gtc.out'];
file2=[path,'gtc.out0'];
if(exist(file1,'file'))
    fid = fopen(file1);
elseif(exist(file2,'file'))
    fid = fopen(file2);
else
    msgbox('read_para.m, No gtc.out', 'Error');
end

s = fgetl(fid);

% jend=200;
% for j=1:jend  % read first jend lines
% while ischar(s) % read all lines
jj=1;
while (jj<600 && ischar(s)) % if 'npe' error, larger '600'
    jj=jj+1;
    
    s=regexprep(s, 'MSTEP * =', 'MSTEP='); s(s==',')='';    
    if(strfind(s,'MSTEP='))
        s(1:7)=[]; mstep=str2num(s);
    end
    s=regexprep(s, 'MSNAP * =', 'MSNAP='); s(s==',')='';    
    if(strfind(s,'MSNAP='))
        s(1:7)=[]; msnap=str2num(s);
    end
    s=regexprep(s, 'NDIAG * =', 'NDIAG='); s(s==',')='';    
    if(strfind(s,'NDIAG='))
        s(1:7)=[]; ndiag=str2num(s);
    end
    s=regexprep(s, 'TSTEP * =', 'TSTEP='); s(s==',')='';    
    if(strfind(s,'TSTEP='))
        s(1:7)=[]; dt0=str2num(s);
    end
    s=regexprep(s, 'MPSI * =', 'MPSI='); s(s==',')='';    
    if(strfind(s,'MPSI='))
        s(1:6)=[]; mpsi=str2num(s);
    end
    s=regexprep(s, 'MTHETAMAX * =', 'MTHETAMAX='); s(s==',')='';    
    if(strfind(s,'MTHETAMAX='))
        s(1:11)=[]; mthetamax=str2num(s);
    end
    s=regexprep(s, 'ELOAD * =', 'ELOAD='); s(s==',')='';    
    if(strfind(s,'ELOAD='))
        s(1:7)=[]; eload=str2num(s);
    end
    
    s=regexprep(s, 'MTOROIDAL * =', 'MTOROIDAL='); s(s==',')='';    
    if(strfind(s,'MTOROIDAL='))
        s(1:11)=[]; mtoroidal=str2num(s);
    end
    s=regexprep(s, 'PSI0 * =', 'PSI0='); s(s==',')='';    
    if(strfind(s,'PSI0='))
        s(1:6)=[]; psi0=str2num(s);
    end
    s=regexprep(s, 'PSI1 * =', 'PSI1='); s(s==',')='';    
    if(strfind(s,'PSI1='))
        s(1:6)=[]; psi1=str2num(s);
    end
    s=regexprep(s, 'NUMBERPE * =', 'NUMBERPE='); s(s==',')='';    
    if(strfind(s,'NUMBERPE=')) % number of PE
        s(1:10)=[]; npe=str2num(s);
    end   
    s=regexprep(s, 'AION * =', 'AION='); s(s==',')='';    
    if(strfind(s,'AION=')) 
        s(1:7)=[]; aion=str2num(s);
    end   
    
    s=regexprep(s, 'R0 * =', 'R0='); s(s==',')='';    
    if(strfind(s,'R0='))
        s(1:4)=[]; R0=str2num(s);
    end
    s=regexprep(s, 'B0 * =', 'B0='); s(s==',')='';    
    if(strfind(s,'B0='))
        s(1:4)=[]; B0=str2num(s);
    end
    s=regexprep(s, 'ETEMP0 * =', 'ETEMP0='); s(s==',')='';    
    if(strfind(s,'ETEMP0='))
        s(1:8)=[]; etemp0=str2num(s);
    end
    s=regexprep(s, 'EDEN0 * =', 'EDEN0='); s(s==',')='';    
    if(strfind(s,'EDEN0='))
        s(1:7)=[]; eden0=str2num(s);
    end
    s=regexprep(s, 'RHO0 * =', 'RHO0='); s(s==',')='';    
    if(strfind(s,'RHO0='))
        s(1:6)=[]; rho0=str2num(s);
    end
    
    s=regexprep(s, ' BETAE * =', 'BETAE='); s(s==',')='';  % read betae  
    if(strfind(s,'BETAE='))
        s(1:7)=[]; betae=str2num(s);
    end

    js=strfind(s,'a_minor='); % read a_minor
    if(js)
        s(s==',')=''; s(1:(js+9))=[]; a_minor=str2num(s);
    end
    
    
    if(strfind(s,'nmodes='))
        s(1:8)=[]; nmodes=str2num(s);
    end
    if(strfind(s,'mmodes='))
        s(1:8)=[]; mmodes=str2num(s);
    end
    if(strfind(s,'qiflux='))
        s(1:8)=[]; qiflux=str2num(s);
    end
    if(strfind(s,'rgiflux='))
        s(1:9)=[]; rgiflux=str2num(s);
    end
    if(strfind(s,'iflux='))
        s(1:7)=[]; iflux=str2num(s);
    end

    s = fgetl(fid);
end





fclose(fid);

% % read gtc.in parameters
% fid = fopen([path,'gtc.in']);
% 
% s1 = fgetl(fid);
% while ischar(s1)
% % for j=1:10
%     s1 = fgetl(fid);
%     s=strtok(s1);
%     if(strfind(s,'mstep='))
%         s(1:6)=[]; mstep=str2num(s);
%     end
%     if(strfind(s,'msnap='))
%         s(1:6)=[]; msnap=str2num(s);
%     end
%     if(strfind(s,'ndiag='))
%         s(1:6)=[]; ndiag=str2num(s);
%     end
%     if(strfind(s,'tstep='))
%         s(1:6)=[]; tstep=str2num(s);
%     end
%     if(strfind(s,'mpsi='))
%         s(1:5)=[]; mpsi=str2num(s);
%     end
%     if(strfind(s,'mthetamax='))
%         s(1:10)=[]; mthetamax=str2num(s);
%     end
%     if(strfind(s,'mtoroidal='))
%         s(1:10)=[]; mtoroidal=str2num(s);
%     end
%     if(strfind(s,'psi0='))
%         s(1:5)=[]; psi0=str2num(s);
%     end
%     if(strfind(s,'psi1='))
%         s(1:5)=[]; psi1=str2num(s);
%     end
%     if(strfind(s,'r0='))
%         s(1:3)=[]; r0=str2num(s);
%     end
%     
% %     disp(s);
% %     strfind(s1,'mstep=');
% end
% 
% fclose(fid);


%%
% s=' nmodes=   4   9  14  18  22  27  31  36'; 
% if(strfind(s,'nmodes='))
%     s(1:8)=[]; nmodes=str2num(s);
% end
