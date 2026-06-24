function ant_multi = ant_merge_channels(varargin)
% ANT_MERGE_CHANNELS  Combine several single-channel ant structs into one.
%
% Each input must already carry .channels(1) (e.g., produced by
% ant_attach_channel). Output is one ant struct whose .channels(k) array
% spans every input channel in order. The top-level legacy fields mirror
% channel 1 so any caller that still reads ant.omega / ant.n / ant.m_modes
% / ant.phi_real etc. keeps working - that caller simply sees channel 1.

if nargin == 0
    error('ant_merge_channels:no_input', 'need at least one ant struct.');
end

ant_multi = struct;
for k = 1:nargin
    a = varargin{k};
    if ~isfield(a, 'channels') || isempty(a.channels)
        error('ant_merge_channels:no_channel', ...
            'input %d has no .channels field; call ant_attach_channel first.', k);
    end
    nc_k = numel(a.channels);
    for j = 1:nc_k
        if k == 1 && j == 1
            ant_multi.channels = a.channels(j);
        else
            ant_multi.channels(end+1) = a.channels(j);
        end
    end
end

% Legacy mirror = channel 1
ch1 = ant_multi.channels(1);
ant_multi.omega       = ch1.omega;
ant_multi.n           = ch1.n;
ant_multi.m_modes     = ch1.m_modes;
ant_multi.num_modes   = ch1.num_modes;
ant_multi.psi         = ch1.psi;
ant_multi.phi_real    = ch1.phi_real;
ant_multi.phi_imag    = ch1.phi_imag;
ant_multi.apara_real  = ch1.apara_real;
ant_multi.apara_imag  = ch1.apara_imag;
ant_multi.dptdp_real  = ch1.dptdp_real;
ant_multi.dptdp_imag  = ch1.dptdp_imag;
ant_multi.dapdp_real  = ch1.dapdp_real;
ant_multi.dapdp_imag  = ch1.dapdp_imag;
ant_multi.num_channels = numel(ant_multi.channels);
end
