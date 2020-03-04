function output = waveletAttributesAnalyzer(x, fs, wvltStr)
% output = waveletAttributesAnalyzer(x, fs, wvltStr)
%
% Signal attribute analysis using wavelet filters based on an analyitic signal
% 
% Arguments
%   x       : input signal. One column vector
%   fs      : samplinlg frequency (Hz)
%   wvltStr : structure variable consisting of the analytic signals
%             Fileds used:
%     fc_list  : a vector of carrier frequencies (Hz)
%     wvlt     : a set of structure variables. 
%                Each variable defines an analytic signal. Two fields
%       w      : complex column vector defining an analytic signal
%       bias   : half length of the vector (2 * bias + 1) is the length
%
% Output
%   output : structure variabl with the following fields
%     rawWavelet                : Each column has filtered output
%     inst_freq_map             : Each column has sample-wise inst. freq.
%     amp_squared_map           : Each column has product of the absolute
%                                 values of succeeding samples
%     group_delay_map           : Each column has sample-wise group delay
%     elapsedTimeForFiltering   : Elapsed time for filtering (s)
%     elapsedTimeForPostProcess : Elapsed time for postprocess (s)
%     elapsedTime               : Total elapsed time (s)

% Designed and implemented by Hideki Kawahara
%
%Copyright 2018 Hideki Kawahara
%
%Licensed under the Apache License, Version 2.0 (the "License");
%you may not use this file except in compliance with the License.
%You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
%Unless required by applicable law or agreed to in writing, software
%distributed under the License is distributed on an "AS IS" BASIS,
%WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%See the License for the specific language governing permissions and
%limitations under the License.


start_tic = tic;
n_channels = length(wvltStr.fc_list);
n_data = length(x);
rawWavelet = zeros(n_data, n_channels);
n_buffer = n_data;
x = [zeros(wvltStr.wvlt(1).bias, 1);x;zeros(wvltStr.wvlt(1).bias, 1)];

tic
buffer_index = 1:n_data;
for ii = 1:n_channels
    y = fftfilt(wvltStr.wvlt(ii).w, x);
    rawWavelet(:, ii) = y(wvltStr.wvlt(ii).bias + wvltStr.wvlt(1).bias + buffer_index);
end
elapsedTimeForFiltering = toc;
%% 
tic
fc_list = wvltStr.fc_list;
amp_map = abs(rawWavelet);
amp_squared_map = amp_map([2:n_buffer n_buffer], :) .* amp_map;
inst_freq_map = angle(rawWavelet([2:n_buffer n_buffer], :) ./ rawWavelet) * fs / 2 / pi;
inst_freq_map(end, :) = inst_freq_map(end - 1, :);
fc_list_extend = [fc_list fc_list(end) * fc_list(end) / fc_list(end - 1)];
freq_step = diff(fc_list_extend);
group_delay_map = -angle(rawWavelet(:, [2:end end]) ./ rawWavelet) * diag(1.0 ./ freq_step) / 2 / pi;
group_delay_map(:, end) = group_delay_map(:, end - 1);
elapsedTimeForPostProcess = toc;
output.rawWavelet = rawWavelet;
output.amp_squared_map = amp_squared_map;
output.inst_freq_map = inst_freq_map;
output.group_delay_map = group_delay_map;
output.elapsedTimeForFiltering = elapsedTimeForFiltering;
output.elapsedTimeForPostProcess = elapsedTimeForPostProcess;
output.elapsedTime = toc(start_tic);
end