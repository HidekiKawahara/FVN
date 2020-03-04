function output = sourceAttributesAnalysis(x, fs, varargin)
% This function analysis source information using an analytic signal
% with the six-term cosine series envelope
%
% sourceAttributesAnalysis(x, fs)
% sourceAttributesAnalysis(x, fs, range)
% sourceAttributesAnalysis(x, fs, range, outputStruct) % After initialization
% sourceAttributesAnalysis(x, fs, range, low_frequency, high_freuency)
% sourceAttributesAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave);
% sourceAttributesAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave, dsOn)
% sourceAttributesAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave, dsOn, stretching_factor)
% sourceAttributesAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave, dsOn, stretching_factor, wintype)
% sourceAttributesAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave, dsOn, stretching_factor, wintype, ...
%     integration_time)
% sourceAttributesAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave, dsOn,  stretching_factor, wintype, ...
%     integration_time, sampling_multiplier)
%
% Arguments
%   x       : input signal. One column vector
%   fs      : sampling frequency (Hz)
%   range   : index of the start and end point of the input
%   outputStruct : output of the structure variable generated in the first
%                  run. Using this skips initialization
%   low_frequency : Lower limit of periodicity check (Hz)
%   high_freuency : Higher limit of periodicity check (Hz)
%   channels_in_octave  : Number of filters in each octave
%   dsOn : Switch for internal downsampling. 1:downsample, 0:no downsampling
%   stretching_factor : Envelope stretching factor: 1: isometric in the
%                       time-frequency plane. Larger the longer in time
%   wintype           : string variable representing envelope functionu
%     'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'
%   integration_time  : integration time for periodicity measure (ms)
%   sampling_multiplier : coefficient to determine sampling rate for
%                         downsampling
%
% Output
%   output  : structure variable with the following field
%     estPeriod   : indicator of periodicity 1:periodic 0:random (sd hoc!)
%     wvltStrDs   : structure consisting of used analytic signal bank
%     rawWavelet  : Each column vector consits of each filter output
%     downSamplinbgRate : A mnumber representing downsampling ratio
%     fftlds   : FFT buffer length for debug
%     half_aaf_length : for debug
%     w_aaf           : for debug
%     time_axis_wavelet : time axis for wavelet analysis results
%     signal_time_axis  : time axis for input signal
%     gd_dev_map        : output deviation based on group delay
%     dgd_dev_map       : output deviation based on diffed group delay
%     mix_rev_measure   : mixed periodicity meeasure based on gd and dgd
%     fixed_points_freq : frequency of fo candidates
%     fixed_points_measure : mixed measure of fo candidates
%     fixed_points_amp  : signal amplitude of fo condidates
%     elapsedTime       : total elapsed time (s)

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
output = struct;
output.narg = nargin;
fftl = 8192;
low_frequency = 55;
high_freuency = 1200;
stretching_factor = 1.0;
channels_in_octave = 24;
isInitialized = 0;
dsOn = 1;
wintype = 'sixterm';
integration_time = 40; % ms
sampling_multiplier = 6;
switch nargin
    case {5, 6, 7, 8, 9, 10, 11}
        ixrange = varargin{1};
        if length(ixrange) ~= 2
            help sourceInformationAnalysis
            disp('Range is the start and the end index of the data');
        end
        low_frequency = varargin{2};
        high_freuency = varargin{3};
    case 4
        ixrange = varargin{1};
        if length(ixrange) ~= 2
            help sourceInformationAnalysis
            disp('Range is the start and the end index of the data');
        end
        outputStruct = varargin{2};
        if ~isstruct(outputStruct)
            help sourceInformationAnalysis
            disp('The last argument should be an initialized wavelet structure');
        else
            isInitialized = 1;
        end
end
switch nargin
    case 2
        ixrange = zeros(2, 1);
        ixrange(1) = 1;
        ixrange(2) = length(x);
    case 3
        ixrange = varargin{1};
        if length(ixrange) ~= 2
            help sourceInformationAnalysis
            disp('Range is the start and the end index of the data');
        end
    case 4
    case 5
    case 6
        channels_in_octave = varargin{4};
    case 7
        channels_in_octave = varargin{4};
        dsOn = varargin{5};
    case 8
        channels_in_octave = varargin{4};
        dsOn = varargin{5};
        stretching_factor = varargin{6};
    case 9
        channels_in_octave = varargin{4};
        dsOn = varargin{5};
        stretching_factor = varargin{6};
        wintype = varargin{7};
    case 10
        channels_in_octave = varargin{4};
        dsOn = varargin{5};
        stretching_factor = varargin{6};
        wintype = varargin{7};
        integration_time = varargin{8};
    case 11
        channels_in_octave = varargin{4};
        dsOn = varargin{5};
        stretching_factor = varargin{6};
        wintype = varargin{7};
        integration_time = varargin{8};
        sampling_multiplier = varargin{9};
    otherwise
        help sourceInformationAnalysis
end
ixrange = round(ixrange);
x_trim = x(ixrange(1):ixrange(2));
if ~isInitialized
    downSamplinbgRate = 1;
    if dsOn
        downSamplinbgRate = max(1, floor (fs / (high_freuency * sampling_multiplier)));
    end
    fsd = fs / downSamplinbgRate;
    half_aaf_length = round(max(1, downSamplinbgRate * 3.5 - 1));
    w_aaf = blackman(round(2 * half_aaf_length + 1));
    fftlds = fftl / 2;
    wvltStrDs = designAnalyticWavelet(fsd, low_frequency, high_freuency, ...
        channels_in_octave, stretching_factor, wintype);
else
    downSamplinbgRate = outputStruct.downSamplinbgRate;
    wvltStrDs = outputStruct.wvltStrDs;
    fftlds = outputStruct.fftlds;
    half_aaf_length = outputStruct.half_aaf_length;
    w_aaf = outputStruct.w_aaf;
    fsd = fs / downSamplinbgRate;
    integration_time = outputStruct.integration_time;
    sampling_multiplier = outputStruct.sampling_multiplier;
    low_frequency = outputStruct.wvltStrDs.input_parameters.lower_frequency;
    high_freuency = outputStruct.wvltStrDs.input_parameters.higher_frequency;
end
if downSamplinbgRate > 1
xd = fftfilt(w_aaf/sum(w_aaf), x_trim);
xd = xd(half_aaf_length:downSamplinbgRate:end);
else
    xd = x;
end

n_channels = length(wvltStrDs.fc_list);
fc_list = wvltStrDs.fc_list;

outputDs = waveletAttributesAnalyzer(xd, fsd, wvltStrDs);
%%

n_samples = length(outputDs.rawWavelet(:, 1));
sample_distance = round(integration_time / 1000 * fsd);
selector = sample_distance:n_samples - sample_distance;
bias = floor(sample_distance / 2);
amp_sq = outputDs.amp_squared_map;
wweight = hanning(sample_distance);
%wweight = ones(sample_distance, 1);
amp_sm_tmp = fftfilt(wweight, amp_sq .* amp_sq(:, [1 1:end-1]));
gd = abs(outputDs.group_delay_map);
gd_sm_tmp = fftfilt(wweight, amp_sq .* amp_sq(:, [1 1:end-1]) .* gd .* gd(:, [1 1:end-1]));
amp_sm_tmp(amp_sm_tmp <= 0) = min(amp_sm_tmp(amp_sm_tmp > 0));
gd_sm = gd_sm_tmp ./ amp_sm_tmp;
gd_dev_map = gd_sm(selector - bias, :) * diag(fc_list .^ 2);

amp_sm_tmp_if = fftfilt(wweight, amp_sq .* amp_sq([1 1:end-1], :));
amp_sm_tmp_if(amp_sm_tmp_if <= 0) = min(amp_sm_tmp_if(amp_sm_tmp_if > 0));
dif_map = abs(outputDs.inst_freq_map([1 1:end-1], :) - outputDs.inst_freq_map);
if_sm_tmp = fftfilt(wweight, amp_sq .* amp_sq(:, [1 1:end-1]) .* dif_map .* dif_map(:, [1 1:end-1]));
if_sm_map = if_sm_tmp ./ amp_sm_tmp_if;
if_dev_map = if_sm_map(selector - bias, :) * diag(1 ./ fc_list .^ 4) * fsd ^2;
mix_measure = sqrt((gd_dev_map + if_dev_map) / 2);
mix_measure(mix_measure <= 0) = min(mix_measure(mix_measure > 0));
mix_rev_measure = 1.0 ./ mix_measure;

if_smoothed = fftfilt(wweight, outputDs.inst_freq_map .* outputDs.amp_squared_map);
amp_smoothed = fftfilt(wweight, outputDs.amp_squared_map);
if_smoothed = if_smoothed ./ amp_smoothed;
if_smoothed = if_smoothed(selector - bias, :);

if_norm_map = if_smoothed * diag(1 ./ fc_list);
log_if_dev_map = log(max(0.5, if_norm_map));
fixed_points = zeros(n_samples, n_channels);
fixed_points_freq = zeros(n_samples, n_channels);
fixed_points_measure = zeros(n_samples, n_channels);
fixed_points_amp = zeros(n_samples, n_channels);
amplitude = abs(outputDs.rawWavelet);
channel_idx = 1:n_channels;
for ii = selector
    buffer_id = ii - selector(1) + 1;
    tmp_fixp = channel_idx(log_if_dev_map(buffer_id, :) .* log_if_dev_map(buffer_id, [2:end, end]) < 0 & ...
        log_if_dev_map(buffer_id, :) > 0 & if_smoothed(buffer_id, :) > low_frequency & ...
        if_smoothed(buffer_id, :) < high_freuency);
    if ~isempty(tmp_fixp)
        fixed_points(buffer_id, 1:length(tmp_fixp)) = tmp_fixp(:)';
        for kk = 1:length(tmp_fixp)
            r = log_if_dev_map(buffer_id, tmp_fixp(kk)) / ...
                (log_if_dev_map(buffer_id, tmp_fixp(kk)) - log_if_dev_map(buffer_id, tmp_fixp(kk) + 1));
            fixed_points_freq(buffer_id, kk) = ...
                fc_list(tmp_fixp(kk)) + (fc_list(tmp_fixp(kk) + 1) - fc_list(tmp_fixp(kk))) * r;
            fixed_points_measure(buffer_id, kk) = (1 - r) * mix_rev_measure(buffer_id, max(1, tmp_fixp(kk) - 1)) ...
                + r * mix_rev_measure(buffer_id, tmp_fixp(kk));
            fixed_points_amp(buffer_id, kk) = (1 - r) * amplitude(buffer_id, max(1, tmp_fixp(kk) - 1)) ...
                + r * amplitude(buffer_id, tmp_fixp(kk));
        end
        [sortedv, sortIdx] = sort(fixed_points_measure(buffer_id, 1:length(tmp_fixp)), 'descend');
        fixed_points_measure(buffer_id, 1:length(tmp_fixp)) = sortedv;
        fixed_points_freq(buffer_id, 1:length(tmp_fixp)) = fixed_points_freq(buffer_id, sortIdx);
        fixed_points_amp(buffer_id, 1:length(tmp_fixp)) = fixed_points_amp(buffer_id, sortIdx);
    end
end
estSNR = revMeasure2SNR(fixed_points_measure, wintype);
output.estPeriod = 1 ./ (1 + exp(-(estSNR - 15) / 6));
output.sampling_frequency = fs;
output.wvltStrDs = wvltStrDs;
output.rawWavelet = outputDs.rawWavelet;
output.inst_freq_map = outputDs.inst_freq_map;
output.if_smooth_map = if_smoothed;
output.downSamplinbgRate = downSamplinbgRate;
output.fftlds = fftlds;
output.half_aaf_length = half_aaf_length;
output.w_aaf = w_aaf;
output.integration_time = integration_time;
output.sampling_multiplier = sampling_multiplier;
output.time_axis_wavelet = (1:n_samples) / fsd + ixrange(1) / fs;
output.signal_time_axis = (ixrange(1):ixrange(2)) / fs;
output.gd_dev_map = gd_dev_map;
output.if_dev_map = if_dev_map;
output.amp_smoothed = amp_smoothed(selector - bias, :);
output.n_effective = size(if_dev_map, 1);
output.mix_rev_measure = mix_rev_measure;
output.mix_measure = mix_measure;
output.fixed_points_freq = fixed_points_freq;
output.fixed_points_measure = estSNR;
output.fixed_points_amp = fixed_points_amp;
output.elapsedTime = toc(start_tic);
end

function dBsnr = revMeasure2SNR(revmeasure, wintype)
%winthpe_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
%caliblation_list =  [-14.6137 -10.3863 6.8289 -14.4067 -14.6077 -14.6057 -14.6175];
tmpdB = 20 * log10(abs(revmeasure));
switch wintype
    case 'sixterm'
        dBsnr = tmpdB - 14.6137;
    case 'hanning'
        dBsnr = tmpdB - 10.3863;
    case 'hamming'
        dBsnr = tmpdB + 6.8289;
    case 'blackman'
        dBsnr = tmpdB - 14.4067;
    case 'nuttall12'
        dBsnr = tmpdB - 14.6077;
    case 'kaiser'
        dBsnr = tmpdB - 14.6057;
    case 'dpss'
        dBsnr = tmpdB - 14.6175;
end
end