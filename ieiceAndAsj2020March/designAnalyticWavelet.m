function output = designAnalyticWavelet(fs, fl, fh, varargin)
% output = designAnalyticWavelet(fs, fl, fh)
% output = designAnalyticWavelet(fs, fl, fh, channels_oct)
% output = designAnalyticWavelet(fs, fl, fh, channels_oct, mag)
% output = designAnalyticWavelet(fs, fl, fh, channels_oct, mag, wintype)
% Generate a set of anaoytic signal wavelets using given envelope the
%   default envelope is the six-term cosine series proosed in ref[1].
%
%
% Arguments
%   fs    : sampling frequency (Hz)
%   fl    : lower bound of the carrier frequency (Hz)
%   fh    : upper bound of the carrier frequency (Hz)
%   channels_oct : number of wavelets in one octave
%                  default 12
%   mag   : envelope stretching factor: mag=1 places zero at 0 and 2*fc
%           default 1
%   wintype : Envelope defining function
%             default six-term cosine series in ref [1]. Available types:
%             'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'
%
% Output
%   output   : structure with the following fields : nch: number of wavelets
%    input_parameters: [1×1 struct] copy of input argument
%             fc_list: [1×nch double] list of carrier frequencies (Hz)
%                wvlt: [1×nch struct] a set of wavelets with fields:
%                  w    : complex column variable, impulse response
%                  bias : -bias:bias gives the index of w
%
% Reference
% [1] Kawahara, H., Sakakibara, K., Morise, M., Banno, H., Toda, T., Irino, T. (2017)
% A New Cosine Series Antialiasing Function and its Application to Aliasing-Free
% Glottal Source Models for Speech and Singing Synthesis. Proc. Interspeech 2017,
% 1358-1362, DOI: 10.21437/Interspeech.2017-15.
% [2] ﻿Nuttall, A. H. (1981). Some windows with very good sidelobe behavior. 
% IEEE Trans. Audio Speech and Signal Processing, 29(1), 84–91.
% 

%
% Copyright 2018 Hideki Kawahara
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

mag = 1;
channels_oct = 12;
wintype = 'sixterm';
switch nargin
    case 3
    case 4
        channels_oct = varargin{1};
    case 5
        channels_oct = varargin{1};
        mag = varargin{2};
    case 6
        channels_oct = varargin{1};
        mag = varargin{2};
        wintype = varargin{3};
    otherwise
        help designAnalyticWavelet
        output = [];
        return;
end
fc_list = fl * 2.0 .^ (0:1/channels_oct:log2(fh/fl));
wvlt = struct;
for ii = 1:length(fc_list)
    fc = fc_list(ii);
    switch wintype
        case 'sixterm'
            [n_to, tx, w] = wv_sixterm(fs, fc, mag);
        case 'hanning'
            [n_to, tx, w] = wv_hanning(fs, fc, mag);
        case 'hamming'
            [n_to, tx, w] = wv_hamming(fs, fc, mag);
        case 'blackman'
            [n_to, tx, w] = wv_blackman(fs, fc, mag);
        case 'nuttall12'
            [n_to, tx, w] = wv_nuttall12(fs, fc, mag);
        case 'kaiser'
            [n_to, tx, w] = wv_kaiser(fs, fc, mag);
        case 'dpss'
            [n_to, tx, w] = wv_dpss(fs, fc, mag);
        otherwise
            help designAnalyticWavelet
            output = [];
            return;
    end
    wvlt(ii).w = w;
    wvlt(ii).bias = n_to;
    wvlt(ii).t_axis = tx;
end

%---
input_parameters.sampling_frequency = fs;
input_parameters.lower_frequency = fl;
input_parameters.higher_frequency = fh;
input_parameters.stretching_factor = mag;
input_parameters.channels_per_octave = channels_oct;
input_parameters.wintype = wintype;
%---
output.input_parameters = input_parameters;
output.fc_list = fc_list;
output.wvlt = wvlt;
end

function [n_to, tx, w] = wv_sixterm(fs, fc, mag)
ak = [0.2624710164;0.4265335164;0.2250165621;0.0726831633;0.0125124215;0.0007833203];
n_to = round(mag * 3 * fs / fc);
tx_we = (-n_to:n_to)'/ (2 * n_to + 1);
tx = (-n_to:n_to)'/ fs;
we = cos(tx_we * [0 1 2 3 4 5] * pi * 2) * ak;
wr = we .* cos(2 * pi * tx * fc);
wi = we .* sin(2 * pi * tx * fc);
w = (wr+1i*wi)/sum(we);
end

function [n_to, tx, w] = wv_hanning(fs, fc, mag)
ak = [0.5;0.5];
n_to = round(mag * fs / fc);
tx_we = (-n_to:n_to)'/ (2 * n_to + 1);
tx = (-n_to:n_to)'/ fs;
we = cos(tx_we * [0 1] * 2 * pi) * ak;
wr = we .* cos(2 * pi * tx * fc);
wi = we .* sin(2 * pi * tx * fc);
w = (wr+1i*wi)/sum(we);
end

function [n_to, tx, w] = wv_hamming(fs, fc, mag)
ak = [0.53836;0.46164];
n_to = round(mag * fs / fc);
tx_we = (-n_to:n_to)'/ (2 * n_to + 1);
tx = (-n_to:n_to)'/ fs;
we = cos(tx_we * [0 1] * 2 * pi) * ak;
wr = we .* cos(2 * pi * tx * fc);
wi = we .* sin(2 * pi * tx * fc);
w = (wr+1i*wi)/sum(we);
end

function [n_to, tx, w] = wv_blackman(fs, fc, mag)
ak = [0.42;0.5;0.08];
n_to = round(mag * 1.5 * fs / fc);
tx_we = (-n_to:n_to)'/ (2 * n_to + 1);
tx = (-n_to:n_to)'/ fs;
we = cos(tx_we * [0 1 2] * 2 * pi) * ak;
wr = we .* cos(2 * pi * tx * fc);
wi = we .* sin(2 * pi * tx * fc);
w = (wr+1i*wi)/sum(we);
end

function [n_to, tx, w] = wv_nuttall12(fs, fc, mag)
ak = [0.355768;0.487396;0.144232;0.012604];
n_to = round(mag * 2 * fs / fc);
tx_we = (-n_to:n_to)'/ (2 * n_to + 1);
tx = (-n_to:n_to)'/ fs;
we = cos(tx_we * [0 1 2 3] * 2 * pi) * ak;
wr = we .* cos(2 * pi * tx * fc);
wi = we .* sin(2 * pi * tx * fc);
w = (wr+1i*wi)/sum(we);
end

function [n_to, tx, w] = wv_kaiser(fs, fc, mag)
n_to = round(mag * 2 * fs / fc);
tx = (-n_to:n_to)'/ fs;
we = kaiser(length(tx), 15.03); % 15.03 has the closest max sidelobe to 6-term
wr = we .* cos(2 * pi * tx * fc);
wi = we .* sin(2 * pi * tx * fc);
w = (wr+1i*wi)/sum(we);
end

function [n_to, tx, w] = wv_dpss(fs, fc, mag)
n_to = round(mag * 2 * fs / fc);
tx = (-n_to:n_to)'/ fs;
[ee, ~] = dpss(length(tx), 4.72); % 4.72 has the closest max sidelobe to 6-term
we = ee(:, 1);
wr = we .* cos(2 * pi * tx * fc);
wi = we .* sin(2 * pi * tx * fc);
w = (wr+1i*wi)/sum(we);
end
