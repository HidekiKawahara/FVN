function output = loudnessWithA(x, fs)
% sound pressure level with A-weighting
%   output = loudnessWithA(x, fs)
% Argument
%   x    : input signal. double 1-d vector
%   fs   : sampling frequency (Hz) recommend to use 44100 Hz or higher
%
% Output
%   output   : structure with the following field
%      (note the sound pressure levels are not calibrated)
%        fast: double 1-d vector with sound pressu level in dB
%        slow: double 1-d vector with sound pressu level in dB
%    filtered: double 1-d vector with filtered signal
%        time: double 1-d vector time axis (s)

% copyright 2019 Hideki Kawahara
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

aWeightFir = designAweightingFIRFilter(fs);
n_fir = length(aWeightFir);
y = fftfilt(aWeightFir, [x;zeros(n_fir, 1)]);
y = y(round(n_fir / 2) + (1:length(x)));
r_fast = exp(-1/0.125/fs); % time constant 125 ms
r_slow = exp(-1/fs); % time constant 1 s
y_fast = filter(1, [1 -r_fast] / (1 - r_fast), y.^2);
y_slow = filter(1, [1 -r_slow] / (1 - r_slow), y.^2);
output.fast = 10 * log10(y_fast);
output.slow = 10 * log10(y_slow);
output.filtered = y;
output.time = (0:length(x) - 1)' / fs;
end

function aWeightFir = designAweightingFIRFilter(fs)
%   aWeightFir = designAweightingFIRFilter(fs)

%   by Hideki Kawahara
%   28/May/2014

%
% ----- definition of A-weight -----
%---- JIS C1516:2014
aWeightTable = [...
10 -70.4; ...
12.5 -63.4; ...
16 -56.7; ...
20	-50.5; ...	
25	-44.7	; ...
31.5	-39.4	; ...
40	-34.6	; ...
50	-30.2	; ...
63	-26.2	; ...
80	-22.5	; ...
100	-19.1	; ...
125	-16.1	; ...
160	-13.4	; ...
200	-10.9	; ...
250	-8.6	; ...
315	-6.6	; ...
400	-4.8	; ...
500	-3.2	; ...
630	-1.9; ...
800	-0.8; ...
1000	0; ...
1250	0.6; ...
1600	1; ...
2000	1.2; ...
2500	1.3; ...
3150	1.2; ...
4000	1; ...
5000	0.5; ...
6300	-0.1; ...
8000	-1.1; ...
10000	-2.5; ...
12500	-4.3; ...
16000	-6.6; ...
20000	-9.3];
fxBase = aWeightTable(:,1);
aWeightBase = aWeightTable(:,2);
%
fftl = 4096;
fx = (0:fftl-1)'/fftl*fs;
tmpGain = interp1(fxBase,aWeightBase-10*log10(fxBase/1000),fx,'linear','extrap');
tmpAmp = 10.0.^(tmpGain/20);
tmpAmp(fftl/2+1:fftl) = tmpAmp(fftl/2+1:-1:2);
aWeightFirTmp = fftshift(real(ifft(tmpAmp)));
[~, maxIndex] = max(aWeightFirTmp);
aWeightFir = aWeightFirTmp(maxIndex+(-500:500));
end