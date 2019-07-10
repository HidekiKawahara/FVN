function output = oneThirdSpecDisplayForFVN(x, fs, spl, varargin)
% display A-weighted sound pressure level and spectral analysis
%   output = oneThirdSpecDisplayForFVN(x, fs, spl)
%   output = oneThirdSpecDisplayForFVN(x, fs, spl, displayOn)
%   output = oneThirdSpecDisplayForFVN(x, fs, spl, displayOn, msgText)
% Argument
%   x          : double 1-d vector, signal
%   fs         : sampling frequency (Hz)
%   spl        : sound pressure level of the signal with A-weight  (dB)
%   displayOn  : 1:draw figures, 0: no figure
%   msgText    : text string for the figure title
%
% Output
%   output : structure consisting of the follwing fields
%                 sigWavePower: double, raw signal power
%                    time_axis: double 1-d vector
%                    fastSPLdB: double 1-d vector
%                    slowSPLdB: double 1-d vector
%    rawLevelAnalysisStructure: [1~1 struct] calibration information
%          signalPowerSpectrum: double 1-d vector
%      backgroundPowerSpectrum: double 1-d vector
%       frequencyAxisPowerSpec: double 1-d vector
%                         fftl: 1048576
%            samplingFrequency: 44100
%             signalOneThirdDB: double 1-d vector
%         backgroundOneThirdDB: double 1-d vector
%        frequencyAxisOneThird: double 1-d vector
%            calibrationFactor: raw signal level to SPL (dB)

% Copyright 2019 Hideki Kawahara
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

switch nargin
    case 3
        displayOn = 1;
        msgText = ' ';
    case 4
        displayOn = varargin{1};
        msgText = ' ';
    case 5
        displayOn = varargin{1};
        msgText = varargin{2};        
    otherwise
        help oneThirdSpecDisplayForFVN
        output = [];
        return;
end
%--- calibrate level
levelStr = loudnessWithA(x, fs);
t_axis = levelStr.time;
bgLevel = mean(levelStr.fast(t_axis > t_axis(end) * 0.6));
th_level = (max(levelStr.fast) + bgLevel) / 2;
signalOn = levelStr.fast > th_level;
onSetLoc = min(t_axis(signalOn > signalOn([1 1:end - 1])));
offSetLoc = max(t_axis(signalOn < signalOn([1 1:end - 1])));
sigLevel = mean(levelStr.fast(t_axis > onSetLoc + 1 & t_axis < offSetLoc - 1));
calibrationFactor = spl - sigLevel;
%---- time course
if displayOn
    figure;
    plot(t_axis, levelStr.fast + calibrationFactor, t_axis, levelStr.slow + calibrationFactor);
    grid on;
    axis([t_axis(1) t_axis(end) [bgLevel + calibrationFactor - 10, spl + 10]]);
    xlabel('time (s)');
    ylabel('sound pressure level (dB, A-weight)');
    legend('fast', 'slow', 'location', 'best');
    title(msgText);
end
%----- linear spectrum
xspl = x * 10 ^ (calibrationFactor / 20);
xOn = xspl(t_axis > onSetLoc + 0.3 & t_axis < offSetLoc - 0.3);
n_segment = length(xOn);
xOff = xspl(end - n_segment + 1:end);
fftl = 2 ^ ceil(log2(n_segment) + 1);
fx = (0:fftl - 1)' / fftl * fs;
w = blackman(n_segment);
w = w / sqrt(sum(w .^ 2));
pwOn = abs(fft(xOn .* w, fftl)) .^ 2 / fftl;
pwOff = abs(fft(xOff .* w, fftl)) .^ 2 / fftl;
if displayOn
    figure;
    semilogx(fx, 10 * log10(pwOn), fx, 10 * log10(pwOff));
    grid on;
    maxLvl = 10 * log10(max(max(pwOn), max(pwOff)));
    set(gca, 'xlim', [10 fs / 2], 'ylim', maxLvl + [-100 5]);
    xlabel('frequency (Hz)');
    ylabel('sound pressure level (dB, for bin)');
    legend('signal', 'background', 'location', 'best');
    title(msgText);
end
%---- one third octave
fxH = fx(1:fftl / 2) * 2 ^ (1/6);
fxL = fx(1:fftl / 2) * 2 ^ (-1/6);
cumPwOn = cumsum(pwOn) * 2;
cumPwOff = cumsum(pwOff) * 2;
pwOnH = interp1(fx, cumPwOn, fxH, 'linear', 'extrap');
pwOnL = interp1(fx, cumPwOn, fxL, 'linear', 'extrap');
pwOffH = interp1(fx, cumPwOff, fxH, 'linear', 'extrap');
pwOffL = interp1(fx, cumPwOff, fxL, 'linear', 'extrap');
onethOn = pwOnH - pwOnL;
onethOff = pwOffH - pwOffL;
if displayOn
    figure;
    semilogx(fx(1:fftl / 2), 10 * log10(onethOn), fx(1:fftl / 2), 10 * log10(onethOff));
    grid on;
    set(gca, 'xlim', [10 fs / 2]);
    xlabel('frequency (Hz)');
    ylabel('sound pressure level (dB, 1/3)');
    legend('signal', 'background', 'location', 'best');
    title(msgText);
end
%---- output copy
output.sigWavePower = sum((xOn .* w) .^ 2);
output.time_axis = t_axis;
output.fastSPLdB = levelStr.fast + calibrationFactor;
output.slowSPLdB = levelStr.slow + calibrationFactor;
output.rawLevelAnalysisStructure = levelStr;
output.signalPowerSpectrum = pwOn;
output.backgroundPowerSpectrum = pwOff;
output.frequencyAxisPowerSpec = fx;
output.fftl = fftl;
output.samplingFrequency = fs;
output.signalOneThirdDB = 10 * log10(onethOn);
output.backgroundOneThirdDB = 10 * log10(onethOff);
output.frequencyAxisOneThird = fx(1:fftl / 2);
output.calibrationFactor = calibrationFactor;
end
