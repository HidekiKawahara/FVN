function output = generateSparseFVNset(fs, sigmaT, varargin)
% Generate a set of close to orthogonal unit FVNs for acoustic measurements
% output = generateSparseFVNset(fs, sigmaT)
% output = generateSparseFVNset(fs, sigmaT, seed)
% output = generateSparseFVNset(fs, sigmaT, seed, displayOn)
% output = generateSparseFVNset(fs, sigmaT, seed, displayOn, fileOut)
%
% Argument
%   fs         : sampling frequency (Hz)
%   sigmaT     : duration of the FVN (s)
%   seed       : seed for random number generator
%   displayOn  : display statistics of the generated FVNs (default: 0)
%   fileOut    : file output (1: output, 0: no output: defaut)
%
% Output
%   output   : structure variable with the following field
%                   xFVN: generated sample FVN
%               timeAxis: double 1d vector (sample)
%          frequencyAxis: double 1d vector (sample)
%             mixedPhase: double 1d vector (sample)
%       frequencySegment: double (Hz)
%    normalizedBandwidth: double (Hz)
%                   fftl: double (contents 2^N)
%      samplingFrequency: double (Hz)
%                 fvnSet: generated sparse FVN set (fftl times 16)
%                 sigmaT: double nominal duration (s)
%       totalElapsedTime: (s)

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

startTic = tic;
switch nargin
    case 2
        rng(1234);
        displayOn = 0;
        fileOut = 0;
    case 3
        seedN = varargin{1};
        rng(seedN);
        displayOn = 0;
        fileOut = 0;
    case 4
        seedN = varargin{1};
        rng(seedN);
        displayOn = varargin{2};
        fileOut = 0;
    case 5
        seedN = varargin{1};
        rng(seedN);
        displayOn = varargin{2};
        fileOut = varargin{3};
    otherwise
        help generateSparseFVNset
        output = [];
        return
end
output = generateStandardFVN(fs, sigmaT);
fftl = output.fftl;

n_test = 100;
fvnSamples = zeros(fftl, n_test);
for ii = 1:n_test
    output = generateStandardFVN(fs, sigmaT);
    fvnSamples(:, ii) = output.xFVN;
end

fvnXcorr = zeros(n_test, n_test);
for ii = 1:n_test
    for jj = 1:n_test
        fvnXcorr(ii, jj) = max(abs(xcorr(fvnSamples(:, ii), fvnSamples(:, jj))));
    end
end

% search for the differentiest set

n_max = 16;
n_count = 0;
idxList = zeros(n_max, 1);

minVec = min(fvnXcorr);
[minV, minIdx] = min(minVec);
n_count = n_count + 1;
idxList(n_count) = minIdx;
minvList(n_count) = minV;
while n_count < n_max
    minV = 1;
    for jj = 1:n_count
        %refii = idxList(jj);
        for ii = 1:n_test
            if ~ismember(ii, idxList(1:n_count))
                if minV >= fvnXcorr(jj, ii)
                    minV = fvnXcorr(jj, ii);
                    minIdx = ii;
                end
            end
        end
    end
    n_count = n_count + 1;
    idxList(n_count) = minIdx;
    minvList(n_count) = minV;
end

[~, sortedIdx] = sort(minvList);

% distribition of cross correlation

if displayOn
    gg = fvnXcorr(:);
    gg = gg(gg < 0.9);
    figure;plot(sort(gg), (1:length(gg))/length(gg));grid on;
    hold all;
    stem(sort(minvList), (1:16)/16, 'o-');
    set(gca, 'fontsize', 14);
    xlabel('correlation');
    ylabel('probability');
    title(['fs:' num2str(fs) ' Hz  sigmaT:' num2str(sigmaT)]);
    legend('all (100) FVNs', 'selected (16) FVNs', 'location', 'best');
end

if fileOut
    switch sigmaT
        case 0.01
            fvnMin10ms = fvnSamples(:, idxList(sortedIdx));
            save fvnMin10ms.mat fvnMin10ms fs
        case 0.025
            fvnMin25ms = fvnSamples(:, idxList(sortedIdx));
            save fvnMin25ms.mat fvnMin25ms fs
        case 0.05
            fvnMin50ms = fvnSamples(:, idxList(sortedIdx));
            save fvnMin50ms.mat fvnMin50ms fs
        case 0.1
            fvnMin100ms = fvnSamples(:, idxList(sortedIdx));
            save fvnMin100ms.mat fvnMin100ms fs
        case 0.2
            fvnMin200ms = fvnSamples(:, idxList(sortedIdx));
            save fvnMin200ms.mat fvnMin200ms fs
        case 0.4
            fvnMin400ms = fvnSamples(:, idxList(sortedIdx));
            save fvnMin400ms.mat fvnMin400ms fs
    end
else
    disp('Please use contents of output.');
end
output.fvnSet = fvnSamples(:, idxList(sortedIdx));
output.sigmaT = sigmaT;
output.samplingFrequency = fs;
output.totalElapsedTime = toc(startTic);
end