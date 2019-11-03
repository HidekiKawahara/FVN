function output = allInOneFVNmeasurement(memoText, varargin)
% Acoustic impulse response and nonlinear component measurement
% How to call
% output = allInOneFVNmeasurement(memoText)
% output = allInOneFVNmeasurement(memoText, number_of_FVNs)
% output = allInOneFVNmeasurement(memoText, number_of_FVNs, response_length)
% output = allInOneFVNmeasurement(memoText, number_of_FVNs, response_length, spl)
% output = allInOneFVNmeasurement(memoText, number_of_FVNs, response_length, spl, exmode)
%
% Argument
%  memoText        : text string
%  number_of_FVNs  : number of FVN sequences for nonlinearity measurement
%                    must be equal or larger than 4.
%                    The default value is 4.
%  response_length : observation length of the impulse response (s)
%                    The default value is 0.3 (s)
%  spl             : sound pressure level at microphone in A-weighting (dB)
%                    The default value is 80 dB
%  exmode          : execution mode selector
%                    'measurement', 'measurement2ch', 'calibration', 'diagnosis'
%                    default 'measurement'
%
% Output
%   output         : structure with the following fields
%                    Field name explains itself
%         averagedThirdOctaveResponse: 1-d vector double
%    averagedThirdOctaveNonlinearComp: 1-d vector double
%       averagedThirdOctaveBackground: 1-d vector double
%            averagedResponseWaveform: 1-d vector double
%          individualResponseWaveform: matrix
%         averagedBackgroundWavelform: 1-d vector double
%        individualBackgroundWaveform: matrix
%                            memoText: text string
%                      number_of_FVNs: integer
%                      responseLength: double (s)
%                  soundPressureLevel: double (dB)
%                         signalRange: indecies of start and end
%                     backgroundRange: indecies of start and end
%                   samplingFrequency: double (Hz)
%                     modeOfOperation: text string
%                        creationDate: text string
%                      levelStructure: structure with calibration info.
%                         elapsedTime: double (s)
%                    recordedWaveform: 1-d vector double
%
%  This function displays a figure showing the average responses
%  It also save three files in 'measurement' mode
%  They are EPSF file of the figure, 'mat' file consisting of the output
%  and the WAVE file consisting of the recorded waveform.
%  Note that the 'mat' file does not consists of the recorded waveform
%  field. The files have the same unique file name other than the extension.

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

mypath = mfilename('fullpath');
myDirectory = mypath(1:strfind(mypath, 'allInOneFVNmeasurement')-1);
switch nargin
    case 1
        n_mix = 4;
        to = 0.3;
        spl = 80;
        switch memoText
            case 'calibration'
                exmode = 'calibration';
            otherwise
                exmode = 'measurement';
        end
    case 2
        n_mix = varargin{1};
        switch n_mix
            case {2,3,4,5,6,7,8,9,10}
            otherwise
                disp('The number of FVNs must be from 2 to 10.');
                output.errorMessage = 'The number of FVNs must be from 2 to 10.';
                return;
        end
        switch memoText
            case 'calibration'
                exmode = 'calibration';
            otherwise
                exmode = 'measurement';
        end
        to = 0.3;
        spl = 80;
    case 3
        n_mix = varargin{1};
        switch n_mix
            case {2,3,4,5,6,7,8,9,10}
            otherwise
                disp('The number of FVNs must be from 2 to 10.');
                output.errorMessage = 'The number of FVNs must be from 2 to 10.';
                return;
        end
        to = varargin{2};
        spl = 80;
        switch memoText
            case 'calibration'
                exmode = 'calibration';
            otherwise
                exmode = 'measurement';
        end
    case 4
        n_mix = varargin{1};
        to = varargin{2};
        spl = varargin{3};
        switch memoText
            case 'calibration'
                exmode = 'calibration';
            otherwise
                exmode = 'measurement';
        end
    case 5
        n_mix = varargin{1};
        to = varargin{2};
        spl = varargin{3};
        exmode = varargin{4};
    otherwise
        help allInOneFVNmeasurement
        output.errorMessage = 'Please check how to use this.';
        return;
end

%load basicFVNparamsFor44k.mat a signMat signPeriod
tmpCommand = ['load ' myDirectory 'basicFVNparamsFor44k.mat a signMat signPeriod'];
eval(tmpCommand);
fvn_mat_files = {'fvnMin10ms','fvnMin25ms','fvnMin50ms','fvnMin100ms','fvnMin200ms','fvnMin400ms'};
%durationList = [10, 25, 50, 100, 200, 400] / 1000;
if to < 0.06
    %load fvnMin10ms fvnMin10ms fs
    tmpCommand = ['load ' myDirectory 'fvnMin10ms fvnMin10ms fs'];
    eval(tmpCommand);
    fvnID = 1;
    fvnSet = fvnMin10ms;
elseif to < 0.09
    %load fvnMin25ms fvnMin25ms fs
    tmpCommand = ['load ' myDirectory 'fvnMin25ms fvnMin25ms fs'];
    eval(tmpCommand);
    fvnID = 2;
    fvnSet = fvnMin25ms;
elseif to < 0.177
    %load fvnMin50ms fvnMin50ms fs
    tmpCommand = ['load ' myDirectory 'fvnMin50ms fvnMin50ms fs'];
    eval(tmpCommand);
    fvnID = 3;
    fvnSet = fvnMin50ms;
elseif to < 0.3
    %load fvnMin100ms fvnMin100ms fs
    tmpCommand = ['load ' myDirectory 'fvnMin100ms fvnMin100ms fs'];
    eval(tmpCommand);
    fvnID = 4;
    fvnSet = fvnMin100ms;
elseif to < 0.58
    %load fvnMin200ms fvnMin200ms fs
    tmpCommand = ['load ' myDirectory 'fvnMin200ms fvnMin200ms fs'];
    eval(tmpCommand);
    fvnID = 5;
    fvnSet = fvnMin200ms;
else
    %load fvnMin400ms fvnMin400ms fs
    tmpCommand = ['load ' myDirectory 'fvnMin400ms fvnMin400ms fs'];
    eval(tmpCommand);
    fvnID = 6;
    fvnSet = fvnMin400ms;
end
output.fvnFile = fvn_mat_files{fvnID};
fftlFVN = size(fvnSet, 1);

start_tic = tic;
nto = round(to * fs);
fftl = 32768;
%
buffer = zeros(2 * fftlFVN + 2 * nto * 8 * signPeriod(n_mix), 1);
base_idx = 1:fftlFVN;
for ii = 1:n_mix
    nHead = 0;
    xFVN = fvnSet(:, ii);%fvnStr(ii).xFVN;
    for jj = 1: signPeriod(n_mix) * 8
        buffer(nHead + base_idx) = buffer(nHead + base_idx) + xFVN * signMat(jj, ii);
        nHead = nHead + nto;
    end
    nEnd = nHead + nto + base_idx(end);
end
pinkBuffer = filter(1, a, buffer);
%dataLength = length(pinkBuffer) / fs;
switch exmode
    case 'calibration'
        endPoint = fftlFVN / 2 + nto * 8 * signPeriod(n_mix);
        selector = fftlFVN / 4:endPoint;
        disp(['Test signal plays ' num2str(length(selector) / fs + 2, '%5.0f') ' seconds.']);
        disp(['Adjust the test signal level to ' num2str(spl) ' dB at the microphone']);
        player = audioplayer([pinkBuffer(selector) pinkBuffer(selector) * 0] / max(abs(pinkBuffer)) * 0.8, fs);
        %playblocking(player);
        recObj = audiorecorder(fs, 24, 1);
        record(recObj);
        pause(0.5);
        play(player);
        while strcmp(get(player, 'runnin'), 'on')
            pause(0.5);
        end
        pause(1);
        xrecord = getaudiodata(recObj);
        stop(recObj);
        if max(abs(xrecord)) > 0.9
            disp(['Reduce the mic sensitivity. MaxValue:' num2str(max(abs(xrecord)))]);
        elseif max(abs(xrecord)) < 0.1
            disp(['Increase the mic sensitivity. MaxValue:' num2str(max(abs(xrecord)))]);
        else
            disp(['OK! MaxValue:' num2str(max(abs(xrecord)))])
        end
        output.buffer = buffer;
        output.pinkBuffer = pinkBuffer;
        return;
    case 'measurement'
        endPoint = nEnd;%fftlFVN / 2 + nto * 8 * signPeriod(n_mix) + fs;
        %disp(['Measurement starts. Be quiet for ' num2str((length(buffer) * 1.1) / fs, '%5.0f') ' seconds. Please.']);
        disp(['Measurement starts. Be quiet for ' num2str(endPoint / fs, '%5.0f') ' seconds. Please.']);
        %player = audioplayer([pinkBuffer pinkBuffer * 0] / max(abs(pinkBuffer)) * 0.8, fs);
        player = audioplayer(pinkBuffer(1:nEnd) / max(abs(pinkBuffer)) * 0.8, fs);
        recObj = audiorecorder(fs, 24, 1);
        record(recObj);
        pause(1);
        play(player);
        while strcmp(get(player, 'runnin'), 'on')
            pause(0.5);
        end
        pause(2);
        xrecord = getaudiodata(recObj);
        stop(recObj);
    case 'measurement2ch'
        disp(['Measurement starts. Be quiet for ' num2str((length(buffer) * 1.1) / fs, '%5.0f') ' seconds. Please.']);
        player = audioplayer([pinkBuffer pinkBuffer * 0] / max(abs(pinkBuffer)) * 0.8, fs);
        recObj = audiorecorder(fs, 24, 2);
        record(recObj);
        pause(1);
        play(player);
        while strcmp(get(player, 'runnin'), 'on')
            pause(0.5);
        end
        pause(1);
        xrecord = getaudiodata(recObj);
        stop(recObj);
    case 'diagnosis'
        disp('This mode is not ready yet.');
        help allInOneFVNtest
        output.errorMessage = 'This mode is not ready yet.';
        return;
end

levelStr = oneThirdSpecDisplayForFVN(xrecord(:, 1), fs, spl, 0);

xr = filter(a, 1, xrecord);

fx =(0:fftl - 1)' / fftl * fs;
tmpResp = fftfilt(fvnSet(end:-1:1, 1), xr);
switch exmode
    case 'measurement'
        gg = fftfilt(hanning(101),tmpResp .^2);
    case 'measurement2ch'
        gg = fftfilt(hanning(101),tmpResp(:, 2) .^2);
end
maxgg = max(gg);
ggtrim = gg .* (gg > maxgg * 0.7);
idx = 1:length(gg);
peaksd = idx(ggtrim > ggtrim([1, 1:end-1]) & ggtrim >= ggtrim([2:end, end]));
peaksd = peaksd - 50;
if ~ismember(length(peaksd), signPeriod * 8)
    disp(['the size is wrong:' num2str(length(peaksd))])
    output.xrecord = xrecord;
    output.equalizedRecord = xr;
    output.memoText = memoText;
    output.n_mix = n_mix;
    output.responseLength = to;
    output.errorMessage = ['the size is wrong:' num2str(length(peaksd))];
    return
end
n_mix = round(log2(length(peaksd) / 2));
biasIdx = peaksd(1);
silIdx = peaksd(end) + 2 * fs;
cf = 100 / std(xrecord(peaksd(2):peaksd(end - 2)));
lresp = round(nto * 1.2);
tmpIdx = 1:lresp;
indivResp = zeros(lresp, n_mix);
indivSil = zeros(lresp, n_mix);
indivRef = zeros(lresp, n_mix);
copiedSilent = xr * 0;
xsilSegment = xr(1: fs);
for jj = 1:fs:length(xr) - fs
    copiedSilent((jj - 1) + (1:fs)) = xsilSegment;
end
for jj = 1:n_mix
    tmpResp = fftfilt(fvnSet(end:-1:1, jj), xr) * cf;
    tmpRespSil = fftfilt(fvnSet(end:-1:1, jj), copiedSilent) * cf;
    avResp = zeros(lresp, 1);
    avSil = zeros(lresp, 1);
    avRef = zeros(lresp, 1);
    for ii = signPeriod(n_mix) * 3 + 1:signPeriod(n_mix) * 5
        avResp = avResp + ...
            tmpResp(biasIdx - 100 + tmpIdx + (ii - 1) * nto, 1) * signMat(ii, jj);
%        avSil = avSil + ...
%            tmpResp(silIdx - 100 + tmpIdx + (ii - 1) * nto, 1) * signMat(ii, jj);
        avSil = avSil + ...
            tmpRespSil(biasIdx - 100 + tmpIdx + (ii - 1) * nto, 1) * signMat(ii, jj);
        switch exmode
            case 'measurement2ch'
                avRef = avRef + ...
                    tmpResp(biasIdx - 100 + tmpIdx + (ii - 1) * nto, 2) * signMat(ii, jj);
        end
    end
    indivResp(:, jj) = avResp / signPeriod(n_mix);
    indivSil(:, jj) = avSil / signPeriod(n_mix);
    switch exmode
        case 'measurement2ch'
            indivRef(:, jj) = avRef / signPeriod(n_mix);
    end
end
avResp = mean(indivResp, 2);
avSil = mean(indivSil, 2);

lrespview =  nto;
av_resp_th = ZoneThirdPowerSpectrum(mean(avResp(1:lrespview,:),2), fftl, fs);
sil_th_av = ZoneThirdPowerSpectrum(avSil(1:lrespview,:), fftl, fs);
referencePower = ZbiasSet(av_resp_th, fftl, fs);
figHandel = figure;
semilogx(fx(1:fftl / 2), 10 * log10(av_resp_th / referencePower), 'linewidth', 2);
grid on;
hold all
semilogx(fx(1:fftl / 2), 10 * log10(sil_th_av / referencePower), 'linewidth', 2);
%
if n_mix >= 4
    av_dresp_th = sil_th_av * 0;
    for ii = 1:n_mix
        dresp_th = ZoneThirdPowerSpectrum(indivResp(1:lrespview, ii) - avResp(1:lrespview), fftl, fs);
        av_dresp_th = av_dresp_th + dresp_th;
    end
    av_dresp_th = av_dresp_th / n_mix;
    semilogx(fx(1:fftl / 2), 10 * log10(av_dresp_th / referencePower), '--', 'linewidth', 2);
end
axis([20 fs / 2 -65 20]);
title([memoText '  n-FVN:' num2str(n_mix) '  tp:' num2str(to) ' (s)  SPL(A):' num2str(spl) ' dB']);
if n_mix >= 4
    legend('LIN.', 'equiv BG', 'NON-LIN.', 'location', 'best');
else
    legend('LIN.', 'equiv BG', 'location', 'best');
end
xlabel('frequency (Hz)');
ylabel('level (dB: rel. to average)');

recordNameBody = ['fvn44k' datestr(now, 30)];
printFileName = [recordNameBody '.eps'];
print(figHandel, '-depsc', printFileName);
waveFileName = [recordNameBody '.wav'];
audiowrite(waveFileName, xrecord, fs, 'bitspersample', 24);

if ~isempty(memoText)
    switch exmode
        case 'diagnosis'
            output.buffer = buffer;
            output.pinkBuffer = pinkBuffer;
            output.equalizedRecord = xr;
    end
end
output.averagedThirdOctaveResponse = av_resp_th;
if n_mix >= 4
    output.averagedThirdOctaveNonlinearComp = av_dresp_th;
end
output.averagedThirdOctaveBackground = sil_th_av;
output.averagedResponseWaveform = avResp;
output.individualResponseWaveform = indivResp;
output.averagedBackgroundWaveform = avSil;
output.individualBackgroundWaveform = indivSil;
switch exmode
    case 'measurement2ch'
        output.averagedReferenceWaveform = avRef;
        output.individualReferenceWaveform = indivRef;
end
output.memoText = memoText;
output.number_of_FVNs = n_mix;
output.responseLength = to;
output.soundPressureLevel = spl;
output.signalRange = biasIdx - 100 + signPeriod(n_mix) * 3 + [1 lresp];
output.backgroundRange = silIdx - 100 + signPeriod(n_mix) * 3 + [1 lresp];
output.samplingFrequency = fs;
output.modeOfOperation = exmode;
output.creationDate = datestr(now);
output.levelStructure = levelStr;
output.elapsedTime = toc(start_tic);
saveFileName = [recordNameBody '.mat'];
save(saveFileName, 'output');
output.recordedWaveform = xrecord;
end
%%

function output = ZoneThirdPowerSpectrum(x, fftl, fs)
fx = (0:fftl - 1) / fftl * fs;
fxH = fx(1:fftl / 2) * 2 ^ (1 / 6);
fxL = fx(1:fftl / 2) * 2 ^ (-1 / 6);
cumAvResp = cumsum(abs(fft(x, fftl)) .^ 2);
avRespH = interp1(fx, cumAvResp, fxH, 'linear', 'extrap');
avRespL = interp1(fx, cumAvResp, fxL, 'linear', 'extrap');
output = (avRespH - avRespL) ./ (fxH - fxL);
end

function output = ZbiasSet(oneThirdPoer, fftl, fs)
fx = (0:fftl - 1) / fftl * fs;
fxC = fx(1:fftl / 2);
fxCheck = 2 .^ (log2(10):1/6:log2(20000));
tmpPower = interp1(fxC, oneThirdPoer, fxCheck, 'linear', 'extrap');
output = sum(tmpPower) / length(fxCheck);
end