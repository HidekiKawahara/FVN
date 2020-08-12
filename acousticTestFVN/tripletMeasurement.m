function output = tripletMeasurement(varargin)
% FVN based acoustic analysis
%  output = tripletMeasurement(channelID, calLeveldB, outPath, nChannel, synchSource, compactLog)
%  output = tripletMeasurement(fullPathName)
%
%  Input argument
%    channelID      : output signal channel: "L" or "R"
%    calLeveldB     : calibration data: dB to convert internal level to SPL
%    outPath        : file output directory; for example, "~/tmp"
%    nChannel       : number of channels of AD converter
%    synchSource    : 'Internal' : (default) generate ref. signal internally
%                   : 'R-Ch.' : uses R-channel signal to synchronize
%                   : 'Internal2' : generate ref. signal internally and
%                                   analyzes both (L and R) input channels
%    compactLog     : 1  : (default) compact reusable mat-file output and
%                   :      eps-file of analysis visualization
%                   : 0  : mat-file with all analysis results and
%                   :      eps-file of analysis visualization
%                   : -1 : no file output
%    fullPathName   : Analyse stored data using compactLog = 1
%
%  Output
%    output : structure variable with the following fields
%      recordNameBody   : File name root of output files (.eps, .mat)
%      analysisStr      : structure consisting of analysis results (ch-1 in)
%      output2          : (only 'Internal2' mode) strucrure with fields
%        recordNameBody : same above for (ch-2 in)
%        analysisStr    : same above for (ch-2 in)
%
%  Please refer the function "quatroFVNresponseAnalysisQ" for details of
%  analysis results

%   Copyright 2020 Hideki Kawahara
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
   
   
switch nargin
    case 6
        channelID = varargin{1};
        calLeveldB = varargin{2};
        outPath = varargin{3};
        nChannel = varargin{4};
        synchSource = varargin{5};
        compactLog = varargin{6};
    case 1
        eval(['load ' varargin{1}]);
        compactLog = -1;
        outPath = [];
end
fvnSequenceName = 'testTripleFvn44Pink.wav';
[xRefWav, fs] = audioread(fvnSequenceName);
nRepeat = 44;
nto = round(0.2 * fs);
outputLPC = pinkLPCgenM(fs);
fvnSetFile = 'fvnMin100ms.mat';
tmp = load(fvnSetFile);
fvnMin100ms = tmp.fvnMin100ms;
fvnSet = fvnMin100ms(:, 1:4);
switch channelID
    case 'L'
        testSig = [xRefWav, xRefWav * 0];
    case 'R'
        testSig = [xRefWav * 0, xRefWav];
end
if nargin == 6
    audioPlayer = audioplayer(testSig, fs, 24);
    audioRecorder = audiorecorder(fs, 24, nChannel);
    dataDuration = length(xRefWav) / fs;
    record(audioRecorder);
    pause(0.3);
    play(audioPlayer);
    pause(dataDuration + 1);
    y = getaudiodata(audioRecorder);
    stop(audioRecorder);
end
xr = filter(outputLPC.pinkLPC, 1, y);
calibCf = 10 .^ (calLeveldB / 20);
if nChannel == 2
    switch synchSource
        case 'R-Ch.'
            xRef = xr(:, 2) * calibCf(1);
            xRec = xr(:, 1) * calibCf(1);
            analysisOut = analyzeTriplFVN(fvnSet, xRec, xRef, fs, nRepeat, ...
                nto, y, calibCf(1), channelID, calLeveldB);
            output = reportTheResults(analysisOut, outPath, compactLog);
        case {'Internal','Internal2'}
            xRec = xr(:, 1) * calibCf(1);
            xRef = [];
            analysisOut = analyzeTriplFVN(fvnSet, xRec, xRef, fs, nRepeat, ...
                nto, y, calibCf(1), channelID, calLeveldB);
            output = reportTheResults(analysisOut, outPath, compactLog);
            switch synchSource
                case 'Internal2'
                    xRec2 = xr(:, 2) * calibCf(2);
                    analysisOut2 = analyzeTriplFVN(fvnSet, xRec2, xRef, fs, ...
                        nRepeat, nto , y(:, 2), calibCf(2), channelID, calLeveldB);
                    output2 = reportTheResults(analysisOut2, outPath, compactLog);
                    output.output2 = output2;
            end
    end
else
    xRec = xr(:, 1) * calibCf(1);
    xRef = [];
    analysisOut = analyzeTriplFVN(fvnSet, xRec, xRef, fs, nRepeat, ...
        nto, y, calibCf(1), channelID, calLeveldB);
    output = reportTheResults(analysisOut, outPath, compactLog);
end
if compactLog == 1
    fMatName = output.recordNameBody + "C.mat";
    save([outPath '/' char(fMatName)], "fvnSequenceName", "fvnSetFile", ...
        "synchSource", "nChannel", "y", "outputLPC", "channelID", "calLeveldB");
end
end

function output = reportTheResults(analysisOut, outPath, compactLog)
outputQ = analysisOut.outputQ;
if compactLog >= 0
    tmpFname = ['fvnResp' datestr(now, 30) 'LAeq'];
    foName = string(tmpFname)  + num2str(round(analysisOut.spl)) + "dBtst.eps";
    print(analysisOut.fhgHandle, "-depsc", [outPath '/' char(foName)])
    if ~compactLog
        fMatName = string(tmpFname)  + num2str(round(analysisOut.spl)) + "dBtst.mat";
        save([outPath '/' char(fMatName)], "outputQ");
    end
    output.recordNameBody = [tmpFname num2str(round(analysisOut.spl))  'dBtst'];
else
    output.recordNameBody = [];
end
output.analysisStr = outputQ;
end

function analysisOut = analyzeTriplFVN(fvnSet, xRec, xRef, fs, nRepeat, nto, ...
    y, calibCf, channelID, calLeveldB)
stTic = tic;
tmpOut = aweightTable;
outputQ = quatroFVNresponseAnalysisQ(fvnSet, xRec, xRef, fs, nRepeat, nto);
weightFilt = weightingFilter;
xA = weightFilt(y(:, 1)) * calibCf;
if calLeveldB(1) == 0
    spl = 80;
else
    spl = 20 * log10(std(xA(round(length(xRec) / 2) + (-nto:nto))));
end
thridLvl = interp1(outputQ.frequencyAxis, 10 .^ (outputQ.smavGain / 10), tmpOut(:,1), 'linear',"extrap");
splInnr = 10*log10(sum(thridLvl(:) .* 10 .^ (tmpOut(:,2) / 10)));
calibSpl = spl - splInnr;
fhgHandle = figure;
semilogx(outputQ.frequencyAxis, outputQ.smavGain + calibSpl, 'LineWidth', 2);grid on;
hold all;semilogx(outputQ.frequencyAxisEx, outputQ.smextGain + calibSpl, 'LineWidth', 2);grid on;
hold all;semilogx(outputQ.frequencyAxis, outputQ.smbackGroundSpec + calibSpl, 'LineWidth', 2);grid on;
hold all;semilogx(outputQ.frequencyAxis, outputQ.smprecedingSpec + calibSpl, 'LineWidth', 2);grid on;
hold all;semilogx(outputQ.frequencyAxis, outputQ.smnonLinearSpec + calibSpl, 'LineWidth', 2);grid on;
axis([10 fs/2 10 85])
legend('lin.', 'lin.Ex', 'BG.est', 'BG.pre', 'nonLin.', 'Location', "south", "numcolumns", 5)
set(gca, "FontSize", 15, 'LineWidth', 2)
xlabel("frequency (Hz)")
ylabel("Sound pressure level (dB, ref. 20\mu Pa)")
if calLeveldB(1) == 0
    title("Presumed LAeq 80dB with 1/3 octave smoothing")
else
    title("LAeq: " + num2str(spl,'%6.2f') + " dB, with 1/3 octave smoothing")
end
drawnow
% ---- add analysis conditions to results ---
outputQ.fvnConditions.samplingFrequency = fs;
outputQ.fvnConditions.nRepeat = nRepeat;
outputQ.fvnConditions.repetitionInterval = nto;
outputQ.fvnConditions.fvnSet = fvnSet;
outputQ.channelID = channelID;
outputQ.dB2convertSPL = calLeveldB;
analysisOut.elapsedTime = toc(stTic);
analysisOut.spl = spl;
analysisOut.outputQ = outputQ;
analysisOut.fhgHandle = fhgHandle;
end