function output = tripletMeasurement(channelID, calLeveldB, outPath, nChannel, synchSource)
%
[xRef, fs] = audioread("testTripleFvn44Pink.wav");
nRepat = 44;
nto = round(0.2 * fs);
outputLPC = pinkLPCgenM(fs);
tmp = load('fvnMin100ms.mat');
fvnMin200ms = tmp.fvnMin100ms;
fvnSet = fvnMin200ms(:, 1:4);
switch channelID
    case 'L'
        testSig = [xRef, xRef * 0];
    case 'R'
        testSig = [xRef * 0, xRef];
end
audioPlayer = audioplayer(testSig, fs, 24);
audioRecorder = audiorecorder(fs, 24, nChannel);
dataDuration = length(xRef) / fs;
record(audioRecorder);
pause(0.3);
play(audioPlayer);
pause(dataDuration + 1);
y = getaudiodata(audioRecorder);
stop(audioRecorder);
xr = filter(outputLPC.pinkLPC, 1, y);
calibCf = 10 .^ (calLeveldB / 20);
if nChannel == 2
    switch synchSource
        case 'R-Ch.'
            xRef = xr(:, 2) * calibCf(1);
            xRec = xr(:, 1) * calibCf(1);
            analysisOut = analyzeTriplFVN(fvnSet, xRec, xRef, fs, nRepat, nto, y, calibCf(1), channelID, calLeveldB);
            output = reportTheResults(analysisOut, outPath);
        case {'Internal','Internal2'}
            xRec = xr(:, 1) * calibCf(1);
            xRef = [];
            analysisOut = analyzeTriplFVN(fvnSet, xRec, xRef, fs, nRepat, nto, y, calibCf(1), channelID, calLeveldB);
            output = reportTheResults(analysisOut, outPath);
            switch synchSource
                case 'Internal2'
                    xRec2 = xr(:, 2) * calibCf(2);
                    analysisOut2 = analyzeTriplFVN(fvnSet, xRec2, xRef, fs, nRepat, nto , y(:, 2), calibCf(2), channelID, calLeveldB);
                    output2 = reportTheResults(analysisOut2, outPath);
                    output.output2 = output2;
            end
    end
else
    xRec = xr(:, 1) * calibCf(1);
    xRef = [];
    analysisOut = analyzeTriplFVN(fvnSet, xRec, xRef, fs, nRepat, nto, y, calibCf(1), channelID, calLeveldB);
    output = reportTheResults(analysisOut, outPath);
end
end

function output = reportTheResults(analysisOut, outPath)
outputQ = analysisOut.outputQ;
tmpFname = ['fvnResp' datestr(now, 30) 'LAeq'];
foName = string(tmpFname)  + num2str(round(analysisOut.spl)) + "dBtst.eps";
print(analysisOut.fhgHandle, "-depsc", [outPath '/' char(foName)])
fMatName = string(tmpFname)  + num2str(round(analysisOut.spl)) + "dBtst.mat";
save([outPath '/' char(fMatName)], "outputQ");
output.recordNameBody = [tmpFname num2str(round(analysisOut.spl))  'dBtst'];
output.analysisStr = outputQ;
end

function analysisOut = analyzeTriplFVN(fvnSet, xRec, xRef, fs, nRepat, nto, y, calibCf, channelID, calLeveldB)
stTic = tic;
tmpOut = aweightTable;
outputQ = quatroFVNresponseAnalysisQ(fvnSet, xRec, xRef, fs, nRepat, nto);
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
outputQ.fvnConditions.nRepeat = nRepat;
outputQ.fvnConditions.repetitionInterval = nto;
outputQ.fvnConditions.fvnSet = fvnSet;
outputQ.channelID = channelID;
outputQ.dB2convertSPL = calLeveldB;
analysisOut.elapsedTime = toc(stTic);
analysisOut.spl = spl;
analysisOut.outputQ = outputQ;
analysisOut.fhgHandle = fhgHandle;
end