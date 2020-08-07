function output = tripletMeasurement(channelID, calLeveldB, outPath, nChannel)
%%
[xRef, fs] = audioread("testTripleFvn44Pink.wav");
nRepat = 44;
nto = round(0.2 * fs);
tmpOut = aweightTable;
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
calibCf = 10 ^ (calLeveldB / 20);
if nChannel == 2
xRef = xr(:, 2) * calibCf;
xRec = xr(:, 1) * calibCf;
else
    xRec = xr(:, 1) * calibCf;
    xRef = [];
end
outputQ = quatroFVNresponseAnalysisQ(fvnSet, xRec, xRef, fs, nRepat, nto);
weightFilt = weightingFilter;
xA = weightFilt(y(:, 1)) * calibCf;
spl = 20 * log10(std(xA(round(length(xRec) / 2) + (-nto:nto))));
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
title("LAeq: " + num2str(spl) + " dB, with 1/3 octave smoothing")
set(gca, "FontSize", 15, 'LineWidth', 2)
xlabel("frequency (Hz)")
ylabel("level (dB)")
drawnow
% ---- add analysis conditions to results ---
outputQ.fvnConditions.samplingFrequency = fs;
outputQ.fvnConditions.nRepeat = nRepat;
outputQ.fvnConditions.repetitionInterval = nto;
outputQ.fvnConditions.fvnSet = fvnSet;
outputQ.channelID = channelID;
outputQ.dB2convertSPL = calLeveldB;
% ----
tmpFname = ['fvnResp' datestr(now, 30) 'LAeq'];
foName = string(tmpFname)  + num2str(round(spl)) + "dBtst.eps";
print(fhgHandle, "-depsc", [outPath '/' char(foName)])
fMatName = string(tmpFname)  + num2str(round(spl)) + "dBtst.mat";
save([outPath '/' char(fMatName)], "outputQ");
output.recordNameBody = [tmpFname num2str(round(spl))  'dBtst'];
output.analysisStr = outputQ;
end

