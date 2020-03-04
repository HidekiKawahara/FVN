function output = fvnTAFResponseAnalysysSG(fName, id1, id2, frep, nRep, fvnSet, displayOn)
%output = fvnTAFResponseAnalysysG(fName, fc, fv, id1, id2, frep, nRep, fvnSet, displayOn)
%% Test Scrpt for TAF-SIGMUS demo
%% Preparation for test
% fvnTAFResponseAnalysysG(fName, fc, fv, id1, id2, frep, nRep, fvnSet, displayOn)

%clear variables
%close all

%{
fName = 'tafSIGMUSSample.wav';
%fc = 200; % this one may not be necessary (can be calculated from the signal)
%fv = 200; % this one also may not be necessary
id1 = 9;
id2 = 13;
frep = 0.5; % repetion length in second
nRep = 8; % number of unit FVN repetitions
fvnSet = 'fvnMin200ms'; %16 sets of unit FVN responses
displayOn = 2; % 1: debug, 2:report generation
% Fetching necessary data from files.
%}

startTic = tic;
%[xx, fss] = audioread('taf220mod25fvn1a2aHRec4.wav');
[xx, fs] = audioread(fName);
tmp = load(fvnSet);
eval(['fvnMin200ms = tmp.' fvnSet ';']);
nto = round(frep * fs);
% Chect for the fundamental frequency of the stimulus
% First, I check the rms level for each 0.1 segment of the acquired voice and 
% the test stimuli, which is the loopback 

xref = xx(:, 2); % This is the loopback signal applied to the headphone
if displayOn == 3
    disp('missing fundamental')
    displayOn = 2;
    xref = [0; diff(xref .^ 2)];
end
xSp = xx(:, 1); % This is the acquired speech signal
nData = length(xref);
lSegment = round(0.1 * fs);
nSegment = floor(nData / lSegment);
rmsRef = zeros(nSegment, 1);
rmsSp = zeros(nSegment, 1);
for ii = 1:nSegment
    rmsRef(ii) = std(xref((ii - 1) * lSegment + (1:lSegment)));
    rmsSp(ii) = std(xSp((ii - 1) * lSegment + (1:lSegment)));
end
tSegment = (1:nSegment) * lSegment / fs;
if displayOn == 1
    figure
elseif displayOn == 2
    figureHandle = figure;
    set(figureHandle, 'position', [289  169  1111  816])
    subplot(321)
end
if displayOn > 0
    plot(tSegment, 20 * log10(rmsSp), 'linewidth', 2);
    hold all
    plot(tSegment, 20 * log10(rmsRef), 'linewidth', 2);
    set(gca, 'fontsize', 13, 'linewidth', 2);
    xlabel('time (s)')
    ylabel('level (dB rel. MSB)')
    grid on;
    axis([0 tSegment(end) -100 0])
    legend('voice', 'test stimulus', 'location', 'south')
    title(fName)
end
%% 
% Then, I select reliable initial 0.5 second of the test stimuli for determining 
% the corresponding note

thresholdRef = median(rmsRef) / 4;
tmpRmsRef = double(rmsRef > thresholdRef);
onPoints = tSegment(tmpRmsRef > tmpRmsRef([1 1:end - 1]));
offPoints = tSegment(tmpRmsRef > tmpRmsRef([2:end end]));
segmentLength = offPoints - onPoints;
[maxSegment, maxIdx] = max(segmentLength);
refSegment = [onPoints(maxIdx) offPoints(maxIdx)];
if maxSegment > 10
    refSample = xref(round((onPoints(maxIdx) + 0.2) * fs) + (1:round(fs / 2)));
    refStr = sourceAttributesAnalysis(refSample, fs, [1 length(refSample)], ...
        100, 1200, 4);
    rawRefFo = median(refStr.fixed_points_freq(:,1));
    noteIDNominal = round(log2(rawRefFo / 27.5) * 12);
    nominalRefFo = 27.5 * 2 ^ (noteIDNominal / 12);
    nominalOctaveID = floor((noteIDNominal -3 + 12) / 12);
    if displayOn == 1
        figure;
        plot(refStr.time_axis_wavelet, refStr.fixed_points_freq(:,1));
        grid on;
        hold all
        xlabel('segment time (s)')
        ylabel('fundamental frequency (Hz)')
    end
else
    error('Something wrong. Check!');
end
%% 
% Then also check for the fo for the voice

thresholdSp = median(rmsSp) / 4;
thresholdSp = median(rmsSp(rmsSp > thresholdSp)) / 4;
tmpRmsSp = double(rmsSp > thresholdSp);
onPoints = tSegment(tmpRmsSp > tmpRmsSp([1 1:end - 1]));
offPoints = tSegment(tmpRmsSp > tmpRmsSp([2:end end]));
segmentLength = offPoints - onPoints;
[maxSegment, maxIdx] = max(segmentLength);
spSegment = [onPoints(maxIdx) offPoints(maxIdx)];
if maxSegment > 10
    spSample = xSp(round((onPoints(maxIdx) + 0.2) * fs) + (1:round(fs / 2)));
    spStr = sourceAttributesAnalysis(spSample, fs, [1 length(spSample)], ...
        100, 1200, 4);
    rawSpFo = median(spStr.fixed_points_freq(:,1));
    if displayOn == 1
        plot(spStr.time_axis_wavelet, spStr.fixed_points_freq(:,1));
        grid on;
        xlabel('segment time (s)')
        ylabel('fundamental frequency (Hz)')
        axis([0 0.5 nominalRefFo * [2^(-1/3) 2^(1/3)]])
        legend('voice', 'test stimulus', 'location', 'best')
    end
else
    error('Something wrong. Check!');
end

noteNames = {'A', 'A', 'B', 'C', 'C', 'D', 'D', 'E', 'F', 'F', 'G', 'G'};
noteModList = {' ', '#', ' ', ' ', ' #', ' ', '#', ' ', ' ', '#', ' ', '#'};
noteName = [noteNames{rem(noteIDNominal, 12) + 1} num2str(nominalOctaveID) ...
    noteModList{rem(noteIDNominal, 12) + 1}];

fc = nominalRefFo;
fv = rawSpFo;
%%
fl = fc * 2^(-1/12);
fh = fc * 2^(1/12);
outputWavelet = designAnalyticWavelet(fs, fl, fh);
fvl = fv * 2^(-1/12);
fvh = fv * 2^(1/12);
outputWaveletv = designAnalyticWavelet(fs, fvl, fvh);

xxc = fftfilt(outputWavelet.wvlt(2).w, xx(:, 2));
xxcv = fftfilt(outputWaveletv.wvlt(2).w, xx(:, 1));
%%
idx = 1:length(xxc);
ndata = length(xxc);
xxrIf = angle(xxc(min(ndata, idx + 1)) ./ xxc) * fs / 2 / pi;
xxlIf = angle(xxcv(min(ndata, idx + 1)) ./ xxcv) * fs / 2 / pi;
idxP = min(ndata, idx + 1);

ttx = idx / fs;
xxlIfCentOrg = log2(abs(xxlIf) / fc) * 1200;
xxrIfCentOrg = log2(abs(xxrIf) / fc) * 1200;
if displayOn == 2 
    subplot(323);
elseif displayOn == 1
    figure
end
if displayOn == 1 || displayOn == 2
    tMargin = 0.4;
    refValid = ttx > refSegment(1) + tMargin & ttx < refSegment(2) - tMargin;
    spValid = ttx > spSegment(1) + tMargin & ttx < spSegment(2) - tMargin;
    plot(ttx(spValid), xxlIfCentOrg(spValid), 'linewidth', 2);
    hold all
    plot(ttx(refValid) , xxrIfCentOrg(refValid) , 'linewidth', 2);
    set(gca, 'xlim', ...
        [min(refSegment(1), spSegment(end)) max(refSegment(1), spSegment(end))])
    grid on
    set(gca,'fontsize', 14, 'linewidth', 2);
    ylabel('deviation (cent)');
    xlabel('time (s)');
    legend('voice', 'test stimulus', 'location', 'best')
    title(['target:' num2str(fc) ' Hz note:' noteName])
    drawnow
end
%%
xxlIfCent = xxlIfCentOrg .* spValid(:);
xxlIfCent(isnan(xxlIfCent)) = 0;
lFvn = length(fvnMin200ms(:, id1));
lFvnHlf = floor(lFvn / 2);
mon1L220 = fftfilt(fvnMin200ms(end:-1:1, id1), [xxlIfCent; zeros(lFvn, 1)]);
mon2L220 = fftfilt(fvnMin200ms(end:-1:1, id2), [xxlIfCent; zeros(lFvn, 1)]);
xxrIfCent = xxrIfCentOrg .* refValid(:);
xxrIfCent(isnan(xxrIfCent)) = 0;
mon1R220 = fftfilt(fvnMin200ms(end:-1:1, id1), [xxrIfCent; zeros(lFvn, 1)]);
mon2R220 = fftfilt(fvnMin200ms(end:-1:1, id2), [xxrIfCent; zeros(lFvn, 1)]);
mon1L220 = mon1L220(lFvnHlf + idx);
mon2L220 = mon2L220(lFvnHlf + idx);
mon1R220 = mon1R220(lFvnHlf + idx);
mon2R220 = mon2R220(lFvnHlf + idx);

basepx = min(ndata, idx + nto);

if displayOn == 1 || displayOn == 2
    if displayOn == 1
        figure
    elseif displayOn == 2
        subplot(422);
    end
    plot(ttx, (mon1L220 + mon1L220(basepx)) / 2);grid on;
    hold all
    plot(ttx, (mon2L220 - mon2L220(basepx)) / 2);
    set(gca, 'xlim', ttx([1 end]), 'linewidth', 2);
    set(gca,'fontsize', 14);
    ylabel('deviation (cent)');
    if displayOn == 1
        figure
    elseif displayOn == 2
        subplot(424);
    end
    plot(ttx, (mon1R220 + mon1R220(basepx)) / 2);grid on;
    hold all
    plot(ttx, (mon2R220 - mon2R220(basepx)) / 2);
    set(gca, 'xlim', ttx([1 end]), 'linewidth', 2);
    set(gca,'fontsize', 14);
    ylabel('deviation (cent)');
    xlabel('time (s)');
    drawnow
end
%%
initialSync1 = mon1R220 + mon1R220(basepx);
initialSync2 = mon2R220 - mon2R220(basepx);
maxV = max(initialSync1 - initialSync2);
peakThreshold = maxV * 0.9;
cleanSync = (initialSync1 - initialSync2) .* ((initialSync1 - initialSync2) > peakThreshold);
idxM = max(1, idx - 1);
maxIdx = min(idx(cleanSync > cleanSync(idxP) & cleanSync > cleanSync(idxM)));
%%
%t0 = 14.5;
%nt0 = round(t0 * fs);
nt0 = maxIdx - 4410;
idxTmp = 1:round(fs * 2 * frep);

csignList = [-1 1 -1 1 -1 1];
stdRes = zeros(4, 2);
responseBody = ttx < 0.99;
for jj = 1:4
    avResp1 = zeros(round(fs * 2 * frep), 1);
    avResp2 = zeros(round(fs * 2 * frep), 1);
    avRef1 = zeros(round(fs * 2 * frep), 1);
    avRef2 = zeros(round(fs * 2 * frep), 1);
    csign = csignList(jj);
    for ii = (1:nRep) + jj - 1
        avResp1 = avResp1 + mon1L220(nt0 + (ii - 1) * nto + idxTmp);
        avResp2 = avResp2 + mon2L220(nt0 + (ii - 1) * nto + idxTmp) * csign;
        avRef1 = avRef1 + mon1R220(nt0 + (ii - 1) * nto + idxTmp);
        avRef2 = avRef2 + mon2R220(nt0 + (ii - 1) * nto + idxTmp) * csign;
        csign = -1 * csign;
    end
    stdRes(jj, 1) = std(avResp1(responseBody) + avResp2(responseBody));
    stdRes(jj, 2) = std(avResp1(responseBody) - avResp2(responseBody));
end
if displayOn == 1
    figure
    plot(stdRes(:, 1) ./ stdRes(:, 2))
end
[~, bestIdx] = max(stdRes(:, 1) ./ stdRes(:, 2));
avResp1 = zeros(round(fs * 2 * frep), 1);
avResp2 = zeros(round(fs * 2 * frep), 1);
avRef1 = zeros(round(fs * 2 * frep), 1);
avRef2 = zeros(round(fs * 2 * frep), 1);
csign = csignList(bestIdx);
for ii = (1:nRep) + bestIdx - 1
    avResp1 = avResp1 + mon1L220(nt0 + (ii - 1) * nto + idxTmp);
    avResp2 = avResp2 + mon2L220(nt0 + (ii - 1) * nto + idxTmp) * csign;
    avRef1 = avRef1 + mon1R220(nt0 + (ii - 1) * nto + idxTmp);
    avRef2 = avRef2 + mon2R220(nt0 + (ii - 1) * nto + idxTmp) * csign;
    csign = -1 * csign;
end
biasChkCondition = ttx < 0.15 ;% | ttx > 0.5 & ttx < 0.6;
avResp1Bias = mean(avResp1(biasChkCondition));
avResp2Bias = mean(avResp2(biasChkCondition));

%%
if displayOn == 2 || displayOn == 1
    if displayOn == 1
        figure;
        set(gcf, 'position', [680 558 560 250])
    else
        subplot(325);
    end
    plot(idxTmp / fs, (avResp1 - avResp1Bias) / nRep, 'linewidth', 2);grid on;
    hold all
    plot(idxTmp / fs, (avResp2 - avResp2Bias) / nRep, 'linewidth', 2)
    set(gca, 'xlim', idxTmp([1 end]) / fs, 'linewidth', 2);
    set(gca,'fontsize', 14);
    xlabel('response time (s)')
    ylabel('deviation (cent)')
    if displayOn == 1
        figure;
    else
        subplot(224);
    end
    plot(idxTmp / fs, (avResp2 + avResp1 - avResp1Bias - avResp2Bias) / nRep / 2, ...
        'linewidth', 2); grid on;
    hold all
    plot(idxTmp / fs, (avRef2 + avRef1) / nRep / 2, 'linewidth', 2);
    set(gca, 'xlim', idxTmp([1 end]) / fs, 'linewidth', 2);
    set(gca,'fontsize', 14);
    xlabel('response time (s)')
    ylabel('deviation (cent)')
    drawnow
end
%%
cfBin = [1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 1 1 -1 -1];
avRefRest = zeros(round(fs * 2 * frep), 1);
avRespRest = zeros(round(fs * 2 * frep), 1);
%figure;
for kk = 1:16
    if kk ~= id1 && kk~= id2
        mon3R220 = fftfilt(fvnMin200ms(end:-1:1, kk), [xxlIfCent; zeros(lFvn, 1)]);
        mon3L220 = fftfilt(fvnMin200ms(end:-1:1, kk), [xxlIfCent; zeros(lFvn, 1)]);
        mon3R220 = mon3R220(lFvnHlf + idx);
        mon3L220 = mon3L220(lFvnHlf + idx);
        avRef3 = zeros(round(fs * 2 * frep), 1);
        avResp3 = zeros(round(fs * 2 * frep), 1);
        for ii = (1:nRep) + bestIdx - 1
            avRef3 = avRef3 + mon3R220(nt0 + (ii - 1) * nto + idxTmp) * cfBin(ii - 1 + 2);
            avResp3 = avResp3 + mon3L220(nt0 + (ii - 1) * nto + idxTmp) * cfBin(ii - 1 + 2);
        end
        avRefRest = avRefRest + (avRef3 / nRep) .^ 2;
        avRespRest = avRespRest + (avResp3 / nRep) .^ 2;
        %plot(avRef3/4)
        %hold all
    end
end
avRefRest = sqrt(avRefRest / 14);
avRespRest = sqrt(avRespRest / 14);

if displayOn == 2 || displayOn == 1
    if displayOn == 1
        hold all
    else
        subplot(224);
        hold all
    end
    baseResp = (avResp2 + avResp1 - avResp1Bias - avResp2Bias) / nRep / 2;
    plot(idxTmp / fs, baseResp + avRespRest, 'k:', 'linewidth', 2);
    plot(idxTmp / fs, baseResp - avRespRest, 'k:', 'linewidth', 2);
    set(gca, 'xlim', idxTmp([1 end]) / fs, 'linewidth', 2);
    legend('response', 'stimulus', '+sigma', '-sigma');
    if displayOn == 2
    print(figureHandle, '-dpng', '-r200', [fName(1:end - 4) '.png']);
    end
end
%%
output.avResp1 = avResp1 / nRep;
output.avResp2 = avResp2 / nRep;
output.avRef1 = avRef1 / nRep;
output.avRef2 = avRef2 / nRep;
output.avRefRest = avRefRest;
output.avRespRest = avRespRest;
output.respTimeAxis = idxTmp / fs;
output.signalTimeAxis = ttx;
output.nIteration = 6;
output.fName = fName;
output.fc = fc;
output.elapsedTime = toc(startTic);
end