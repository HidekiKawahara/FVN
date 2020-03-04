%% Test event localization
% by Hideki Kawahara
% 25/Jan./2020

%matName = 'fvn44k20200123T224712.mat';
matName = 'fvn44k20200125T053950.mat';
tmp = load(matName);

hh = tmp.output.averagedResponseWaveform;
fs = tmp.output.samplingFrequency;
hl = length(hh);

fc1 = 10000;
n1 = round(fs * 0.001);
w1e = hanning(2 * n1 + 1);
tt = (-n1:n1)' / fs;
w1 = w1e .* exp(1i * tt * 2 * pi * fc1);
tmpout = fftfilt(w1, hh);
tmpoute = abs(tmpout(n1 + 1:end)) .^ 2;
[~, idxMax] = max(tmpoute);

%%
w2 = blackman(2 * n1 + 1);
baseIdx = (-n1:n1)';
hhseg = hh(idxMax + baseIdx);
fftl = 1024;
hhF = fft(hhseg .* w2, fftl);
fx = (0:fftl - 1) / fftl * fs;

idxC = 1:fftl;
idxF = [2:fftl fftl];
idxB = [1 1:fftl - 1];
dphase = fx * (n1 / fs * 2 * pi);
tmpgd = -(angle(hhF(idxF) ./ hhF) + angle(hhF ./ hhF(idxB))) / 2 / (2 * pi * fx(2));

%%
hhACep = ifft(log(abs(hhF)));
hhCCep = hhACep;
hhCCep(2:fftl/2) = 2 * hhCCep(2:fftl/2);
hhCCep(fftl/2+1:end) = 0;
hhFC = exp(fft(hhCCep));
tmpPgd = -(angle(hhFC(idxF) ./ hhFC) + angle(hhFC ./ hhFC(idxB))) / 2 / (2 * pi * fx(2));

fxSel = fx < 20000;
gdFix = tmpgd * fs - tmpPgd * fs - n1;
figure;plot(fx, (tmpgd - tmpPgd  - n1 / fs) * 10^6);grid on
axis([0 fs/2 [min(gdFix/fs) max(gdFix/fs)] * 10^6]);
set(gca, 'fontsize', 15);
xlabel('frequency (Hz)');
ylabel('group delay (micro second)');

figure;semilogx(fx, (tmpgd  - tmpPgd - n1/fs) * 10^6);grid on
axis([100 fs/2 [min(gdFix/fs) max(gdFix/fs)] * 10^6]);
set(gca, 'fontsize', 15);
xlabel('frequency (Hz)');
ylabel('group delay (micro second)');
title('after causality compensation');

figure;semilogx(fx, (tmpgd -n1/fs)*10^6); grid on;
hold all;
semilogx(fx, tmpPgd * 10^6);
axis([100 fs/2 [ min(min(tmpgd(fxSel) -n1/fs),  min(tmpPgd(fxSel))), ...
    max(max(tmpgd(fxSel) -n1/fs), max(tmpPgd(fxSel)))] * 10^6 ]);
set(gca, 'fontsize', 15);
xlabel('frequency (Hz)');
ylabel('group delay (micro second)');
legend('measured', 'causality');

tt = ((1:length(hh)) - idxMax) / fs;
pw = abs(hhF) .^ 2;
selectorr = fx > 6000 & fx < 14000;
meanPos = sum(pw(selectorr) .* gdFix(selectorr)) / sum(pw(selectorr)) / fs;
figure;
plot(tt * 10^6, hh/ max(abs(hh)));grid on;
hold all
plot(baseIdx / fs * 10^6, w2 * max(hh)/ max(abs(hh)));
plot(meanPos * [1 1] * 10^6, 1.1 * [min(hh) max(hh)]/ max(abs(hh)), 'g')
axis([[-n1/fs, 0.0015] * 10^6 1.1 * [min(hh) max(hh)] / max(abs(hh))])
set(gca, 'fontsize', 15);
xlabel('time (micro second)');
legend('impulse response', 'time window', 'event origin');

%%

tth = (-n1:n1)' / fs;
w1p = w1e .^ 2;
w1pn = w1p/sum(w1p);
w1w = -(w1p .* tth) / sum(w1p);
hhp = abs(fftfilt(w1, hh)) .^ 2;
mthhp = fftfilt(tth, hhp);
pphhp = fftfilt(ones(length(w1), 1), hhp);

fixMthhp = fftfilt(tth / (2*n1+1), hhp);


