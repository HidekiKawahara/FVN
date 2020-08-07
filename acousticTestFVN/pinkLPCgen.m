function output = pinkLPCgen(fs)
%fs = 44100;
sg = 0.4; % 200 ms period is 400 ms
fftl = 2 ^ ceil(log2(fs * 15 * sg));
lowlimitt = 40;
fx = (0:fftl - 1)' / fftl * fs;
pinkgain = 1.0 ./ sqrt(fx);
pinkgain(fx < lowlimitt) = max(pinkgain(fx >= lowlimitt));
pinkgain(fftl:-1:fftl / 2 + 1) = pinkgain(2:fftl / 2 + 1);
pinkgain = pinkgain / pinkgain(1);
rawImpulsePink = real(ifft(pinkgain)); % it is real already but for safety
% all pole approximation
rawAutoCorr = real(ifft(pinkgain .^ 2));
pList = 20:5:45;
figure;
for ii = 1:length(pList)
    [a, ee] = levinson(rawAutoCorr, pList(ii));
    lpcGain = 1 ./ abs(fft(a/sum(a), fftl));
    semilogx(fx, 20*log10(lpcGain(:) ./ pinkgain));grid on
    hold all
end
output.pinkLPC = a;
end

