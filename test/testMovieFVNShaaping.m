%% movie generation for EA2019-Oct
% by Hideki Kawahaaraa
% 29/Oct./2019

clear variables
close all

fs = 44100;
fftl = 32768 * 64;

rng(13214)
pdiv = 2;
fdiv_av = 3;%13;
baseIdx = 1:fftl;
centerIdx = fftl / 2 + (-4:4);
fringeIdx = baseIdx((baseIdx < fftl / 2 - 4 & baseIdx > fftl / 2 - 10) ...
    | (baseIdx > fftl / 2 + 4 & baseIdx < fftl / 2 + 10));
bw_c = fdiv_av * 2;
output = generateFVN6real(fs, fftl, bw_c, fdiv_av, pdiv);
xFVN = output.xFVN;

dph = output.mixedPhase([(2:end) 1]) - output.mixedPhase([end (1:end - 1)]);
xtmp = (-1:0.001:1) * max(abs(dph));
fx = (0:fftl-1)'/fftl*fs;
figure;
for slopec = [-50, -10, -25]
    outdir = ['mov' datestr(now,30)];
    mkdir(outdir)
    icount = 1;
    for biasc = -2.5:0.01:2.5
        wdph = spsigmoid(dph, slopec, max(abs(dph)) * biasc); %(1.0 ./ (1 + exp(-110 * dph))) .^ 0.3;
        nufreq = cumsum(wdph);
        nufreq = nufreq / nufreq(end) * length(nufreq);
        exPhase = interp1(nufreq, output.mixedPhase, 1:length(nufreq), 'linear', 'extrap');
        wrpH = real(ifft(exp(1i * exPhase)));
        wrpH_reg = wrpH / max(abs(wrpH));
        tDisp = ((1:fftl) - fftl / 2) / fs;

        hfvn = fftshift(wrpH_reg);
        cumHfvn = cumsum(abs(hfvn)+0.00000000001)/sum(abs(hfvn)+0.00000000001);
        xlim = round(interp1(cumHfvn, 1:length(cumHfvn), [0.00003 1-0.00003]));

        subplot(221);
        plot(tDisp(xlim(1):xlim(2)), 20*log10(abs(hfvn(xlim(1):xlim(2))/max(abs(wrpH_reg)))));grid on
        set(gca, 'fontsize', 14)
        axis([tDisp(xlim(1)) tDisp(xlim(2)) -120 0])
        xlabel('time (s)');
        ylabel('amplitude (rel. dB)');
        title(['slope:' num2str(slopec) '  bias:' num2str(biasc, '%5.2f')])
        subplot(223);
        plot(xtmp, spsigmoid(xtmp, slopec, max(abs(xtmp)) * biasc));grid on
        set(gca, 'ylim', [min(spsigmoid(xtmp, slopec, max(abs(xtmp)) * biasc)), ...
            max(spsigmoid(xtmp, slopec, max(abs(xtmp)) * biasc))]);
        subplot(222);
        plot(fx,exPhase);grid on;
        set(gca, 'xlim', [900 1100]);grid on;
        set(gca, 'ylim', [min(exPhase(fx>900 & fx<1100)) max(exPhase(fx>900 & fx<1100))]);
        ylabel('phase (radian)');
        subplot(224);
        plot(fx(1:end-1),-diff(exPhase) / fx(2) / 2 / pi);grid on;
        ylabel('group delay (s)')
        xlabel('frequency (Hz)');
        set(gca, 'xlim', [900 1100]);grid on;
        set(gca, 'ylim', [min(-diff(exPhase(fx>900 & fx<1100))) max(-diff(exPhase(fx>900 & fx<1100)))]/ fx(2) / 2 / pi);
        drawnow
        print('-dpng','-r100', [outdir '/' num2str(icount, '%04d') '.png']);
        icount = icount + 1;
    end
end
