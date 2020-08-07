function output = quatroFVNresponseAnalysisQ(fvnSet, xRec, xRef, fs, nRepat, nto)
%%
if length(xRec) > fs
    startTic = tic;
    B4 = [1  1  1  1  1  1  1  1; ...
        1 -1  1 -1  1 -1  1 -1; ...
        1  1 -1 -1  1  1 -1 -1; ...
        1  1  1  1 -1 -1 -1 -1];
    if length(xRef) < fs
        outputF = generateQuadFVNseq(fvnSet, nRepat, nto);
        tmpxRef = outputF.xTestR;
        pwRef = tmpxRef .^2;
        smpwRef = sqrt(abs(fftfilt(hanning(nto*2)/sum(hanning(nto*2)), pwRef)));
        pwRec = xRec .^ 2;
        smpwRec = sqrt(abs(fftfilt(hanning(nto*2)/sum(hanning(nto*2)), pwRec)));
        smpwRef = smpwRef / max(smpwRef);
        smpwRec = smpwRec / max(smpwRec);
        
        nRec = length(smpwRec);
        nRef = length(smpwRef);
        if nRef > nRec
            nRef = round(nRef/2);
        end
        diffList = zeros(nRec - nRef, 1);
        for ii = 1:nRec - nRef
            diffList(ii) = sum(abs(smpwRec((1:100:nRef) - 1 + ii) - smpwRef(1:100:nRef)) .^2);
        end
        [~, biasIdx] = min(diffList);
        xRef = xRec * 0;
        xRef((1:nRef) - 1 + biasIdx - 2000) = tmpxRef;    
    end
    fftl = size(fvnSet, 1);
    comRec = fftfilt(fvnSet(end:-1:1, :), [xRec(:, 1);zeros(fftl, 1)]);
    bufRec = zeros(size(comRec, 1) + nto * nRepat, 4);
    comRef = fftfilt(fvnSet(end:-1:1, :), [xRef(:, 1);zeros(fftl, 1)]);
    bufRef = zeros(size(comRec, 1) + nto * nRepat, 4);
    baseIdx = 1:size(comRec, 1);
    for ii = 1:8
        kk = rem((ii - 1), 8) + 1;
        bufRec((ii - 1) * nto + baseIdx, :) = bufRec((ii - 1) * nto + baseIdx, :) + comRec * diag(B4(:, kk));
        bufRef((ii - 1) * nto + baseIdx, :) = bufRef((ii - 1) * nto + baseIdx, :) + comRef * diag(B4(:, kk));
    end
    chkIdx = (1:length(bufRef(:,4)))';
    %maxV = max(bufRef(:,4));
    w=hanning(441);
    qPower = fftfilt(w/sum(w), bufRef(:,4) .^ 2);
    repPoint = length(xRef)/2+fftl/2;
    baseLvl = 10*log10(mean(qPower(round(repPoint) + (-fftl/4:fftl/4))));
    maxLvL = 10*log10(max(qPower));
    thLvl = (maxLvL - baseLvl) /4 + baseLvl;
    safeStart = max(chkIdx(chkIdx < repPoint & 10*log10(qPower) > thLvl)) + round(0.3 * fs);
    safeEnd = min(chkIdx(chkIdx > repPoint & 10*log10(qPower) > thLvl)) - round(0.3 * fs);
    %% ---- background noise level check ---
    nK = floor((safeEnd - safeStart) / (8 * nto));
    rRV = zeros(8*nto, 1);
    for k = 0:nK - 1
        rRV = rRV + bufRec(safeStart + (1:8 * nto) + 8 * nto * k, 4) / 8;
    end
    rRV = rRV * sqrt(8 / nK);
    %% ---- averaged impulse response ---
    bufRef = bufRef / max(bufRef(:)) * 8;
    rawPeaks1 = chkIdx(bufRef(:, 1) / 8 > 0.99 & chkIdx > safeStart & chkIdx + 8 * nto < safeEnd) - 44;
    rawPeaks2 = chkIdx(bufRef(:, 2) / 8 > 0.99 & chkIdx > safeStart & chkIdx + 8 * nto < safeEnd) - 44;
    rawPeaks3 = chkIdx(bufRef(:, 3) / 8 > 0.99 & chkIdx > safeStart & chkIdx + 8 * nto < safeEnd) - 44;
    resp1 = zeros(8*nto, 1);
    resp2 = zeros(8*nto, 1);
    resp3 = zeros(8*nto, 1);
    for k = 0:nK - 1
        resp1 = resp1 + bufRec(rawPeaks1(1) + (1:8 * nto) + 8 * nto * k, 1) / 8;
        resp2 = resp2 + bufRec(rawPeaks2(1) + (1:8 * nto) + 8 * nto * k, 2) / 8;
        resp3 = resp3 + bufRec(rawPeaks3(1) + (1:8 * nto) + 8 * nto * k, 3) / 8;
    end
    resp1 = resp1 / nK;
    resp2 = resp2 / nK;
    resp3 = resp3 / nK;
    avResp = (resp1 + resp2 + resp3) / 3;
    resd1 = resp1 - avResp;
    resd2 = resp2 - avResp;
    resd3 = resp3 - avResp;
    respIdx = (1:nto)';
    totalVar = (std(resd1(respIdx) * sqrt(4)) ^ 2 ...
        + std(resd2(respIdx) * sqrt(4)) ^ 2 ...
        + std(resd3(respIdx) * sqrt(4)) ^ 2) / 2;
    %% ---- extended response ----
    rTmp = (bufRef(:, 1) - bufRef(:, 2)) / 32 + bufRef(:, 3) / 16;
    rawPeakRef = chkIdx(rTmp > 0.99 & chkIdx > safeStart & chkIdx + 8 * nto < safeEnd) - 44;
    rExt = zeros(8*nto, 1);
    rTmpResp = (bufRec(:, 1) - bufRec(:, 2)) / 32 + bufRec(:, 3) / 16;
    for k = 0:nK - 1
        rExt = rExt + rTmpResp(rawPeakRef(1) + (1:8 * nto) + 8 * nto * k);
    end
    rExt = rExt / nK;
    %% ---- calculate frequency response
    fftlv = 2 ^ ceil(log2(nto));
    fx = (0:fftlv - 1)' / fftlv * fs;
    fxEx = (0:fftlv*4 - 1)' / fftlv /4 * fs;
    avGain = abs(fft(avResp(1:nto), fftlv));
    extGain = abs(fft(rExt(1:4*nto), fftlv*4));
    exRrv = [rRV;rRV(end-1:-1:2);rRV];
    exRrv = exRrv(1:fftlv);
    bgSpec = abs(fft(exRrv));
    prepSpec = abs(fft(xRec(nto + (1:fftlv))));
    exResd1 = [resd1(respIdx);resd1(nto-1:-1:2);resd1] * sqrt(4);
    exResd2 = [resd2(respIdx);resd2(nto-1:-1:2);resd2] * sqrt(4);
    exResd3 = [resd3(respIdx);resd3(nto-1:-1:2);resd3] * sqrt(4);
    exResd1 = exResd1(1:fftlv);
    exResd2 = exResd2(1:fftlv);
    exResd3 = exResd3(1:fftlv);
    exResPwspec = (abs(fft(exResd1)) .^ 2 + abs(fft(exResd2)) .^ 2 ...
        + abs(fft(exResd3)) .^ 2) / 2;
    nlVar = (std(exResd1)^2 +  std(exResd2)^2 + std(exResd3)^2) / 2;
    % ---- smoothed response ----
    fxH = fx * 2 ^ (1/6);
    fxL = fx * 2 ^ (-1/6);
    fBw = fxH - fxL;
    fBw(1) = fBw(2);
    fxExH = fxEx * 2 ^ (1/6);
    fxExL = fxEx * 2 ^ (-1/6);
    fBwEx = fxExH - fxExL;
    fBwEx(1) = fBwEx(2);
    cumAvGain = cumsum(avGain .^2 * fx(2));
    cumExtGain = cumsum(extGain .^2 * fxEx(2));
    cumBGspec = cumsum(bgSpec .^2 * fx(2));
    cumPrepspec = cumsum(prepSpec .^2 * fx(2));
    cumResPwSpec = cumsum(exResPwspec * fx(2));
    smAvGain = (interp1(fx, cumAvGain, fxH, 'linear', 'extrap') ...
        - interp1(fx, cumAvGain, fxL, 'linear', 'extrap')) ./ fBw;
    smExtGain = (interp1(fxEx, cumExtGain, fxExH, 'linear', 'extrap') ...
        - interp1(fxEx, cumExtGain, fxExL, 'linear', 'extrap')) ./ fBwEx;
    smBGspec = (interp1(fx, cumBGspec, fxH, 'linear', 'extrap') ...
        - interp1(fx, cumBGspec, fxL, 'linear', 'extrap')) ./ fBw;
    smPrepspec = (interp1(fx, cumPrepspec, fxH, 'linear', 'extrap') ...
        - interp1(fx, cumPrepspec, fxL, 'linear', 'extrap')) ./ fBw;
    smResPwSpec = (interp1(fx, cumResPwSpec, fxH, 'linear', 'extrap') ...
        - interp1(fx, cumResPwSpec, fxL, 'linear', 'extrap')) ./ fBw;
    % ---- copy output -----
    output.xRec = xRec;
    output.xRef = xRef;
    output.avRespLevel = 20 * log10(std(avResp(1:nto)));
    output.bgLevel = 20 * log10(std(rRV));
    output.totalLevel = 10 * log10(totalVar);
    output.nlLevel = 10 * log10(nlVar);
    output.avResponse = avResp;
    output.extResponse = rExt;
    output.resd1 = resd1;
    output.resd2 = resd2;
    output.resd3 = resd3;
    output.bufRec = bufRec;
    output.bufRef = bufRef;
    output.frequencyAxis = fx;
    output.frequencyAxisEx = fxEx;
    output.avGain = 20 * log10(avGain);
    output.extGain = 20 * log10(extGain);
    output.backGroundSpec = 20 * log10(bgSpec);
    output.precedingSpec = 20 * log10(prepSpec);
    output.nonLinearSpec = 10 * log10(exResPwspec);
    output.smavGain = 10 * log10(smAvGain);
    output.smextGain = 10 * log10(smExtGain);
    output.smbackGroundSpec = 10 * log10(smBGspec);
    output.smprecedingSpec = 10 * log10(smPrepspec);
    output.smnonLinearSpec = 10 * log10(smResPwSpec);
    output.okInd = 1;
    output.elapsedTime = toc(startTic);
else
    output.okInd = 0;
end
end

