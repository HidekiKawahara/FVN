function x = generateMissingFundamental(fc, id1, id2 )

fs = 44100;
nRepat = 14;
nto = round(0.5 * fs);
w = nuttallwin(4410 * 2);
tmp = load('fvnMin200ms.mat');
fvnMin200ms = tmp.fvnMin200ms;
fftl = size(fvnMin200ms, 1);
idx = 1:fftl;
buffer = zeros(fftl + nRepat * nto + 1, 1);
csign = 1;
for ii = 1:nRepat
    buffer(idx + (ii - 1) * nto) = fvnMin200ms(:, id1) ...
        + buffer(idx + (ii - 1) * nto);
    buffer(idx + (ii - 1) * nto) = fvnMin200ms(:, id2) * csign ...
        + buffer(idx + (ii - 1) * nto);
    csign = csign * (-1);
end
y = fftfilt(w, buffer);
rmsMod = std(y(round(fs*4):round(fs * 8)));
foBase = zeros(fftl + nRepat * nto + 1, 1) + fc;
foMod = 2 .^ ((log2(foBase) * 1200 + y / rmsMod * 25) / 1200);
x = sin(cumsum(2 * pi * foMod / fs)) * 0.2;
for ii = 2:10
    x = x + sin(ii * cumsum(2 * pi * foMod / fs)) / ii;
end
end

