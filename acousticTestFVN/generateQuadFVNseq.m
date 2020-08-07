function output = generateQuadFVNseq(fvnSet, nRepat, nto)
startTic = tic;
B4 = [1  1  1  1  1  1  1  1; ...
    1 -1  1 -1  1 -1  1 -1; ...
    1  1 -1 -1  1  1 -1 -1; ...
    1  1  1  1 -1 -1 -1 -1];
fftl = size(fvnSet, 1);
idx = 1:fftl;
buffer= zeros(fftl + nRepat * nto, 1);
buffer1 = zeros(fftl + nRepat * nto, 1);

for ii = 1:nRepat
    for jj = 1:4
        kk = rem((ii - 1), 8) + 1;
        buffer(idx + (ii - 1) * nto) = fvnSet(:, jj) * B4(jj, kk)...
            + buffer(idx + (ii - 1) * nto);
        if jj < 4
            buffer1(idx + (ii - 1) * nto) = fvnSet(:, jj) * B4(jj, kk)...
                + buffer1(idx + (ii - 1) * nto);
        end
    end
end
output.xTest = buffer;
output.xTestR = buffer1;
output.elapsedTime = toc(startTic);
end

