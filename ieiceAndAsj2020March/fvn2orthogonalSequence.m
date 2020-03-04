function output = fvn2orthogonalSequence(fvnSet, nSequence, fvnIds, fs, tResponse, unitRepetition)

cOrthogonalBinary = zeros(unitRepetition * 2 ^ (nSequence - 1), nSequence);
lFvn = size(fvnSet, 1);
%tResponse = 0.5;
lResponse = round(fs * tResponse);
nRepetition = unitRepetition * 2 ^ (nSequence - 1);
lSequence = lFvn + lResponse * nRepetition;

cOrthogonalBinary(:, 1) = 1;
for iSequence = 2:nSequence
    tmpCoefficient = 1;
    for iPlase = 1:size(cOrthogonalBinary, 1)
        cOrthogonalBinary(iPlase, iSequence) = tmpCoefficient;
        if rem(iPlase, 2 ^ (iSequence - 2)) == 0
            tmpCoefficient = -tmpCoefficient;
        end
    end
end
fvnOrthogonalSequence = zeros(lSequence, nSequence);
basePlace = 0;
for iPlace = 1:nRepetition
    currentIndex = basePlace + (1:lFvn);
    for iSequence = 1:nSequence
        fvnId = fvnIds(iSequence);
        fvnOrthogonalSequence(currentIndex, iSequence) = ...
            fvnOrthogonalSequence(currentIndex, iSequence) ...
            + fvnSet(:, fvnId) * cOrthogonalBinary(iPlace, iSequence);
    end
    basePlace = basePlace + lResponse;
end
output.fvnOrthogonalSequence = fvnOrthogonalSequence;
fvnMixedSequence = sum(fvnOrthogonalSequence, 2);
output.fvnMixedSequence = fvnMixedSequence;
output.cOrthogonalBinary = cOrthogonalBinary;
output.fvnIds = fvnIds;
output.maxPeakValue = max(abs(fvnMixedSequence));
fvnMixSequenceBody = fvnMixedSequence(round(lFvn / 2) + lResponse *2 ...
    + (1:nSequence * 2 * lResponse));
output.lResponse = lResponse; 
output.kurtosisValue = mean(fvnMixSequenceBody .^ 4) / std(fvnMixSequenceBody) ^ 4;
output.nRepetition = nRepetition;
output.unitRepetition = unitRepetition;
output.usedFvns = fvnSet(:, fvnIds);
end