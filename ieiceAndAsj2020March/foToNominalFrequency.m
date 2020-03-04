function nominalFrequency = foToNominalFrequency(fo)
% nominalFrequency = foToNominalFrequency(fo)
% This function quantizes fo into cromatic scale

NOTEA0 = 27.5;
HALFSEMITONENUDGE = 2 ^ (1 / 24);
octaveNo = floor(log2(fo / NOTEA0 * HALFSEMITONENUDGE));
noteNo = rem(round((log2(fo / NOTEA0)) * 12), 12) + 1;
nominalFrequency = 27.5 * 2 ^ ((noteNo - 1) / 12 + octaveNo);
end