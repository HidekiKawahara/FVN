function scoreViewPosition = foToScorePosition(fo)
% convert fo to position on the musical score
scoreLocation = 0:32;
chromaticNoteId = [0 2 3 5 7 8 10  [0 2 3 5 7 8 10] + 12  ...
    [0 2 3 5 7 8 10] + 24 [0 2 3 5 7 8 10] + 36 [0 2 3 5 7] + 48];
scoreViewPosition = interp1(chromaticNoteId, scoreLocation, ...
    log2(fo / 55) * 12, "linear", "extrap");
end

