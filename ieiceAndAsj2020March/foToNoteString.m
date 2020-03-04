function noteString = foToNoteString(fo)
% noteString = foToNoteString(fo)
% This function convertes a frequency to a corresponding musical note name
% example 440 Hz --> A4
%    fo = 440;
%    noteString = foToNoteString(fo);
%  then:
%    contents of the 'noteString' is a string 'A4'
%%
baseNames = {'A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#'};
NOTEC0 = 27.5 / 2 * 2 ^ (3/12);
NOTEA0 = 27.5;
HALFSEMITONENUDGE = 2 ^ (1 / 24);
octaveNo = floor(log2(fo / NOTEC0 * HALFSEMITONENUDGE));
noteNo = rem(round((log2(fo / NOTEA0)) * 12), 12) + 1;
noteString = [baseNames{noteNo} num2str(octaveNo)];
end