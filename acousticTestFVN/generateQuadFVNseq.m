function output = generateQuadFVNseq(fvnSet, nRepeat, nto)
%  generate orthogonal FVN sequences using four unit FVNs
%
%  Input argument
%    fvnSet  : matrix with four unit FVNs as its column
%    nRepeat : number of repetitions of the unit FVN in eaxh sequence
%    nto     : allocation interval of unit FVNs (samples)
%
%  Output
%    output  : structure variable with the following fields
%      xTest       : mixture of all FVN sequences
%      xTestR      : mixture of the first three FVN sequences
%      elapsedTime : elapsed time (s)

%  For details, please refer the following article
%
%   Kawahara, H., Sakakibara, K. I., Mizumachi, M., Morise, M., & Banno, H. (2020). 
%   Simultaneous measurement of time-invariant linear and nonlinear, 
%   and random and extra responses using frequency domain variant of velvet noise. 
%   arXiv preprint arXiv:2008.02439.

%   Copyright 2020 Hideki Kawahara
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

startTic = tic;
B4 = [1  1  1  1  1  1  1  1; ...
    1 -1  1 -1  1 -1  1 -1; ...
    1  1 -1 -1  1  1 -1 -1; ...
    1  1  1  1 -1 -1 -1 -1];
fftl = size(fvnSet, 1);
idx = 1:fftl;
buffer= zeros(fftl + nRepeat * nto, 1);
buffer1 = zeros(fftl + nRepeat * nto, 1);

for ii = 1:nRepeat
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

