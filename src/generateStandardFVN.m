function output = generateStandardFVN(fs, sigmaT)
% Generate a unit FVN with standard parameter setting
%   output = generateStandardFVN(fs, sigmaT)
% Argument
%   fs       : sampling frequency (Hz)
%   sigmaT   : duration of the FVN (s)
%
% Output
%   output   : structure variable with the following field
%     xFVN          : FVN signal located in the center of DFT buffer
%     timeAxis      : time axis for xFVN (s)
%     frequencyAxis : frequency axis for phase for generating FVN (Hz)
%     mixedPhase    : phase for generating FVN (Hz)
%     frequencySegment : average frequency interval between centers (Hz)
%     normalizedBandwidth : nominal bandwidth (Hz)
%     fftl          : length of DFT (Discrete Fourier Transform) buffer (bins)
%     samplingFrequency : sampling frequency (Hz)
%
% Reference
%   [1] Kawahara, H.: Application of the velvet noise and its variant for 
%   synthetic speech and singing, SIGMUS Tech. Report. IPSJ, Vol. 118, 
%   No. 8 (2018).
%   [2] Kawahara, H., Sakakibara, K., Morise, M., Banno, H., Toda, T., 
%   Irino, T. (2017) A New Cosine Series Antialiasing Function and its 
%   Application to Aliasing-Free Glottal Source Models for Speech and 
%   Singing Synthesis. Proc. Interspeech 2017, 1358-1362, 
%   DOI: 10.21437/Interspeech.2017-15.

% copyright 2019 Hideki Kawahara
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

fftl = 2 ^ ceil(log2(fs * 15 * sigmaT));
output = generateFVN6real(fs, fftl, 2/(5 * sigmaT), 1/(5 * sigmaT), 2);
end
