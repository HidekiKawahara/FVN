function output = generateFVN6real(fs, fftl, bw_c, fdiv_av, pdiv)
%  output = generateFVN6real(fs, fftl, bw_c, fdiv_av, pdiv)
%  This function generates one segment of FVN (Frequency domain Velvet
%  Noise) using the six-term cosine series of reference [2].
%
% Argument
%   fs      : sampling frequency (Hz)
%   fftl    : length of DFT (Discrete Fourier Transform) buffer (bins)
%   bw_c    : nominal bandwidth (Hz)
%   fdiv_av : average frequency interval between centers (Hz)
%   pdiv    : denominator to divide pi
%
% Output
%   output: structure with the following fields
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
%
% copyright 2018 Hideki Kawahara
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

fx = (0:fftl - 1)' / fftl * fs;

bw_n = bw_c / fx(2);
base_idx_nuttall = (-ceil(3 * bw_n):ceil(3 * bw_n));

buffer = zeros(fftl, 1);

for fc = 1:fdiv_av / fx(2):fftl / 2
  tmp = fc + (fdiv_av / fx(2)) * rand(1, 1);
  k = floor(tmp);
  fraction = tmp - k;
  idx = rem(k + base_idx_nuttall - 1 + fftl, fftl) + 1;
  w_sixterm = sixtermCosFrac(ceil(3 * bw_n) * 2 + 1, fraction - 0.5);
  buffer(idx) = buffer(idx) + pi / pdiv * w_sixterm * sign(rand(1,1) - 0.5);
end
txnn = ((1:fftl)' - fftl / 2 - 1) / fs;
rnew_buffer = -buffer([1 end:-1:2]);
final_buffer = buffer + rnew_buffer;
rMix = real(fftshift(ifft(exp(1i * final_buffer))));
output.xFVN = rMix;
output.timeAxis = txnn;
output.frequencyAxis = fx;
output.mixedPhase = final_buffer;
output.frequencySegment = fdiv_av;
output.normalizedBandwidth = bw_c;
output.fftl = fftl;
output.samplingFrequency = fs;
end

function output = sixtermCosFrac(n_length, fraction)
base_index = ((0:n_length - 1)' - n_length / 2 + 0.5 + fraction) / n_length;
BB = [0.2624710164;0.4265335164;0.2250165621;0.0726831633;0.0125124215;0.0007833203];
output = cos(base_index * [0 1 2 3 4 5] * 2 * pi) * BB;
end