function varargout = standingWvViewer(varargin)
% STANDINGWVVIEWER MATLAB code for standingWvViewer.fig
%      STANDINGWVVIEWER, by itself, creates a new STANDINGWVVIEWER or raises the existing
%      singleton*.
%
%      H = STANDINGWVVIEWER returns the handle to a new STANDINGWVVIEWER or the handle to
%      the existing singleton*.
%
%      STANDINGWVVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STANDINGWVVIEWER.M with the given input arguments.
%
%      STANDINGWVVIEWER('Property','Value',...) creates a new STANDINGWVVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before standingWvViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to standingWvViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help standingWvViewer

% Last Modified by GUIDE v2.5 18-Aug-2019 16:37:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @standingWvViewer_OpeningFcn, ...
    'gui_OutputFcn',  @standingWvViewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end


% --- Executes just before standingWvViewer is made visible.
function standingWvViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to standingWvViewer (see VARARGIN)

% Choose default command line output for standingWvViewer
handles.output = hObject;
handles = initializeGUI(handles);
handles.player1 = audioplayer(handles.testSignal, handles.setting.samplingFrequency, ...
    24);
handles.recorder1 = audiorecorder(handles.setting.samplingFrequency, 24, 2);
set(handles.recorder1, ...
    'TimerFcn', @update_display, 'TimerPeriod', 0.2, 'userdata', handles);

% Update handles structure
guidata(hObject, handles);
play(handles.player1);
record(handles.recorder1);
set(handles.startButton, 'enable', 'off');
set(handles.reportButton, 'enable', 'off');
set(handles.stopButton, 'enable', 'on');
set(handles.radiobutton70, 'enable', 'on', 'value', 0);
set(handles.radiobutton75, 'enable', 'on', 'value', 0);
set(handles.radiobutton80, 'enable', 'on', 'value', 0);

% UIWAIT makes standingWvViewer wait for user response (see UIRESUME)
% uiwait(handles.RoomResponseGUI);
end


% --- Outputs from this function are returned to the command line.
function varargout = standingWvViewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

%--- Private functions
function handles = initializeGUI(handles)
setting = struct;
% --- Initialize FVN
tmp = load(which('fvnMin100ms.mat'));
fs = tmp.fs;
BUFFER_DURATION = 100; % Audio input buffer length in s
RESPONSE_LENGTH = 300; % Assumed effective response length in ms
nResponse = round(RESPONSE_LENGTH / 1000 * fs);
xFVN1 = tmp.fvnMin100ms(:, 1);
xFVN2 = tmp.fvnMin100ms(:, 2);
testSignal = zeros(fs * BUFFER_DURATION, 2);
lengthFvn = length(xFVN1);
startIndex = 0;
baseIndex = 1:lengthFvn;
while startIndex + lengthFvn + nResponse < fs * BUFFER_DURATION
    testSignal(startIndex + baseIndex, 1) = ...
        testSignal(startIndex + baseIndex, 1) + xFVN1;
    testSignal(startIndex + baseIndex, 2) = ...
        testSignal(startIndex + baseIndex, 2) + xFVN2;
    xFVN2 = -xFVN2;
    startIndex = startIndex + nResponse;
end
%%  Spectral shaping using pink noise shape

LOWER_LIMIT = 40;
FFT_LENGTH = 32768;
fx = (0:FFT_LENGTH - 1)' / FFT_LENGTH * fs;
pinkgain = 1.0 ./ sqrt(fx);
pinkgain(fx < LOWER_LIMIT) = max(pinkgain(fx >= LOWER_LIMIT));
pinkgain(FFT_LENGTH:-1:FFT_LENGTH / 2 + 1) = pinkgain(2:FFT_LENGTH / 2 + 1);
pinkgain = pinkgain / pinkgain(1);
% --- all pole approximation
rawAutoCorr = real(ifft(pinkgain .^ 2));
SHAPER_ORDER = 45;
[a1, ~] = levinson(rawAutoCorr, SHAPER_ORDER);

testSignal = filter(1, a1, testSignal); % make the signal pink
testSignal = testSignal / max(abs(testSignal(:))) * 0.85;
sumSignal = testSignal(:, 1) + testSignal(:, 2);
tt = ((1:lengthFvn) - lengthFvn / 2) / fs;
fvnFilt1 = xFVN1(abs(tt) < 0.5);
fvnFilt2 = xFVN2(abs(tt) < 0.5);
fvnFilt1 = fvnFilt1(end:-1:1);
fvnFilt2 = fvnFilt2(end:-1:1);
fx = (0:lengthFvn - 1)' / lengthFvn * fs;
% ---
axes(handles.spectrumAxis);
segIdx = (1:3 * fs);
x1 = (sumSignal(3 * fs + segIdx) + sumSignal(3 * fs + segIdx - nResponse)) / 2;
x2 = (sumSignal(3 * fs + segIdx) - sumSignal(3 * fs + segIdx - nResponse)) / 2;
x3 = ((sumSignal(3 * fs + segIdx) + sumSignal(3 * fs + segIdx - nResponse)) / 2 - ...
    (sumSignal(3 * fs + segIdx - 2 * nResponse) + ...
    sumSignal(3 * fs + segIdx - nResponse)) / 2) / 2;
w = blackman(3 * fs);
w = w / sum(w);
cumpw1 = cumsum(abs(fft(x1 .* w, lengthFvn)) .^ 2);
cumpw2 = cumsum(abs(fft(x2 .* w, lengthFvn)) .^ 2);
cumpw3 = cumsum(abs(fft(x3 .* w, lengthFvn)) .^ 2);
fo = fs / nResponse;
fxHarmonics = (fo:fo:350)';
pw1Harmonics = interp1(fx, cumpw1, fxHarmonics + fo / 2) - interp1(fx, cumpw1, fxHarmonics - fo / 2);
pw2Harmonics = interp1(fx, cumpw2, fxHarmonics + fo / 2) - interp1(fx, cumpw2, fxHarmonics - fo / 2);
pw3Harmonics = interp1(fx, cumpw3, fxHarmonics + fo / 2) - interp1(fx, cumpw3, fxHarmonics - fo / 2);
handles.lchSpectrumHandle = plot(fxHarmonics, 10 * log10(pw1Harmonics), 'linewidth', 2);
hold all
handles.rchSpectrumHandle = plot(fxHarmonics - fo / 2, 10 * log10(pw2Harmonics), 'linewidth', 2);
handles.errSpectrumHandle = plot(fxHarmonics, 10 * log10(pw3Harmonics), 'linewidth', 2);
grid on;
%max(abs(x3))
set(gca, 'xlim', [50 350]);
handles.fx_harmonics = fxHarmonics;
handles.fo = fo;
handles.fx = fx;
handles.w = w;
handles.lengthFvn = lengthFvn;
handles.nResponse = nResponse;
xlabel('frequency (Hz)');
handles.splylabel = ylabel('level (dB, rel.)');

% --- Impulse response plot
axes(handles.iresoponseAxis);
l_resp = zeros(round(nResponse * 0.85), 1);
r_resp = zeros(round(nResponse * 0.85), 1);
ref_resp = zeros(round(nResponse * 0.85), 1);
t_resp = (1:length(l_resp))' / fs * 1000 - 2; % ms with 2 ms offset
zero_idx = round(2 / 1000 * fs);
l_resp(zero_idx + 22) = 1;
r_resp(zero_idx + 10) = 1;
ref_resp(zero_idx) = 1;
handles.lch_iresp_handle = plot(t_resp, l_resp, 'linewidth', 2);
hold all
handles.rch_iresp_handle = plot(t_resp, r_resp, 'linewidth', 2);
handles.refch_iresp_handle = plot(t_resp, ref_resp * NaN, 'linewidth', 1);
grid on;
axis([0 14 -1 1]);
handles.resp_length = length(l_resp);
handles.zero_idx = zero_idx;
xlabel('time (ms)');

% --- long impulse response plot
axes(handles.longRespAxis);
l_respL = zeros(round(nResponse * 1.2), 1);
r_respL = zeros(round(nResponse * 1.2), 1);
ref_respL = zeros(round(nResponse * 1.2), 1);
t_respL = (1:length(l_respL))' / fs * 1000 - 2; % ms with 2 ms offset
zero_idx = round(2 / 1000 * fs);
l_respL(zero_idx + 22) = 1;
r_respL(zero_idx + 10) = 1;
ref_respL(zero_idx) = 1;
handles.back_iresp_handleL = plot(t_respL, 20 * log10(l_respL), 'linewidth', 1);
handles.lch_color = get(handles.back_iresp_handleL, 'color');
hold all
handles.front_iresp_handleL = plot(t_respL, 20 * log10(r_respL), 'linewidth', 1);
handles.rch_color = get(handles.front_iresp_handleL, 'color');
handles.refch_iresp_handleL = plot(t_respL, 20 * log10(ref_respL), 'linewidth', 1);
grid on;
axis([-2 nResponse * 1.1 / fs * 1000 -60 0]);
handles.resp_lengthL = length(l_respL);
handles.zero_idx = zero_idx;
xlabel('time (ms)');

% --- log-linear spectrum
axes(handles.logSpecAxis);
fftlresp = 2 ^ ceil(log2(length(l_resp)));
fx_resp = (0:fftlresp - 1) / fftlresp * fs;
handles.lch_logspec = semilogx(fx_resp, 20 * log10(abs(fft(l_resp,fftlresp))), 'linewidth', 2);
hold all
handles.rch_logspec = semilogx(fx_resp, -5 + 20 * log10(abs(fft(r_resp,fftlresp))), 'linewidth', 2);
handles.err_logspec = semilogx(fx_resp, -10 + 20 * log10(abs(fft(ref_resp,fftlresp))), 'linewidth', 2);
grid on;
set(handles.logSpecAxis, 'xlim', [20 fs / 2], 'ylim', [-60 10]);
xlabel('frequency (Hz)');
ylabel('level (dB, rel.)');
handles.fftlresp = fftlresp;
handles.fx_resp = fx_resp;

% --- level indicator
axes(handles.levelAxis);
handles.l_bar = bar(1, -20 + 100, 'facecolor', [0 0.8 0]);
hold all
handles.r_bar = bar(2, -20 + 100, 'facecolor', [0 0.8 0]);
handles.l_peak = plot([0.5 1.5], [-10 -10] + 100, 'linewidth', 3, 'color', [0.8 0 0]);
handles.r_peak = plot([1.5 2.5], [-11 -11] + 100, 'linewidth', 3, 'color', [0.8 0 0]);
grid on;
axis([0.5 2.5 0 100])
set(handles.levelAxis, 'ytick', 0:10:100, 'ytickLabel', num2str((0:10:100)' - 100));
set(handles.levelAxis, 'xtick', [1 2], 'xticklabel', ['Mic'; 'Ref']);

% --- setting
setting.samplingFrequency = fs;
setting.signalDuration = BUFFER_DURATION;
setting.responseLengthInSample = nResponse;
setting.xFVN1 = xFVN1;
setting.xFVN1 = tmp.fvnMin100ms(:, 2);
setting.fvnFilt1 = fvnFilt1;
setting.fvnFilt2 = fvnFilt2;
setting.pinkCoeff = a1;
handles.testSignal = testSignal;
handles.setting = setting;
handles.calLevel = 0;
handles.levelBias = 0;
end

%--- private timer functions
function update_display(hObject, eventdata, handlesDummy)
handles = get(hObject, 'userdata');
y = getaudiodata(hObject);
fs = handles.setting.samplingFrequency;
lengthFvn = handles.lengthFvn;
w = handles.w;
fo = handles.fo;
fx_harmonics = handles.fx_harmonics;
fx = handles.fx;
nResponse = handles.nResponse;
a1 = handles.setting.pinkCoeff;
if length(y) > 2 * nResponse + 1
    set(handles.l_bar, 'ydata', 100 + 20 * log10(std(y(end-2*nResponse:end, 1))));
    set(handles.l_peak, 'ydata', [0 0] + 100 + 20 * log10(max(abs(y(end-2*nResponse:end, 1)))));
    set(handles.r_bar, 'ydata', 100 + 20 * log10(std(y(end-2*nResponse:end, 2))));
    set(handles.r_peak, 'ydata', [0 0] + 100 + 20 * log10(max(abs(y(end-2*nResponse:end, 2)))));
end
if length(y) > 7 * fs + 6 * nResponse %length(y) > 7 * fs + 4 * nResponse
    marginy = 7 * fs;
    sum_signal = y(end - marginy:end, :);
    sum_signal1 = fftfilt(handles.setting.fvnFilt1, sum_signal(length(sum_signal) - fs * 6:length(sum_signal), 1));
    sum_signal2 = fftfilt(handles.setting.fvnFilt2, sum_signal(length(sum_signal) - fs * 6:length(sum_signal), 1));
    sum_signal_ref = fftfilt(handles.setting.fvnFilt2, sum_signal(length(sum_signal) - fs * 6:length(sum_signal), 2));
    segIdx = (1:6 * nResponse);
    w = blackman(length(segIdx));
    w = w / sum(w);
    start_idx = length(sum_signal1) - length(segIdx) - fs;
    x1 = (sum_signal1(start_idx + segIdx) + sum_signal1(start_idx + segIdx - nResponse)) / 2;
    x2 = (sum_signal2(start_idx + segIdx) - sum_signal2(start_idx + segIdx - nResponse)) / 2;
    x3 = ((sum_signal1(start_idx + segIdx) + sum_signal1(start_idx + segIdx - nResponse)) / 2 ...
        - (sum_signal1(start_idx + segIdx - 2 * nResponse) + sum_signal1(start_idx + segIdx - nResponse)) / 2) / 2;
    x_ref = sum_signal_ref(start_idx + segIdx);
    cumpw1 = cumsum(abs(fft(x1 .* w, lengthFvn)) .^ 2);
    cumpw2 = cumsum(abs(fft(x2 .* w, lengthFvn)) .^ 2);
    cumpw3 = cumsum(abs(fft(x3 .* w, lengthFvn)) .^ 2);
    pw1_harmonics = interp1(fx, cumpw1, fx_harmonics + fo / 2) - interp1(fx, cumpw1, fx_harmonics - fo / 2);
    pw2_harmonics = interp1(fx, cumpw2, fx_harmonics + fo / 2 - fo / 2) - interp1(fx, cumpw2, fx_harmonics - fo / 2 - fo / 2);
    pw3_harmonics = interp1(fx, cumpw3, fx_harmonics + fo / 2) - interp1(fx, cumpw3, fx_harmonics - fo / 2);
    set(handles.lchSpectrumHandle, 'ydata', 10 * log10(pw1_harmonics) + handles.levelBias);
    set(handles.rchSpectrumHandle, 'ydata', 10 * log10(pw2_harmonics) + handles.levelBias);
    set(handles.errSpectrumHandle, 'ydata', 10 * log10(pw3_harmonics) + handles.levelBias);
    set(handles.spectrumAxis, 'ylim', [-120 -20] + handles.levelBias);
    measuredData.recordedWave = sum_signal;
    measuredData.x1Org = x1;
    measuredData.x2Org = x2;
    measuredData.x3Org = x3;

    
    % --- impulse response display
    x_ref = filter(a1 / sum(a1), 1, x_ref);
    x1 = filter(a1 / sum(a1), 1, x1);
    x2 = filter(a1 / sum(a1), 1, x2);
    x3 = filter(a1 / sum(a1), 1, x3);
    resp_length = handles.resp_length;
    resp_lengthL = handles.resp_lengthL;
    %resp_length = round(27 / 1000 * fs);
    zero_idx = handles.zero_idx;
    
    maxRef = max((x_ref));
    peak_loc = segIdx((x_ref) > 0.7 * maxRef);
    ref_loc = max(peak_loc(peak_loc < 3 * nResponse - zero_idx - resp_length));
    %{
    peak_loc = ref_idx(abs(x_ref) > 0.7 * maxRef & (ref_idx < length(x_ref) - resp_length - zero_idx));
    ref_loc = max(peak_loc);
    set(handles.lch_iresp_handle, 'ydata', x1(ref_loc - zero_idx + (1:resp_length)) / max(abs(x1(ref_loc - zero_idx + (1:resp_length)))));
    %}
    max_lch = max(abs(x1(ref_loc - zero_idx + (1:resp_length))));
    max_rch = max(abs(x2(ref_loc - zero_idx + (1:resp_length))));
    set(handles.lch_iresp_handle, 'ydata', x1(ref_loc - zero_idx + (1:resp_length)) / max_lch);
    set(handles.rch_iresp_handle, 'ydata', x2(ref_loc - zero_idx + (1:resp_length)) / max_rch);
    
    %---- save the latest data
    measuredData.x1 = x1;
    measuredData.x2 = x2;
    measuredData.x3 = x3;
    measuredData.x_ref = x_ref;
    measuredData.reference_location = ref_loc;
    measuredData.zero_idx = zero_idx;
    measuredData.calLevel = handles.calLevel;
    measuredData.levelBias = handles.levelBias;
    measuredData.samplingFrequency = fs;
    handles.measuredData = measuredData;
    
    % --- log spectrum display
    fftlresp = handles.fftlresp;
    fx_resp = handles.fx_resp;
    cummpw1 = cumsum(abs(fft(x1(ref_loc + (1:resp_length)), fftlresp)) .^ 2);
    cummpw2 = cumsum(abs(fft(x2(ref_loc + (1:resp_length)), fftlresp)) .^ 2);
    cummpw3 = cumsum(abs(fft(x3(ref_loc + (1:resp_length)), fftlresp)) .^ 2);
    fx_resp_H = fx_resp * 2 ^ (1/6);
    fx_resp_L = fx_resp * 2 ^ (-1/6);
    pws1 = abs(interp1(fx_resp, cummpw1, fx_resp_H, 'linear', 'extrap') ...
        - interp1(fx_resp, cummpw1, fx_resp_L, 'linear', 'extrap')) ./ (fx_resp_H - fx_resp_L);
    pws2 = abs(interp1(fx_resp, cummpw2, fx_resp_H, 'linear', 'extrap') ...
        - interp1(fx_resp, cummpw2, fx_resp_L, 'linear', 'extrap')) ./ (fx_resp_H - fx_resp_L);
    pws3 = abs(interp1(fx_resp, cummpw3, fx_resp_H, 'linear', 'extrap') ...
        - interp1(fx_resp, cummpw3, fx_resp_L, 'linear', 'extrap')) ./ (fx_resp_H - fx_resp_L);
    selectr = 100 * 2 .^ (0:1/24:log2(10000 / 100));
    selectrIdx = round(selectr / fx_resp(2));
    maxPwr = max(mean(pws1(selectrIdx)), mean(pws2(selectrIdx)));
    set(handles.lch_logspec, 'ydata', 10 * log10(pws1 / maxPwr));
    set(handles.rch_logspec, 'ydata', 10 * log10(pws2 / maxPwr));
    set(handles.err_logspec, 'ydata', 10 * log10(pws3 / maxPwr));
    
    % --- long impulse response
    %{
    if get(handles.RchannelRadio, 'value')
        set(handles.front_iresp_handleL, 'visible', 'on');
    else
        set(handles.front_iresp_handleL, 'visible', 'off');
    end
    if get(handles.LchannelRadio, 'value')
        set(handles.back_iresp_handleL, 'visible', 'on');
    else
        set(handles.back_iresp_handleL, 'visible', 'off');
    end
    %}
    lx1_max = length(x1);
    if mean(pws1(selectrIdx)) < mean(pws2(selectrIdx))
    set(handles.back_iresp_handleL, 'ydata', 20*log10(abs(x1(min(lx1_max, ref_loc - zero_idx + (1:resp_lengthL))) / max_lch)), ...
        'color', handles.lch_color);
    set(handles.front_iresp_handleL, 'ydata', 20*log10(abs(x2(min(lx1_max, ref_loc - zero_idx + (1:resp_lengthL))) / max_rch)), ...
        'color', handles.rch_color);
    else
    set(handles.back_iresp_handleL, 'ydata', 20*log10(abs(x2(min(lx1_max, ref_loc - zero_idx + (1:resp_lengthL))) / max_rch)), ...
        'color', handles.rch_color);
    set(handles.front_iresp_handleL, 'ydata', 20*log10(abs(x1(min(lx1_max, ref_loc - zero_idx + (1:resp_lengthL))) / max_lch)), ...
        'color', handles.lch_color);
    end

    %end
    if handles.calLevel == 0
        if get(handles.radiobutton80, 'value');handles.calLevel = 80; end
        if get(handles.radiobutton75, 'value');handles.calLevel = 75; end
        if get(handles.radiobutton70, 'value');handles.calLevel = 70; end
        set(handles.SPLText, 'string', num2str(handles.calLevel ,'%4.0f'));
        if handles.calLevel ~= 0
            sploutput = loudnessWithA(y(end - marginy:end, 1), fs);
            handles.levelBias = handles.calLevel - sploutput.slow(end);
            set(handles.splylabel, 'string', 'sound pressure level (dB)');
        end
    else
        set(handles.SPLbuttonGroup, 'visible', 'off');
    end
end % if length(y) > 7 * fs + 4 * nResponse
% --- time display
set(handles.soundText, 'string', num2str(handles.setting.signalDuration - length(y) / fs , ...
    '%5.1f'));

if length(y) / fs > handles.setting.signalDuration
    switch get(handles.recorder1, 'running')
        case 'on'
            stop(handles.recorder1);
    end
    switch get(handles.player1, 'running')
        case 'on'
            stop(handles.player1);
    end
    set(handles.startButton, 'enable', 'on');
    set(handles.stopButton, 'enable', 'off');
end
set(hObject, 'userdata', handles);
end

% --- Executes on button press in startButton.
function startButton_Callback(hObject, eventdata, handles)
% hObject    handle to startButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(handles.recorder1, 'running')
    case 'off'
        record(handles.recorder1);
end
switch get(handles.player1, 'running')
    case 'off'
        play(handles.player1);
end
set(handles.startButton, 'enable', 'off');
set(handles.stopButton, 'enable', 'on');
set(handles.reportButton, 'enable', 'off');
guidata(hObject, handles);
end


% --- Executes on button press in stopButton.
function stopButton_Callback(hObject, eventdata, handles)
% hObject    handle to stopButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(handles.recorder1, 'running')
    case 'on'
        stop(handles.recorder1);
end
switch get(handles.player1, 'running')
    case 'on'
        stop(handles.player1);
end
set(handles.startButton, 'enable', 'on');
set(handles.stopButton, 'enable', 'off');
set(handles.reportButton, 'enable', 'on');
guidata(hObject, handles);
end


% --- Executes on button press in quitButton.
function quitButton_Callback(hObject, eventdata, handles)
% hObject    handle to quitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(handles.recorder1, 'running')
    case 'on'
        stop(handles.recorder1);
end
switch get(handles.player1, 'running')
    case 'on'
        stop(handles.player1);
end
close(handles.RoomResponseGUI)
end


% --------------------------------------------------------------------
function SPLbuttonGroup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to SPLbuttonGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --- Executes on button press in LchannelRadio.
function LchannelRadio_Callback(hObject, eventdata, handles)
% hObject    handle to LchannelRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LchannelRadio
end


% --- Executes on button press in RchannelRadio.
function RchannelRadio_Callback(hObject, eventdata, handles)
% hObject    handle to RchannelRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RchannelRadio
end


% --- Executes on button press in reportButton.
function reportButton_Callback(hObject, eventdata, handles)
% hObject    handle to reportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
reportFigName = ['stView' datestr(now, 30) '.eps'];
print(handles.RoomResponseGUI, '-noui', '-depsc', reportFigName);
reportDataName = ['stData' datestr(now, 30) '.mat'];
userdata = get(handles.recorder1, 'userdata');
measuredData = userdata.measuredData;
save(reportDataName, 'measuredData');
end