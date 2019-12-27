%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTER ASSIGNMENT 4
%%%%%%%%%%%%%%%%%%%%%%%
% November 23, 2019
% ./ca4/modulate_2.m
% Bandlimits 2 signals to 2kHz, resamples them to 120x10^3 samples/second,
% and modulates them 10k apart.
% by: Duran, Leomar <https://github.com/lduran2>
% by: Yacouba Bamba
% on: 2019-11-23 T21:51ZR
% submitted for: Temple University, ECE 3512 Signals, Fall 2019 Semester

% clear the namespace
clear;

% parameters of the audio
Fs = 44100;%[samples/s] % sampling frequency
duration = 5;%[s]       % duration of each audio
samples = Fs*duration;  % the total number of samples
srange = [1 samples];   % the range of samples to read

% parameters of the filter
order = 5;%[th]         % order Butter worth filter
wc = 2*pi*2000;%[rad/s] % cutoff frequency

% parameters for upsample
Fr = 120e+03;%[samples/s] % resampling frequency
resamples = Fr*duration;  % the total number of resamples

% parameters for the carrier signals
w_carry1 = 2*pi*20e+03; %[rad/s]
w_carry2 = 2*pi*30e+03; %[rad/s]

% frequency analysis limits
PROD_xlim = [10e+03 40e+03]; % [rad/s]
PROD_ylim = [0 1.6e-03];

%%
% Step 1) Get the audio signal and normalize for use.

% import the audio, ignore fs
[x1] = audioread('5s-mono-bass-Bedbugs.mp3', srange);
[x2] = audioread('5s-mono-bass-Broken.mp3', srange);

% take only channel 1 and normalize the input signals
norm_x1 = myNormalize(x1(:,1));
norm_x2 = myNormalize(x2(:,1));

% time signal, using original Fs
ts = (0 : (samples - 1))/Fs; % scale by 1:Fs.

%%
% Step 2) Bandlimit the two audio signals to 2kHz with a 5th order
% Butterworth filter.

% create the filter
[num den] = butter(5, wc, 's');
H = tf(num, den);
% perform bandlimits
y1 = lsim(H, norm_x1, ts);
y2 = lsim(H, norm_x2, ts);

%%
% Step 3 and 4) Normalize the signals and add 1, to get range [0, 2].
norm_y1 = 1 + myNormalize(y1);
norm_y2 = 1 + myNormalize(y2);

%%
% Step 5) Upsample the audio signals to 120kHz.
up_y1 = resample(norm_y1, Fr, Fs);
up_y2 = resample(norm_y2, Fr, Fs);

%%
% Step 6) Create the carrier signals.

% time signal, using resample Fr
tr = (0 : (resamples - 1))/Fr; % scale by 1:Fr.

% the carrier signals
carrier1 = cos(w_carry1 * tr);
carrier2 = cos(w_carry2 * tr);

%%
% Step 7) Dot multiply each upsample with a carrier
product = 0;
product = product + (up_y1 .* (carrier1'));
product = product + (up_y2 .* (carrier2'));

%%
% Step 8) Analysis of the frequency content.
[PROD f] = myFFT(product, Fr);
figure(1);
clf;
plot(f, abs(PROD));
title('Dot product in frequency domain');
xlabel('frequency [Hz]');
ylabel('magnitude');
% limit the range of plotting
xlim(PROD_xlim);
ylim(PROD_ylim);

%%
% Step 9) Export the product.
data = [ tr' product];
save('my-file.txt', 'data', '-ascii');


%%
% Normalizes the vector so that it's within the range [-1,1]
function result = myNormalize(v)
    result = v/max(abs(v));
end
