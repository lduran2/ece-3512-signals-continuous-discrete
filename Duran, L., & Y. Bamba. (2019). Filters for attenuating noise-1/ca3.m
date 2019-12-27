%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTER ASSIGNMENT 3
%%%%%%%%%%%%%%%%%%%%%%%
% October 24, 2019
% ./ca3.m
% Attenuates a 60 Hz noise signal from an electrocardiogram.
% by: Duran, Leomar <https://github.com/lduran2>
% by: Yacouba Bamba
% on: 2019-10-25 T16:28ZQ
% submitted for: Temple University, ECE 3512 Signals, Fall 2019 Semester

% Clear everything and go to figure 1.
clf;       % clear current figure
figure(1); % go to figure 1
clf;       % clear figure 1
clear;     % clear the namespace
clc;       % clear the command window

% Load the ECG data...
data = load('ecg_data.txt'); % from the file.
t = data(:,1);               % time data
x = data(:,2);               % voltage data
Fs = 1000;%[samples/s]       % sampling rate

% Noise parameters
noise = 60;%[Hz]             % noise in Hertz
wnoise = 2*pi*noise;%[rad/s] % noise in rad/s
tolerance = 0.10*wnoise;     % tolerance = 10% of noise


%%
% Step 0: Plot the ECG signal with respect to time. You should see a very
% noise corrupted signal.
subplot(3,2,1);                     % subplot 1 in 3x2 grid
plot(t, x);                         % plot input voltage vs time
title('time domain of ECG signal'); % title the plot
xlabel('time [s]');                 % label x ...
ylabel('magnitude [V]');            % and y axes


%%
% Step 1: Examine the frequency content of the signal to confirm that the
% noise source is in fact
% 60Hz hum.

% use myFFT to get X in frequency domain
[X f] = myFFT(x,Fs);      % = myFFT(x,fs,varargin)
X = transpose(X);         % X is oriented opposite the other vectors
subplot(3,2,2);           % subplot 2 in 3x2 grid
plot(f, abs(X));          % plot the magnitude versus frequency
title('frequency domain of ECG signal'); % title the plot
xlim([-100, 100]);        % limit to the domain from -100 Hz to 100 Hz
xlabel('frequency [Hz]'); % label the x ...
ylabel('magnitude [V]');  % and y axes


%%
% Step 2: Using the steps outlined in Chapter 4 and in lecture, design a
% filter in Matlab that will remove as much of the 60Hz noise as possible
% while keeping as much of the EEG signal as possible. You are limited to
% a third order filter. You may use commands such as butter or cheby1 to
% design your filter, and the lsim command to filter your ECG signal.

% Create the transfer function.
% H(s) = (sL + 1/sC)/(R + sL + 1/sC)
%      = (s^2 + 1/LC)/(s^2 + (R/L)s + 1/LC)
%      = (s^2 + wc^2)/(s^2 + rs + wc^2)
% wc := frequency to attenuate = sqrt(1/LC)
% r := tolerance = (R/L) = 10% of wc
H = @(s) (s.^2 + (wnoise)^2)./(s.^2 + s.*tolerance + (wnoise)^2);
Harr = H(j*2*pi*f); % the transfer function as a parallel array to
                    % frequency

% R = 0.5655 mOhm metal strip PR series chip resistor
%   <https://www.digikey.com/en/product-highlight/y/yageo/low-ohmic-chip-resistors>
% L ~ 15 uH fixed inductor
%   <https://www.mouser.com/ProductDetail/Murata-Electronics/1264EY-150M%3dP3?qs=sGAEpiMZZMsg%252By3WlYCkU588IWvyRtFzJRVoYaqW0EU%3D>
% C ~ 0.47 F double layer radial capacitor
%   <https://www.newark.com/panasonic/eecf5r5h474n/cap-0-47f-5-5v-double-layer-radial/dp/85Y4670>

% Plot the transfer function.
subplot(3,2,[3 4]);        % subplot 3 and 4 of a 3x2 grid (1x2 subplot)
plot(f, abs(Harr));        % plot the filter response versus frequency
title('filter response vs frequency'); % title the plot
xlim([-100, 100]);         % limit the domain between -100 Hz and 100 Hz
xlabel('frequency [Hz]');  % label the x
ylabel('filter response [V^0]'); % and y axes


%%
% Step 3: Apply the transfer function to find the output in frequency
% domain.
Y = Harr.*X; % This just involves multiplying the input by the transfer
             % function like a gain.

% Plot the output in frequency domain.
subplot(3,2,6);           % subplot 6 of a 3x2 grid (skip 5 for time domain)
plot(f, abs(Y));          % plot the magnitude of output vs frequency
title('frequency domain of response');
xlim([-100, 100]);        % limit to domain of -100 Hz to 100 Hz
xlabel('frequency [Hz]'); % label the x ...
ylabel('magnitude [V]');  % and y axes


%%
% Step 4: Create the time domain by adding the harmonics for each frequency.

% Add together all the cosines represented by the output Y and the
% frequencies f.
y = 0; % y is only resized once, so I didn't bother with allocation
for k = 1:size(f,2) % for all frequencies
    y = y + 2*abs(Y(k))* cos(2*pi*f(k)*t + angle(Y(k)));
end % for k = 1:size(f,2);

% Plot the time domain.
subplot(3,2,5);           % subplot 5 of the 3x2 grid
plot(t, y);               % plot output versus the time domain
title('time domain of response'); % title the plot
xlabel('frequency [Hz]'); % label the x ...
ylabel('magnitude [V]');  % and y axes

disp('Done');
