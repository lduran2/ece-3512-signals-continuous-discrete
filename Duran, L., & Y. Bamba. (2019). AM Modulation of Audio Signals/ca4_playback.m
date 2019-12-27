%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTER ASSIGNMENT 4
%%%%%%%%%%%%%%%%%%%%%%%
% November 24, 2019
% ./ca4/playback.m
% Peforms bode graphs of the impedances of the phase 1s of the circuit,
% displays the values of the impedance denominators of these and the time
% constant, tau, of the phase 4, and plays back the demodulated audio.
% by: Duran, Leomar <https://github.com/lduran2>
% by: Yacouba Bamba
% on: 2019-11-24 T18:54ZR
% submitted for: Temple University, ECE 3512 Signals, Fall 2019 Semester

% clear the namespace
clear;

% metric prefixes
kilo = 1e+03;
m = 1e-03;
u = 1e-06;
n = 1e-09;

% frequency analysis limits
X_xlim = [-40 40; ...
           -5  5 ]*kilo; % [rad/s]
nX_xlim = size(X_xlim,1);
X_ylim = [0 2.5]*kilo;
xlim_names = { 'carriers', '0' }


%% Stage 1 (fc = 20 kHz)

% parameters
L = 3.3*m;%[H]
C1 = 19.2*n;%[F]
LC1 = L*C1;
s1 = (j*2*pi*20*kilo);

% check denominator
den1 = (((s1^2)*LC1) + 1)

% impedance
Z1 = tf([L 0],[LC1 0 1]);
figure(1);
bode(Z1);
title("Impedance of parallel LC filter at carrier frequency 20 kHz");

%% Stage 1 (fc = 30 kHz)

% parameters
L; % no change
C2 = 8.53*n;%[F]
LC2 = L*C2;
s2 = (j*2*pi*30*kilo);

% check denominator
den2 = (((s1^2)*LC1) + 1)

% impedance
Z2 = tf([L 0],[LC2 0 1]);
figure(2);
bode(Z2);
title("Impedance of parallel LC filter at carrier frequency 30 kHz");

%% Stage 4

% parameters
R = 1*kilo;%[Ohm]
C3 = 1*u;%[F]
RC = R*C3;

% echo RC time constant
tau = RC

%% play the demodulated audio

% file names
lvm_names  = split('bedbugs,broken', ',');
lvm_prefix = 'my-file-';
lvm_suffix = '.lvm';
% number of files
n_lvms = size(lvm_names, 1);

% for each lvm file
for k=1:n_lvms
    % concatenate file names
    fullname = [lvm_prefix, lvm_names{k}, lvm_suffix];
    
    % import each demodulated files
    [x Fs] = lvm_import(fullname);

    % display the frequency content
    [X f] = myFFT(x, Fs);
    figure(3);
    i_row = (nX_xlim*(k-1));
    % for each audio signal, show both with the frequency to carriers and
    % with magnitude to 0
    for m=1:nX_xlim
        subplot(n_lvms, nX_xlim, (i_row + m));
        plot(f,abs(X));
        title(['Demodulated signal ' lvm_names{k} ' to ' xlim_names{m}]);
        xlabel('frequency [Hz]');
        ylabel('magnitude');
        xlim(X_xlim(m,:));
        ylim(X_ylim);
    end %for m=1:nX_xlim

    % play each
    soundsc(x, Fs);
    % wait length of audio + one second
    length = (size(x, 1)/Fs);
    pause(length + 1);
end %for k=1:n_lvms

disp('Done.')
