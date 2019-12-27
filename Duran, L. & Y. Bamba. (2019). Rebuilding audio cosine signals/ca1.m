%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTER ASSIGNMENT 1
%%%%%%%%%%%%%%%%%%%%%%%
% September 13, 2019
% ./ca1.m
% Rebuilds the signals in 4 files which represent compressed audio
% * by: Leomar Duran <https://github.com/lduran2>
% * by: Yacouba Bamba
% * on: 2019-09-03T2000ZQ
% * to: 2019-09-13T1647ZQ
% * submitted for: Temple University, ECE 3512 Signals, Fall 2019 Semester

% clear workspace and load data from file
clear;
% clear commands
clc;
% start the program
main();

%% Tests the method_play function.
function main()
    full_time = 10; % the full time [s] of every audio

    % audio files to play
    audios = { ...
        'audio_1_method_1_no_overlap.mat', ...
        'audio_1_method_2_no_overlap.mat', ...
        'audio_2_method_1_no_overlap.mat', ...
        'audio_2_method_2_no_overlap.mat', ...
    };
    n_audios = size(audios,2); % number of audio files

    % play all of the audio files
    for i_audio = 1:n_audios
        method_play(audios{i_audio}, full_time);
        pause(1); % pause for 1 second
    end % for i_audio = 1:n_audios
end % function main()

%
%% Method_play rebuilds the audio from a compression method. The two
% supported methods divide the audio into a number of windows given by the
% number of rows in the A array (complex exponential coefficients) loaded
% from the file given.
%
% Method_1 stores in each window the dominant
% frequencies per window and their corresponding complex exponential
% coefficients.
%
% Method_2 first finds evenly spaced global frequencies within a range,
% then after dividing the audio into windows, stores in each window the
% complex exponential coefficients corresponding to those global
% frequencies.
%
% Which method is used depends on the number of groups of frequencies.
% Method_1 uses as many frequencie groups as there are windows. If there
% are less frequency groups than there are windows, the default is to use
% method_2.
function method_play(filename, full_time)
    % load data from the file
    audio = load(filename);
    % This loads the following variables:
    %     * fs [Hz]: the sampling rate (samples per second)
    fs = audio.fs;
    %     * freqs [Hz]: a (200 window)x(15 frequencies/window) matrix of
    %                   frequencies
    freqs = audio.freqs;
    %     * A: the complex coefficient of each exponential corresponding to the
    %          frequencies in freqs
    A = audio.A;

    % Get the number of windows and number of frequencies from the matrix of
    % complex exponential coefficients.
    % NOTE: The matrix of coefficients is used because the frequency may
    % just be a vector.
    [n_wins, n_freqs] = size(A);
    % The number of frequencies groups decides whether to use the
    % method 1 or method 2.
    n_freq_groups = size(freqs, 1);
    % if there are at least as many frequencies groups as windows,
    % use method 1
    if (n_freq_groups >= n_wins)
        freq_func = @window_freq;
    else % otherwise use method 2 with 1 group of global frequencies
        freq_func = @global_freq;
    end %end if (n_freq_groups >= n_wins)

    win_time = (full_time / n_wins); % amount of time [s] for each window
    win_fs = (win_time * fs); % the sampling rate [Hz] per window

    full_samples = (fs*full_time); % total number of samples in the full audio
    % the combined audio for 10 seconds, preallocate to 10[s] * 44100[Hz]
    full_audio = zeros(1, full_samples);

    % create a linear vector representing time [s] for input
    t = linspace(0, win_time, win_fs);

    % used as a moving index for full_audio
    % initialize at the beginning of the first window
    i_win_start = 1; % start of the current window
    i_win_end = 0; % end of the current window,
    % warning: orange line because this gets changed before use

    % loop through the windows
    for i_win = 1:n_wins
        % reset x
        x = 0;
        % loop through each frequency in this window
        for i_freq = 1:n_freqs
            % the current frequency in rad/s
            omega = 2*pi*freq_func(freqs, i_win, i_freq);
            % the current complex exponential coefficient
            coefficient = A(i_win,i_freq);
            % get the amplitude from the coefficient
            amp = 2*abs(coefficient);
            % get the phase angle from the coefficient
            phi = angle(coefficient);
            % build the frequency's cosine
            curr_cos = amp*cos(omega*t+phi);
            % add the current frequency cosine to x
            x = x + curr_cos;
        end % for i_freq

        % concatenate the current window to the full audio
        % update the end of the current window
        i_win_end = (i_win_start + win_fs - 1);
        % copy x into the new range
        full_audio(1,i_win_start:i_win_end) = x;
        % update to the start of the next window
        i_win_start = (i_win_end + 1);
    end % for i_win

    % play the full audio
    soundsc(full_audio, fs);
    % pause while playing, because audios will play at the same time
    % otherwise
    pause(full_time);
end % function method_play(filename, full_time)

%
%% Returns the frequency of the given window, specified by index i_freq.
function freq = window_freq(freqs, i_win, i_freq)
    freq = freqs(i_win, i_freq);
end % function freq = window_freq(freqs, i_win, i_freq)

%
%% Returns the global frequency specified by index i_freq. i_win (~1) is
% ignored.
function freq = global_freq(freqs, ~, i_freq)
    freq = freqs(i_freq);
end % function freq = global_freq(freqs, i_win, i_freq)