%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTER ASSIGNMENT 2
%%%%%%%%%%%%%%%%%%%%%%%
% October 3, 2019
% ./ca2.m
% Decodes a message hidden in an audio signal.
% by: Duran, Leomar <https://github.com/lduran2>
% by: Yacouba Bamba
% on: 2019-10-03 T14:00ZQ
% submitted for: Temple University, ECE 3512 Signals, Fall 2019 Semester

clear; % clear namespace
clc; % clear command window

% run the main program
main();


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main program
function main()
    % instead of looping through the messages,
    % a menu is provided since these methods are slow and it is not obvious
    % that the process is running

    % let the user choose a message
    disp('1. Message 1')
    disp('2. Message 2')
    disp('3. Message 3')
    disp('4. Message 4')
    choice = input('Please choose [1/2/3/4/EXIT]: ','s');

    % if the user didn't choose one of the messages, end the main program
    switch (choice)
        case {'1', '2', '3', '4'}
            ;
        otherwise
            return; % exit
    end % switch (choice)

    % whether to show the plots or not
    % the option to show plots is given because the calculation is slow, and
    % the figure(1) window will keep receiving focus
    show_plots = input('Show plots [yes/NO]? ','s');
    should_show_plots = strcmp('yes',show_plots);

    % the eight frequencies encoding the message
    domain = [ 440, 480, 520, 570, 620, 680, 740, 800 ]; %[Hz]
    domain_len = size(domain,2);

    % the tolerance for finding frequencies, each frequency is at least 20 Hz
    % away from the next and previous
    frequency_tolerance = 20;

    % load the message
    message = load(strcat('message_', choice, '.mat'));
    % message_1.mat
    %   struct with fields:
    % 
    %     char_dur_in_sec: 0.1000
    %                  fs: 44100
    %                tone: [1×44100 double]

    string = decodeMessage(message, domain, domain_len, frequency_tolerance, should_show_plots)
end % function menu()


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decodes the specified message, using the given frequency domain and
% tolerance, and allowing the user to pick whether to show plots.
function string = decodeMessage(message, domain, domain_len, frequency_tolerance, should_show_plots)
    % localize the data from the message
    char_dur = message.char_dur_in_sec;
    fs = message.fs;
    tones = message.tone;
    % number of tones in the message
    nTones = size(tones, 2);

    % the length of the tones is the duration times the sampling rate
    % fix to cast as integer
    tones_win_len = fix(char_dur*fs);

    % mirror the frequency domain with the negative of each frequency
    signed_domain = [flip(-domain) domain];
    signed_domain_len = size(signed_domain,2);

    % rather than using strcat, I chose to build a matrix of bytes, and
    % convert to a string afterwards, because strcat was having issues with
    % whitespace characters
    bytes = zeros(0,ceil(nTones/tones_win_len));
    i_byte = 1; % the current byte in the string

    for k = 1:tones_win_len:nTones
        % the tone window taken as a slice of the tones
        tones_win = tones(k:(k + tones_win_len - 1));

        % use myFFT to create Fourier transforms
        [X f] = myFFT(tones_win,fs); % = myFFT(x,fs,varargin)
        X_len = size(X,2); % the number of Fourier transform points

        % split the Fourier transform by magnitude and phase angle
        mag = abs(X); % magnitude
        phi = angle(X); % phase angle [rad]

        % if the plot should be shown
        if (should_show_plots)
            figure(1);
            % plot the Fourier transform
            % by magnitude
            subplot(1,2,1);
            plot(f, mag);
            title('magnitude');

            % by phase angle
            subplot(1,2,2);
            plot(f, phi);
            title('phase angle');
            xlabel('frequency (Hz)');
            ylabel('rad/s');
        end % if (should_show_plots)

        % find the closest frequencies within frequency_tolerance of each
        % value in the signed domain
        found_frequencies = findFrequencies(f, X_len, signed_domain, signed_domain_len, frequency_tolerance);

        % combine the positive and negative magnitudes
        mag_per_f = combinePositiveNegativeMags(mag, found_frequencies, domain);

        % decode the current byte and add it to the byte matrix
        bytes(i_byte) = decode_byte(mag_per_f, domain, domain_len);
        i_byte = (i_byte + 1);
    end % for k = 1:tones_win_len:nTones
    string = char(bytes); % convert the bytes array to a string and return it
end % function decodeMessage(message, domain, domain_len, frequency_tolerance, should_show_plots)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the closest frequencies in f within frequency_tolerance of each
% value in the domain.
function found_frequencies = findFrequencies(f, flen, domain, domain_len, frequency_tolerance)
    % range delimiters,
    first_in_range = 0; % the first element close to the current target frequency
    result = containers.Map; % map of frequencies to their index in f()
    
    % I could have likely made this a lot faster by assuming that all the
    % time windows contain the same range of requencies, about -2 kHz to 2
    % kHz, and using that to calculate the indices of the domain
    % frequencies, but I ran a quick benchmark and the most time appears
    % to be taken in the FFT function.

    % find the closest frequency within frequency_tolerance of each value
    % in the signed domain
    i_domain = 1; % the index of the current frequency being found
    curr_domain = domain(i_domain); % the current frequency being found
    % for all the frequencies in the FFT data
    for i_f = 1:flen
        found_first = (first_in_range > 0); % check if it's already been found
        % check if it's close
        is_close = isCloseToAbs(curr_domain, f(i_f), frequency_tolerance);
        if (~found_first)
            % if not found yet, but is close enough
            if (is_close)
                % save the current index
                first_in_range = i_f;
            end % end if (is_close)
        else % if (found_first)
            % if found, and not close any more
            if (~is_close)
                last_in_range = (i_f - 1); % save the previous index
                % find the middle index, fixed as an integer
                median_in_range = fix(median([first_in_range, last_in_range]));
                % save the middle index by mapping the frequency to it
                % Maps require string keys, so num2str is used.
                result(num2str(domain(i_domain))) = median_in_range;
                % reset the first index
                first_in_range = 0;
                % move to the next domain frequency by index
                i_domain = (i_domain + 1);
                % if the domain frequencies have all been found, exit the
                % loop
                if (i_domain > domain_len)
                    break;
                end % if (i_domain > domain_len)
                % move to the next domain frequency by reference
                curr_domain = domain(i_domain);
            end % if (~is_close)
        end % if (found_first)
    end % for i_f = 1:flen
	found_frequencies = result;
end % function findFrequencies(f, flen, domain, domain_len, frequency_tolerance)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combines the magnitudes matching the positive and negative frequencies.
function mag_per_f = combinePositiveNegativeMags(mag, frequencies, domain)
    result = containers.Map; % map of absolute frequencies to combined magnitudes
	% for each frequency n the domain
    for f=domain
        % add the positive and negative magnitudes k/2, to make k
        % maps require string keys, so num2str is used
        result(num2str(f)) = ...
            mag(frequencies(num2str(f))) ...
            + mag(frequencies(num2str(-f)));
    end % end for l=domain
    mag_per_f = result;
end % function combinePositiveNegativeMags(frequencies, mag)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decodes the byte represented by the specified magnitudes from the given
% domain of the given length.
function byte = decode_byte(mag_per_f, domain, len)
    % find the tolerance for magnitude by finding the local maximum and
    % dividing it into the 55%ile
    magnitude_tolerance = (max(cell2mat(values(mag_per_f)))*.55);

	% create a byte by cycling through the frequencies
    result = 0;
    for i_bit = 1:len
		% if there is enough energy, add the current bit to the result
        if (mag_per_f(num2str(domain(i_bit))) >= magnitude_tolerance)
            result = result + 2^(i_bit - 1);
        end % if (mag_per_f(num2str(domain(i_bit))) > magnitude_tolerance)
    end % for i_bit = 1:len
    byte = result;
end % function decode_byte(mag_per_f)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns if the candidate is close enough to the target, given the
% specified tolerance.
function bool = isCloseToAbs(target, candidate, tolerance)
    bool = (abs(candidate - target) < tolerance);
end % end function isCloseToAbs(target, candidate, tolerance)
