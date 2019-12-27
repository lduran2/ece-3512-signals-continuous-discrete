%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTER ASSIGNMENT 5
%%%%%%%%%%%%%%%%%%%%%%%
% December 9, 2019
% ./ca5.m
% Tests filters on one of four images.
% by: Duran, Leomar <https://github.com/lduran2>
% by: Yacouba Bamba
% on: 2019-12-09 T20:52ZR
% submitted for: Temple University, ECE 3512 Signals, Fall 2019 Semester

% clear the namespace
clear;

%% Step 1: Load the test images

% flags
do_gaussian = false;
use_highpass = true;
do_fft = false;
% image to filter
i_img = 2;

% filter parameters
HSIZE = 512;
sigma = 0.5;

% image file names
img_names = { 'im_circles.png', 'im_coins.png', 'im_eagles.png', 'im_coins_noisy.png' };

% load the image
[img, cm] = imread(img_names{i_img});
has_cm = (size(cm, 1) ~= 0);

% display the image
figure(1)
surf(peaks); % test image
imshow(img,cm);

%% Step 2: Lowpass filtering the images
if (do_gaussian)
    h = fspecial('gaussian', HSIZE, sigma);
else
    h = fspecial('average', HSIZE);
end % end if (do_gaussian)

% Step 3: Highpass filtering the images
if (~use_highpass)
    % use a filter h as filter i
    i = h;
else
    i = zeros(HSIZE, HSIZE);
    % find midpoint in i
    i_mid_i = (uint64(fix(HSIZE/2 - 1)) + 1);
    % fill the one
    i(i_mid_i,i_mid_i) = 1;
    % subtract
    i = i - h;
end % if (~do_highpass)

if (~do_fft)
    disp('Filtering . . .');
    % show the filtered image
    Y = imfilter(img, i);
    disp('Filtering done.');
    imshow(Y,cm);
else
    % else show the frequency content
    I = myFFT2(i);
    % use color map
    if (has_cm)
        colormap(cm);
    end % if (has_cm)
end % if (~do_fft)

% if no color map, use the default
if (~has_cm)
    colormap default;
end % if (~has_cm)

disp('Done.');
