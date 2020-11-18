function [ gain, bias, out] = ImageDisplay (I, nstd, GrayContrast, colorMap,titleText)
% Display image with auto-scaling for grayscale images

% Convert to double
I = double(I);
% Get size info
[ny, nx] = size(I);
% Set up the default axes
if nargin < 2 || isempty(nstd)
    nstd = 2; % default
end
if nargin < 3 || isempty(GrayContrast)
    GrayContrast = 256; % default
    colorMap = gray(GrayContrast);
end
if nargin < 4 || isempty(colorMap)|| strcmp(colorMap, 'gray')
    colorMap = gray(GrayContrast);
end
if nargin < 5 || isempty(titleText)
    titleText = [];
end
% Open new figure window if indicated


% Treat by mapping nstd standard deviations
% of input range into display range

% number of pixels to use for statistics (mean and std)
max_num = min( [ 10000, ny*nx ] );
I = I(:,:,1);
USE = round( linspace(1,ny*nx,max_num));
inmn = nanmean(I(USE));
instd = nanstd(I(USE));
gain = 128/(instd*nstd);
bias = 128-(inmn*128)/(instd*nstd);
I_out = I*gain+bias;
figure; image(I_out);
set(gca,'dataAspectRatio',[1 1 1])
colormap(colorMap);
title(titleText);
end

