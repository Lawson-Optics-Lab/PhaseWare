% Takes a color, noisy, non-uniform background image and fits a 2D model to it.
% Returns a noise-free, smooth color image where each color plane is independently
% normalized with a max value of 1.  See lengthy comment farther down.
% modeledBackgroundImage is a double image in the range 0-1.
function modeledBackgroundImage = GetModelBackgroundImage(nonuniformBackgroundImage)
	try
	% Get the dimensions of the image.  numberOfColorBands should be = 3.
	% [rows columns numberOfColorBands] = size(nonuniformBackgroundImage);
	% Preallocate output image.
	modeledBackgroundImage = zeros(size(nonuniformBackgroundImage));
	[rows columns numberOfColorChannels] = size(nonuniformBackgroundImage);
	if numberOfColorChannels > 1
		% It's a color image.  Correct each color channel one at a time.
		for colorIndex = 1 : 3
			oneColorChannel = nonuniformBackgroundImage(:, :, colorIndex);
			% Model it to get rid of noise, scratches, smudges, etc.
			noiselessImage = ModelOneColorChannel(oneColorChannel);
			% Divide the actual image by the modeled image.
			maxValue = max(max(noiselessImage));
			% Divide by the max of the modeled image so that you can get a "percentage" image 
			% that you can later divide by to get the color corrected image.
			% Note: each color channel may have a different max value,
			% but you do not want to divide by the max of all three color channel
			% because if you do that you'll also be making all color channels have the same
			% range, in effect "whitening" the image in addition to correcting for non-uniform 
			% background intensity.  This is not what you want to do.
			% You want to correct for background non-uniformity WITHOUT changing the color.  
			% That can be done as a separate step if desired.
			noiselessImage = noiselessImage / maxValue; % This is a percentage image.
			% Now insert that percentage image for this one color channel into the color percentage image.
			modeledBackgroundImage(:, :, colorIndex) = noiselessImage;
		% 	minValue = min(modeledBackgroundImage(:))
		% 	maxValue = max(modeledBackgroundImage(:))
		end
	else
		% It's a gray scale image.  (Much simpler situation.)
		% Model it to get rid of noise, scratches, smudges, etc.
		noiselessImage = ModelOneColorChannel(nonuniformBackgroundImage);
		% Divide the actual image by the modeled image.
		maxValue = max(max(noiselessImage));
		% Divide by the max of the modeled image so that you can get a "percentage" image 
		% that you can later divide by to get the color corrected image.
% 		modeledBackgroundImage = noiselessImage / maxValue; % This is a percentage image.
		% 	minValue = min(modeledBackgroundImage(:))
		% 	maxValue = max(modeledBackgroundImage(:))
        modeledBackgroundImage = noiselessImage; % 
	end
	catch ME
		errorMessage = sprintf('Error in function GetSmoothedBackgroundImage().\n\nError Message:\n%s', ME.message);
		fprintf(1, '%s\n', errorMessage);
		uiwait(warndlg(errorMessage));
	end

end
	
