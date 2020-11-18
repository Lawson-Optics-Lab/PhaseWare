
function modeledColorChannel = ModelOneColorChannel(colorChannel)
	try
		modeledColorChannel = colorChannel; % Initialize.
		rows = size(colorChannel, 1);
		columns = size(colorChannel, 2);
% 		midX = columns / 2;
% 		midY = rows / 2;
		[X,Y] = meshgrid(1:columns, 1:rows);
		z = colorChannel;
		x1d = reshape(X, numel(X), 1);
		y1d = reshape(Y, numel(Y), 1);
		z1d = double(reshape(z, numel(z), 1));
		x = [x1d y1d];
		polynomialOrder = 3;
		p = polyfitn(x, z1d, polynomialOrder);
        % Evaluate on a grid and plot:
		zg = polyvaln(p, x);
		modeledColorChannel = reshape(zg, [rows columns]);
		
	catch ME
		errorMessage = sprintf('Error in ModelOneColorChannel():\n%s', ME.message);
		if strfind(errorMessage, 'HELP MEMORY')
			memory;
			errorMessage = sprintf('%s\n\nCheck console for memory information.', errorMessage);
		end
		uiwait(warndlg(errorMessage));
	end
	% Free up memory.  Should be automatically freed but we have had problems here.
	% (Make sure you don't clear any input or output arguments or you'll get an exception.)
	clear('zg', 'x1d', 'y1d', 'z1d', 'X', 'Y', 'x', 'z', 'p');
	return; % from ModelOneColorChannel

