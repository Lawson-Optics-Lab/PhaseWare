function [I_out, BG] = BackgroundRemoval (I_in, BG0)
% Function to remove background from a grayscale image. 


% Inputs: I_in -> input image to remove background
%         BG0  -> input Background
% when BG0 is not used, the function estimates the background automatically
% in this case, it estimates small variations in intensity as a 2D polynomial  
% Outputs: I_out -> background removed image
%          BG    -> entered or estimated background 


% created by: Parsa Omidi
% last modified on Feb 2017
%%
switch nargin
    case 2
        I = medfilt2(I_in);
        h = fspecial('average',4);
        I = imfilter(I,h);
        
        %                 BG = GetModelBackgroundImage(I); % Estimated Background
        BG = BG0;
        I_out = I-BG;
        h = fspecial ('gaussian' , 40, 4);
        I_out = medfilt2(I_out);
        I_out = imfilter(I_out,h);
        I_out = I_out(10:end-10, 10:end-10);
    otherwise
        I = medfilt2(I_in);
        h = fspecial('average',4);
        I = imfilter(I,h);
        
        BG = GetModelBackgroundImage(I); % Estimated Background
        I_out = I-BG;
        h = fspecial ('gaussian' , 40, 4);
        I_out = medfilt2(I_out);
        I_out = imfilter(I_out,h);
        I_out = I_out(10:end-10, 10:end-10);
end

end