function [I_Out,rect_crop_Out]  = FramePreProcessing(I_In, NS, rect_crop_In)
% Pre-processing 


% Inputs: I_In    -> input fringe pattern
%         NS     -> size of the zero padded image
%         rect_crop_In -> an option to take positions for cropping the
%         frame after zero padding
 
% Outputs: I_Out -> output image
%          rect_crop_Out -> the position used for cropping the input image


% created by: Parsa Omidi
% last modified: 2018
switch nargin
    case 2
        I = I_In;
        % Zero Padding
        [Rows,Cols]=size(I);
        Drows = round((NS - Rows)/2);
        Dcols = round((NS - Cols)/2);
        I = padarray(I,[Drows Dcols],0,'both');
        I = I +1i*zeros(size(I));
        % Rectangle Cropping
        [r,c]=size(I);
        if c > r
            rect_crop=[floor((c-r)/2) 1 r r];
        elseif c < r
            rect_crop=[1 floor((r-c)/2) c c];
        else
            rect_crop=[1 1 c r];
        end
        I_Out=imcrop(I,rect_crop);
        rect_crop_Out = rect_crop;
    case 3
        I = I_In;
        % Zero Padding
        [Rows,Cols]=size(I);
        Drows = round((NS - Rows)/2);
        Dcols = round((NS - Cols)/2);
        I = padarray(I,[Drows Dcols],0,'both');
        I = I +1i*zeros(size(I));
        I_Out = imcrop(I, rect_crop_In);
        rect_crop_Out = rect_crop_In;
end
end