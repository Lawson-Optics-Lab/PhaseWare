
% Technique adapted from:
% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
% 
% Inputs: 1. Complex image in .mat double format
%         2. Binary mask (optional)          
% Outputs: 1. Unwrapped phase image
%          2. Phase quality map

function [imOut, maskOut] = UnwrappingGhigliaPritt(imIn, mask) 
im_mask=ones(size(imIn));                     %Mask (if applicable)

im_phase=imIn;                         %Phase image
imOut=zeros(size(imIn));               %Zero starting matrix for unwrapped phase
adjoin=zeros(size(imIn));                     %Zero starting matrix for adjoin matrix
unwrapped_binary=zeros(size(imIn));           %Binary image to mark unwrapped pixels

%% Calculate phase quality map
im_phase_quality=PhaseDerivativeVariance(im_phase);   

%% Identify starting seed point on a phase quality map
minp=im_phase_quality(2:end-1, 2:end-1); minp=min(minp(:));
maxp=im_phase_quality(2:end-1, 2:end-1); maxp=max(maxp(:));
figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis square, axis off; title('Phase quality map'); 
uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
[xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm

%% Unwrap
colref=round(xpoint); rowref=round(ypoint);
imOut(rowref,colref)=im_phase(rowref,colref);                        %Save the unwrapped values
unwrapped_binary(rowref,colref,1)=1;
if im_mask(rowref-1, colref, 1)==1 adjoin(rowref-1, colref, 1)=1; end       %Mark the pixels adjoining the selected point
if im_mask(rowref+1, colref, 1)==1 adjoin(rowref+1, colref, 1)=1; end
if im_mask(rowref, colref-1, 1)==1 adjoin(rowref, colref-1, 1)=1; end
if im_mask(rowref, colref+1, 1)==1 adjoin(rowref, colref+1, 1)=1; end
imOut=GuidedFloodFill(im_phase, imOut, unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap

figure; imagesc(im_phase), colormap(gray), axis square, axis off; title('Wrapped phase'); 
figure; imagesc(imOut), colormap(gray), axis square, axis off; title('Unwrapped phase'); 
end
