function [NoCarrierWrappedPhase,WrappedPhase, dim_out] = FFTWrapping(I0,R,HLPF,RemoveCarrier,dim_in)

% wrapping an open fringe pattern with Fourier transform method


% Inputs: I0    -> input fringe pattern
%         R     -> default radius for the band pass filter
%         HLPF    -> low pass filter
%         RemoveCarrier -> allows for removing carrier with Fourier method
%         or not (1: enable, otherwise disable)
%         dim_in -> an option to take the exact position of the band pass
%         filter as input
 
% Outputs: NoCarrierWrappedPhase -> No Carrier Wrapped Phase
%          WrappedPhase -> Wrapped Phase
%          dim_out -> the position of the selected filter in Fourier domain 


% created by: Parsa Omidi
% last modified: 2018

I_f = fft2(I0);
I_f_shift = fftshift(I_f);
ImageDisplay(abs(I_f_shift), [])
% Select a location indicates one of the peak area (approximately)
% Generate a mask around the selected peak
switch nargin
    case 5
        vertices = dim_in;
        dim_out = vertices;
        X=vertices(:,1);
        Y=vertices(:,2);
        [m,n] = size(I_f);
        BW1 = poly2mask (X,Y,m,n);
    case 4
        %         [peak_col,peak_row] =  ginput(1);
        %         h = imellipse(gca,[round(peak_col-R/2),round(peak_row-R/2),R,R]);
        %         BW1 = createMask(h);
        h = imellipse(gca,[20,20,R,R]);
        vertices = wait(h);
        X=vertices(:,1);
        Y=vertices(:,2);
        [m,n] = size(I_f);
        BW1 = poly2mask (X,Y,m,n);
        dim_out = vertices;
end

% Mask out the peak from the rest of the signal
I_f_onlyPeakTemp = BW1.*I_f_shift;
% Refined the Peak location and the mask
peakValue = max(max(abs(I_f_onlyPeakTemp)));
[peak_row,peak_col,~] = find(abs(I_f_onlyPeakTemp) == peakValue);
Ry2 = round(max(Y)-min(Y));
Rx2 = round(max(X)-min(X));
h = imellipse(gca,[round(peak_col-Rx2/2),round(peak_row-Ry2/2),Rx2,Ry2]);
BW2 = createMask(h);
I_f_onlyPeak = BW2.*I_f_shift;
WrappedPhase = angle(ifft2(ifftshift(I_f_onlyPeak)));
if RemoveCarrier ==1
    % Shift the peak to the center
    I_f_shift = circshift(I_f_onlyPeak,[(round(size(I0,1)/2 - peak_row)) (round(size(I0,2)/2 - peak_col))]);
    % Smooth the signal edges in frequency domain befor
    % applying inverse FFT
    I_f_shift = I_f_shift .* HLPF;
    I_f_shift = ifftshift(I_f_shift);
    I = ifft2(I_f_shift);
    NoCarrierWrappedPhase = angle(I);
else
    NoCarrierWrappedPhase = [];
    switch nargin
        case 5
            dim_out = dim_in;
        case 4
            dim_out= [];
    end
end
end