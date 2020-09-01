function [phase2, phase1, dim_out] = HTWrapping(I0,NS,R,HLPF,RemoveCarrier, dim_in)
% wrapping an open fringe pattern with Hilbert transform 


% Inputs: I0    -> input fringe pattern
%         NS    -> zero padding size
%         R     -> default radius for the band pass filter
%         HLPF    -> low pass filter
%         RemoveCarrier -> allows for removing carrier with Fourier method
%         or not (1: enable, otherwise disable)
%         dim_in -> an option to take the exact position of the band pass
%         filter as input
 
% Outputs: phase2  -> No Carrier Wrapped Phase
%          phase1  -> Wrapped Phase
%          dim_out -> the position of the selected filter in Fourier domain 


% created by: Parsa Omidi
% last modified: 2019
%%
y = hilbert(I0');
Ih = y';
I1 = imag(Ih);
% I0 = (I0-min(I0(:)))./(max(I0(:))-min(I0(:)));
% I1 = (I1-min(I1(:)))./(max(I1(:))-min(I1(:)));
phase1(:,:) = atan2(I1 , I0);
if RemoveCarrier==1
    [I, Drows, Dcols] = ZeroPadding(phase1, NS);
    I_data = fft2(I(:,:));
    I_f_shift = fftshift(I_data);
    ImageDisplay(abs(I_f_shift));
    switch nargin
        case 6
            vertices = dim_in;
            dim_out = vertices;
            X=vertices(:,1);
            Y=vertices(:,2);
            [m,n] = size(I_data);
            BW1 = poly2mask (X,Y,m,n);
        case 5
            %         [peak_col,peak_row] =  ginput(1);
            %         h = imellipse(gca,[round(peak_col-R/2),round(peak_row-R/2),R,R]);
            %         BW1 = createMask(h);
            h = imellipse(gca,[20,20,R,R]);
            vertices = wait(h);
            X=vertices(:,1);
            Y=vertices(:,2);
            [m,n] = size(I_data);
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
    % Shift the peak to the center
    I_f_shift = circshift(I_f_onlyPeak,[(round(size(BW2,1)/2 - peak_row)) (round(size(BW2,2)/2 - peak_col))]);
    % Smooth the signal edges in frequency domain befor
    % applying inverse FFT
    I_f_shift = I_f_shift .* HLPF;
    I_f_shift = ifftshift(I_f_shift);
    I = ifft2(I_f_shift);
    phase2 = angle(I);
    phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
else
    phase2  = [];
    switch nargin
        case 6
            dim_out = dim_in;
        case 5
            dim_out= [];
    end
end
end