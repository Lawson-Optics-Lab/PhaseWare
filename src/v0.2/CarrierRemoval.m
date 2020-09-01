function [I_out, rect2] = CarrierRemoval (I_in, Radius,HLPF, rect_in)

% Function to remove carrier frequency from a single-carrier fringe pattern image. 
% by using FFT method


% Inputs: I_in    -> input fringe pattern
%         Radius  -> radius of band pass filter
%         HLPF    -> low pass filter
%         rect_in -> an option to take the exact position of the band pass
%         filter as input
 
% Outputs: I_out -> carrier removed pattern
%          rect2 -> the position of the seletced filetr in Fourier domain 


% created by: Parsa Omidi
% last modified: 2017
switch nargin
    case 3
        IF = (fft2(fftshift(I_in)));
        IF_s = fftshift(IF);
        abs_I = abs(IF_s).^2;
        ImageDisplay(abs_I,2,512,[],'Amp FFT of the  Z-Padded Phase');
        % res = bhp(I1F_s,thresh,n);
        % im(abs(res).^2)
        % colormap(gray(256))
        img_size = size (I_in);
        RR = Radius;
        h = imellipse(gca,[20,20,RR,RR]);
        vertices = wait(h);
        X=vertices(:,1);
        Y=vertices(:,2);
        [m,n] = size(IF_s);
        bw2 = poly2mask (X,Y,m,n);
        IF_s_selection=bw2.*IF_s;
        ImageDisplay(abs(IF_s_selection))
        % phase_Rec_image1=bw.*phase_Rec_image1;
        % I1F_s_selection = I1F_s_selection(round(min(Y)):round(max(Y)),round(min(X)):round(max(X)));
        peak_value = max(max(abs(IF_s_selection))); [peak_row,peak_col,~] = find(abs(IF_s_selection) == peak_value); % Finding of peak location.
        % Shifting peak to the center
        F_shift_to_center = circshift(IF_s_selection,[(round(img_size(1)/2) - peak_row) (round(img_size(2)/2) - peak_col)]);
        IF_s_selection = F_shift_to_center.*HLPF;
        % im(abs(I1F_s_selection))
        %                 I1=IF_s_selection;
        % IFFT: transform to space domain
        I_out=fftshift(ifft2(fftshift(IF_s_selection)));
        rect2 = vertices;
    case 4
        IF = (fft2(fftshift(I_in)));
        IF_s = fftshift(IF);
        img_size = size (I_in);
         vertices = rect_in;
        X=vertices(:,1);
        Y=vertices(:,2);
        [m,n] = size(IF_s);
        bw2 = poly2mask (X,Y,m,n);
        IF_s_selection=bw2.*IF_s;
%         im(abs(IF_s_selection))
        % phase_Rec_image1=bw.*phase_Rec_image1;
        % I1F_s_selection = I1F_s_selection(round(min(Y)):round(max(Y)),round(min(X)):round(max(X)));
        peak_value = max(max(abs(IF_s_selection))); [peak_row,peak_col,~] = find(abs(IF_s_selection) == peak_value); % Finding of peak location.
        % Shifting peak to the center
        F_shift_to_center = circshift(IF_s_selection,[(round(img_size(1)/2) - peak_row) (round(img_size(2)/2) - peak_col)]);
        IF_s_selection = F_shift_to_center.*HLPF;
        % im(abs(I1F_s_selection))
        %                 I1=IF_s_selection;
        % IFFT: transform to space domain
        I_out=fftshift(ifft2(fftshift(IF_s_selection)));
        %         I_out = I_out(Drows:end-Drows,Dcols:end-Dcols);
        rect2 = rect_in;
end

