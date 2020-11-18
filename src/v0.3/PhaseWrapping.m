function [phase2, phase1, dim_out] = PhaseWrapping(I_in, wrappingAlgorithm, NS2,RemoveCarrier,FrequencyFilterType,FrequencyFilterSize,WFT_parameters,filterRadious,dim_in)
%% Wrapping Algorithms for open fringe patterns
%Inputs:
% I_in: a single frame fringe pattern
% wrappingAlgorithm: options for algorithms: 1->complex amplitude data,
% 2->Fourier (FFT), 3-> Wavelet, 4-> Hilbert, 5-> Windowed FFT.
% NS2: Zero padding size
% RemoveCarrier: option to remove carrier (0-> donnot remove, 1-> remove carrier)
% FrequencyFilterType: Type of filter to select the main component in frequency domain
% FrequencyFilterSize: Indicate the cutoff frequency and should be a positive number
% filterRadious: The filter is a circular window and this parameter shows the number of pixels for the radius of it
% dim_in: optional input.  Indicates the coordinates of the filter. When its empty, a matlab figure will apear to allow user select a single coordinate which approximately shows the center of the filter and the press enter


%Outputs:
% phase1:  wrapped phase map with carrier
% phase2:  Carrier removed wrapped phase map with Fourier method. When RemoveCarrier=0, phase2 will bring an empty value
% dim_out: coordinates of the selected filter in frequency domain. When the user input a specific coordinate for the dim_in, the dim_out will be equal to dim_in
[ny,nx]=size(I_in);
% filterRadious = round(ny/4);
HLPF = fftshift(lpfilter(FrequencyFilterType, ny, nx,FrequencyFilterSize));%'gaussian', 'ideal', 'btw'
[HLPF, ~, ~] = ZeroPadding(HLPF, NS2);
switch wrappingAlgorithm
    case 0
        [file,path] = uigetfile('*.m');
        if (ischar(path))
            selectedfile = fullfile(path,file);
            addpath(path);
            [~,name,~] = fileparts(selectedfile);
            fh = str2func(name);
            for ii = 1:NoF
                I0 = V0(:,:,ii);
                [phase2,phase1] = double(fh(I0));
            end
        end
    case 1
        [ny,nx] = size(I_in);
        ImageVector = I_in(:);
        ind = find(real(ImageVector)>0);
        I_phase = zeros(1,length(ImageVector))';
        I_phase(ind) = atan(imag(ImageVector(ind))./real(ImageVector(ind)));
        ind = find(real(ImageVector)<=0);
        I_phase(ind) = atan(imag(ImageVector(ind))./real(ImageVector(ind)))+pi*sign(imag(ImageVector(ind)));
        phase1 = reshape(I_phase,ny,nx);
        if RemoveCarrier==1
            [I, Drows, Dcols] = ZeroPadding(real(phase1), NS2);
            I_data = fft2(I);
            I_f_shift = fftshift(I_data);
            ImageDisplay(abs(I_f_shift));
            switch nargin
                case 9
                    %                 [I, dim_out] = CarrierRemoval (I,  filterRadious,HLPF,dim_in);
                    vertices = dim_in;
                    dim_out = vertices;
                    X=vertices(:,1);
                    Y=vertices(:,2);
                    [m,n] = size(I_data);
                    BW1 = poly2mask (X,Y,m,n);
                case 8
                    %                 [I, dim_out] = CarrierRemoval (I,  filterRadious,HLPF);
                    h = imellipse(gca,[20,20,filterRadious,filterRadious]);
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
            I_f_shift = circshift(I_f_onlyPeak,[(round(size(I_f_shift,1)/2 - peak_row)) (round(size(I_f_shift,2)/2 - peak_col))]);
            % Smooth the signal edges in frequency domain befor
            % applying inverse FFT
            I_f_shift = I_f_shift .* HLPF;
            I_f_shift = ifftshift(I_f_shift);
            I = ifft2(I_f_shift);
            phase2 = angle(I);
            phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
        else
            phase2 = [];
            switch nargin
                case 9
                    dim_out = dim_in;
                case 8
                    dim_out= [];
            end
        end
    case 2
        switch nargin
            case 9
                [I, Drows, Dcols] = ZeroPadding(real(I_in), NS2);
                [phase2, phase1, dim_out] = FFTWrapping(I,filterRadious,HLPF,RemoveCarrier,dim_in);
                phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
                phase1 = phase1(Drows:end-Drows,Dcols:end-Dcols);
            case 8
                [I, Drows, Dcols] = ZeroPadding(real(I_in), NS2);
                [phase2, phase1, dim_out] = FFTWrapping(I,filterRadious,HLPF,RemoveCarrier);
                phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
                phase1 = phase1(Drows:end-Drows,Dcols:end-Dcols);
        end
    case 3
        switch nargin
            case 9
%                 [I, Drows, Dcols] = ZeroPadding(real(I_in), NS2);
                [phase2, phase1, dim_out] = WTWrapping(real(I_in),NS2, filterRadious,HLPF,RemoveCarrier, dim_in);
%                 phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
%                 phase1 = phase1(Drows:end-Drows,Dcols:end-Dcols);
            case 8
%                 [I, Drows, Dcols] = ZeroPadding(real(I_in), NS2);
                [phase2, phase1, dim_out] = WTWrapping(real(I_in),NS2, filterRadious,HLPF,RemoveCarrier);
%                 phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
%                 phase1 = phase1(Drows:end-Drows,Dcols:end-Dcols);
        end
    case 4
        switch nargin
            case 9
%                 [I, Drows, Dcols] = ZeroPadding(real(I_in), NS2);
                [phase2, phase1, dim_out] = HTWrapping(I_in,NS2, filterRadious,HLPF,RemoveCarrier, dim_in);
%                 phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
%                 phase1 = phase1(Drows:end-Drows,Dcols:end-Dcols);
            case 8
%                 [I, Drows, Dcols] = ZeroPadding(real(I_in), NS2);
                [phase2, phase1, dim_out] = HTWrapping(I_in,NS2, filterRadious,HLPF,RemoveCarrier);
%                 phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
%                 phase1 = phase1(Drows:end-Drows,Dcols:end-Dcols);
        end
    case 5
        switch nargin
            case 9
%                 [I, Drows, Dcols] = ZeroPadding(real(I_in), NS2);
                [phase2, phase1, dim_out] = WFTWrapping(real(I_in),NS2, filterRadious,HLPF,RemoveCarrier, WFT_parameters,dim_in);
%                 phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
%                 phase1 = phase1(Drows:end-Drows,Dcols:end-Dcols);
            case 8
%                 [I, Drows, Dcols] = ZeroPadding(real(I_in), NS2);
                [phase2, phase1, dim_out] = WFTWrapping(real(I_in),NS2, filterRadious,HLPF,RemoveCarrier, WFT_parameters);
%                 phase2 = phase2(Drows:end-Drows,Dcols:end-Dcols);
%                 phase1 = phase1(Drows:end-Drows,Dcols:end-Dcols);
        end
end
end
