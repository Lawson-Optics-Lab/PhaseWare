function rec_image = ConvolutionPhaseAmpReconstruction(Interferogram, FocusDistance, lambda, PixelSize, Magnification)
% Reference: Nehmetallah, Georges T., Rola Aylo, and Logan Williams...
%     "Analog and digital holography with matlab"
%     Society of Photo-Optical Instrumentation Engineers (SPIE), 2015.
% Author: Georges T. Nehmetallah
% Modified by Parsa Omidi, 2018
% lambda: optical wavelength [m]
[Ny,Nx] = size(Interferogram);
minN=min(Ny,Nx);
I=Interferogram(1:minN,1:minN);% Crop to have rows equal to columns
[Ny,Nx]=size(I);
nx = [-Nx/2:1:Nx/2-1];
ny = [-Ny/2:1:Ny/2-1]';
% lambda = 1000*lambda;
X = nx*PixelSize;
Y = ny*PixelSize;
[XX, YY]=meshgrid(X,Y);
X2 = nx/(Nx*PixelSize);
Y2 = ny/(Ny*PixelSize);
[XX2, YY2]=meshgrid(X2,Y2);
% num = exp(1i*2*pi/lambda*(d^2+XX.^2+YY.^2).^0.5);
% den = (d^2+XX.^2+YY.^2).^0.5;
% g = -1i/lambda*num./den;
% rec_image = fftshift(fft2(fft2(I).*fft2(g)));

%%

dp=FocusDistance*Magnification;
f=(1/FocusDistance+1/dp)^-1; % See Eq. 3.36 in ref [5]
L=exp(1i*pi/(lambda*f)*(XX.^2+YY.^2)); % See Eq. 3.37 in ref [5]
P=exp(1i*pi*lambda*FocusDistance^2/f*(XX2.^2+YY2.^2)); % See Eq. 3.7 in ref [5]
num=exp(-1i*2*pi/lambda*(dp^2+XX.^2+YY.^2).^0.5);
den=(dp^2+XX.^2+YY.^2).^0.5;
g=1i/lambda*num./den;% See Eq. 3.28 in ref [5]
rec_image=fftshift(P.*ifft2(fft2(I.*L).*fft2(g)));
Mag_Rec_image=abs(rec_image);
% figure
% imagesc(Mag_Rec_image);
% colormap gray;
end