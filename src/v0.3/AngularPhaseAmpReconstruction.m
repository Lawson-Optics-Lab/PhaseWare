function [Hr] = AngularPhaseAmpReconstruction(hologram,reconstruction_distance,...
wavelength,pixel_size,phase_mask,scale_factor)

%This function reconstructs a hologram using angular spectrum
%and returns the reconstructed hologram matrix
%This function ONLY works for Gabor (in-line) Holograms
%Off-axis reconstructions will only return the zero-order

%The angular spectrum method, as implemented here, reconstructs image
%pixels that are equal in size to the CCD pixels. 

%reconstruction_distance in METERS
%wavelength in METERS
%pixel size in METERS

%phase_mask defaults to uniform phase of constant intensity (plane wave)
%(i.e. this can be used to reconstruct with some tilt or lens phase, etc.)

%scale_factor will change the size of the original hologram
%this has little effect on the angular spectrum method, but can be used to
%increase the apparent numerical resolution
%scale_factor defaults to 1

%Padding the hologram has no significant effect on this method, and is
%therefore not included in this function

% Reference: Nehmetallah, Georges T., Rola Aylo, and Logan Williams...
%     "Analog and digital holography with matlab"
%     Society of Photo-Optical Instrumentation Engineers (SPIE), 2015.
% Author: Georges T. Nehmetallah
%%
A = size(hologram);

if nargin < 6
    scale_factor = 1;   
end
if nargin < 5
    phase_mask = 1; 
end
if nargin < 4
    error('Not enough input arguments')    
end
if length(phase_mask) == 1
     phase_mask = ones(A).*phase_mask; 
% phase_mask = 1;
end
B = phase_mask;
H = hologram;
d = reconstruction_distance;
w = wavelength;
dx = pixel_size;

%if the image is not square, crop to be square along the smallest dimension
if A(1,1)<A(1,2)
    %crop the x axis of the image to match the y axis
    %crops symmetrically about the center
    H = H(:,round(A(1,2)/2+1-A(1,1)/2):round(A(1,2)/2+A(1,1)/2));
    B = B(:,round(A(1,2)/2+1-A(1,1)/2):round(A(1,2)/2+A(1,1)/2));   
end

if A(1,1)>A(1,2)  
    %crop the y axis of the image to match the x axis
    H = H'; %transpose
    B = B';    
    C = size(H);  %this is the line that fixes this part
    %after transposing, the sizes didn't match anymore    
    H = H(:,round(C(1,2)/2+1-C(1,1)/2):round(C(1,2)/2+C(1,1)/2)); %crop
    B = B(:,round(C(1,2)/2+1-C(1,1)/2):round(C(1,2)/2+C(1,1)/2)); %crop
    H = H'; %transpose back
    B = B';
    %not sure if this works perfectly yet - havent had to try it
end

%use double precision to allow for complex numbers
H = double(H);
B = double(B);

%if the image is too big to compute efficiently, scale it
%scale is from 0 to 1, in percentage 
%to scale above 1, the phase mask can cause problems because the phase
%does not scale well with bicubic resampling
if scale_factor ~= 1
    H = imresize(H, scale_factor,'bicubic');  
    if length(B(:,1)) > 1  %only rescale profile if it's not a scalar
        B = imresize(B, scale_factor,'bicubic'); 
    end
    dx=dx/scale_factor;
end

Image_Pixel_Size = dx;
n = length(H);   %size of hologram matrix nxn
H = double(H)-mean(mean(H));  %must be class double for imaginary #s

%%%%%%Angular Spectrum Method
k0=2*pi/w;              %k-vector
k = (-n/2:1:n/2-1)*dx;  %array same dimentions as hologram
l = (-n/2:1:n/2-1)*dx;
[XX,YY] = meshgrid(k,l);
step = k(2)-k(1);  %step size

k_max = 2*pi/step;  %maximum k vector
k_axis = linspace(0,k_max,n)-k_max/2;  %confine reconstruction to k<k_max 
[Kx Ky]=meshgrid(k_axis,k_axis);

E=B.*exp(-1i*sqrt(k0^2-Kx.^2-Ky.^2)*d-1i*k0*d); %Angular spectrum
Hf=fftshift(fft2(sqrt(H)));
[Hr]=ifft2((Hf.*E));    
Mag_Rec_image=abs(Hr);
% figure
% imagesc(Mag_Rec_image);
% colormap gray;
end %end function
    
 


