% Matlab Program to demonstrate the "high pass Filtering of an image using
% 2D-DFT"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm for filtering in the frequency Domain
% Step1: Given ain input image f(x,y) of size M x N, obtain the pading
% parameters P and Q. Typically, we select P = 2M and Q = 2N
% Step2: Form a padded image fp(x,y) of size P X Q by appending the
% necessary number of zeros to f(x,y).
% Step3: Multiply fp(x,y) by (-1)^(x+y)
% Step4: Compute the DFT, F(u,v) of the image from Step 3
% Step5: Generate a Real, Symmetric Filter Function H(u,v) of size P X Q
% with center at coordinates (P/2,Q/2), 
% Step 6:Form the product G(u,v) = H(u,v)F(u,v) using array multiplication
% Obtain the processed image 
% Step 7: gp(x,y) = {real{inverse DFT[G(u,v)]}(-1)^(x+y)
% Step 8: Obtain the final processed result g(x,y) by extracting the M X N region
% from the top, left quadrant of gp(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorith Implementation
% Step1: Given ain input image f(x,y) of size M x N, obtain the pading
% parameters P and Q. Typically, we select P = 2M and Q = 2N
function IO = FreqFilterDFT(Iin, Lev)
b = im2double(Iin);
[m,n] = size(b);
% creating a null array of size 2m X 2n
c = zeros(2*m,2*n);
[p,q] = size(c);
% Step 2
% appdending the original image with the null array to create a padding
% image hence it is Step 2
for i = 1:p
    for j = 1:q
        if i <= m && j<= n
            c(i,j) = b(i,j);
        else
            c(i,j) = 0;
        end
    end
end

% Step 3
% creating a null array of size p X q 
d = zeros(p,q);
%pre processed image for calculating DFT
% Multiplying the padded image with (-1)^(x+y)
for i = 1:p
    for j = 1:q
        d(i,j) = c(i,j).*(-1).^(i + j);
    end
end

% Step 4 
% Computing the 2D DFT using "fft2" matlab command
e = fft2(d);

%%%%%%%%%%%%%%%%%%%%
% Step 5
% Generating the Real, Symmetric Filter Function
% Here we will implement a "Low Pass Filter" using "freqspace" matlab
% command
[x,y] = freqspace(p,'meshgrid');
z = zeros(p,q);
for i = 1:p
    for j = 1:q
        z(i,j) = sqrt(x(i,j).^2 + y(i,j).^2);
    end
end

% Choosing the Cut off Frequency and hence defining the low pass filter
% mask 
H = zeros(p,q);
for i = 1:p
    for j = 1:q
        if z(i,j) <= Lev  % here 0.4 is the cut-off frequency of the LPF
            H(i,j) = 1;
        else
            H(i,j) = 0;
        end
    end
end

% Step 6:Form the product G(u,v) = H(u,v)F(u,v) using array multiplication
% Obtain the processed image 
% from the previous program lines we know that, 
% e : the 2D DFT output of pre processed image
% H : the mask for Low Pass Filter
% let out is the variable 
h1 = e.*H;

% Step 7: gp(x,y) = {real{inverse DFT[G(u,v)]}(-1)^(x+y)
% calculation of inverse 2D DFT of the "out"
h2 = ifft2(h1);

% post process operation 
h3 = zeros(p,q);
for i = 1:p
    for j = 1:q
        h3(i,j) = h2(i,j).*((-1).^(i+j));
    end
end
% figure;
% imshow(h3);title('Post Processed image');
% Step 8: Obtain the final processed result g(x,y) by extracting the M X N region
% from the top, left quadrant of gp(x,y)
% let the smoothed image or low pass filtered image is "out"
out = zeros(m,n);
for i = 1:m
    for j = 1:n
        out(i,j) = h3(i,j);
    end
end
IO = b-out; 
% figure;
% imshow([b IO]);title('input image                 output image');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program developed by : Jagadeesh Samudrala, Asst. Prof, Dept. of ECE,
% Aditya Engineering College, Surampalem, East Godavari Dist, Andhra
% Pradesh, India.
% for any Questions mail me at: samudrala.naren@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%