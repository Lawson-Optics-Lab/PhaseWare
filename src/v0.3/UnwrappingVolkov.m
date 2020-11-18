function imOut = UnwrappingVolkov(imIn)
imIn = [imIn fliplr(imIn); flipud(imIn) flipud(fliplr(imIn))];
imIn = mod(imIn,2*pi);
[Ny, Nx] = size(imIn);
Z = exp(1i*imIn);

qx = ([0:Nx/2 (-Nx/2+1):-1] + 0.01)*2*pi/Nx;
qy = ([0:Ny/2 (-Ny/2+1):-1] + 0.01)*2*pi/Ny;
[qxg, qyg] = meshgrid(qx,qy);
[dIdx, dIdy] = gradient(imIn);
[dZdx, dZdy] = gradient(Z);
dkdy = (real(dZdy./(1i*Z)) - dIdy)/(2*pi);
dkdx = (real(dZdx./(1i*Z)) - dIdx)/(2*pi);

k = real( (1/(1i))*(ifft2((fft2(dkdx).*qxg + fft2(dkdy).*qyg)./(qxg.^2 + qyg.^2)) ));
k = k - min(k(:));
k = round(k);
k = k(1:end/2, 1:end/2);

imIn = imIn(1:end/2, 1:end/2);
imOut = imIn + 2*pi*k;
imOut = imOut - min(imOut(:));
end