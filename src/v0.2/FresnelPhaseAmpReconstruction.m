function rec_image = FresnelPhaseAmpReconstruction(Interferogram, FocusDistance, lambda, PixelSize)
% Fresnel reconstruction for digital holographic interferometry

% Reference: Nehmetallah, Georges T., Rola Aylo, and Logan Williams...
%     "Analog and digital holography with matlab"
%     Society of Photo-Optical Instrumentation Engineers (SPIE), 2015.
% Author: Georges T. Nehmetallah
% modified by Parsa Omidi

%%
% lambda: optical wavelength [m]
[Ny,Nx] = size(Interferogram);

BeamDiameter = 10e-3; % [m]
k = 2*pi / (lambda./10^-3); % optical wavevector
fX = (-Ny/2 : Ny/2-1) / (Ny*PixelSize);
% observation-plane coordinates
[XX, YY] = meshgrid(lambda./10^-3 * FocusDistance * fX);
% analytic result
planeWave = exp(1i*k/(2*FocusDistance)*(XX.^2+YY.^2))/(1i*(lambda./10^-3)*FocusDistance) * BeamDiameter^2*pi/4 ...
    .* jinc(BeamDiameter*sqrt(XX.^2+YY.^2)/((lambda./10^-3)*FocusDistance));

rec_image = (IFT2Dc(FT2Dc(Interferogram).*planeWave));
%         Mag_RecImage=abs(RecImage).^2;
%         im (Mag_RecImage)

end