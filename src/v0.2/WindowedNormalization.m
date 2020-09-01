function [Inorm] = WindowedNormalization(I_in, Wy, Wx)
% Function to normalize an image using a sliding window 

% Inputs: I_in    -> input image
%         Wy      -> size of the window (in y direction)
%         Wx      -> size of the window (in x direction)
%         
 
% Outputs: Inorm -> normalized


% created by: Parsa Omidi
% last modified: 2017
[ny,nx] = size (I_in);

%% fixed window 
wx=Wx; %window size
wy=Wy;
pminy=floor(wy/2)+1;
pmaxy=ny-floor(wy/2)-1;               
pminx=floor(wx/2)+1;
pmaxx=nx-floor(wx/2)-1; 
Inorm=zeros(size(I_in));
for m=pminy:1:pmaxy
    for n=pminx:1:pmaxx
        window=I_in(m-floor(wy/2):m-floor(wy/2)+wy-1,n-floor(wx/2):n-floor(wx/2)+wx-1);
        meanI=mean(mean(window));
        Imax=max(max(window));
        Imin=min(min(window));
        Inorm(m,n)=0.5+(I_in(m,n)-meanI)/(Imax-Imin);
    end
end

end

