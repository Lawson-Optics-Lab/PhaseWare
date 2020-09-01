% GoldsteinUnwrap2D implements 2D Goldstein branch cut phase unwrapping algorithm.
%
% References::
% 1. R. M. Goldstein, H. A. Zebken, and C. L. Werner, “Satellite radar interferometry:
%    Two-dimensional phase unwrapping,” Radio Sci., vol. 23, no. 4, pp. 713–720, 1988.
% 2. D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
%    Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.


function IM_unwrapped = UnwrappingGoldstein(IM_phase)

IM_mask=ones(size(IM_phase)); 

%%  Set parameters
max_box_radius=4;                           %Maximum search box radius (pixels)
threshold_std=5;                            %Number of noise standard deviations used for thresholding the magnitude image

%% Unwrap
residue_charge=PhaseResidues(IM_phase, IM_mask);                            %Calculate phase residues
branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
[IM_unwrapped, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping

%% Display results
figure; imagesc(residue_charge), colormap(gray), axis square, axis off, title('Phase residues (charged)');
figure; imagesc(branch_cuts), colormap(gray), axis square, axis off, title('Branch cuts');
figure; imagesc(immultiply(IM_phase,IM_mask)), colormap(gray), axis square, axis off, title('Wrapped phase');
tempmin=min(min(IM_unwrapped));          %This bit is just done to create a pleasing display when a mask is used
temp=(IM_unwrapped==0);
temp_IM=IM_unwrapped;
temp_IM(temp)=tempmin;
figure; imagesc(temp_IM), colormap(gray), axis square, axis off, title('Unwrapped phase');



end



