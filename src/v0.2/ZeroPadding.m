function [I_out, Diff_y, Diff_x] = ZeroPadding (I_in, NS2)

[ny,nx] = size(I_in);
Diff_y = floor((NS2 - ny)/2);
Diff_x = floor((NS2 - nx)/2);
if rem(ny,2) ~= 0
    I_in = cat(1,I_in, zeros(1, size(I_in, 2)));
end
if rem(nx,2) ~= 0
    I_in = cat(2,I_in, zeros(size(I_in, 1), 1));
end
I_out = padarray(I_in,[Diff_y Diff_x],0,'both');
end