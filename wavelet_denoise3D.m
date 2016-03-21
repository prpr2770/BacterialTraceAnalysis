function y = wavelet_denoise3D(x,th)
% x is the noisy 3d image
% th is the percentile hight of the noise to be removed

Wpack = wavedec3(x,8,'sym4');
coefs = Wpack.dec;
m = size(coefs,1);

% finding the maximum wavelet coefficient
maxp = 0;
for i = 1 : m
    temp = max(abs(coefs{i}(:)));
    if temp > maxp
        maxp = temp;
    end
end

% making small wavelet coefficients zero
Cd = cell(m,1);
for i = 1 : m
    Cd{i} = double(coefs{i} .*(abs(coefs{i})>th*maxp));
end

Wpack.dec = Cd;
y = waverec3(Wpack);


