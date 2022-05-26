function S = rpfft3(X, options)
% RPFFT3 Compute one random fourier surrogate of the input array X by taking the three dimensional Fourier transform, assigning random phases to all frequency components, then inverse fourier transforming the result.
% The symmetry condition is F(i, j, k) = F(-i, -j, -k) for i, j, k > 0 (?)
arguments
    X (:,:,:) {mustBeNumeric}
    options.subsample logical
end

if any(isnan(X)); error("Input contains NaN's"); end
if any(arrayfun(@(x) ~mod(x, 2), size(X))); error("Input array must be odd-sized in all dimensions"); end

m = mean(X, 'all');
Y = X - m;
F = fftn(Y);
sz = size(Y);
fs = arrayfun(@fftfreqs, sz, 'un', 0);
phi = rand(sz).*2*pi;

function out = idx(i, j, k)
    out = cellfun(@findmin, {fs{1} - i, fs{2} - j, fs{3} - k}, 'un', 0);
    out = cellfun(@(x) x(1), out, 'un', 0);
end

for i = fs{1}; for j = fs{2}; for k = fs{3};
    a = idx(i, j, k);
    b = idx(-i, -j, -k);
    phi(a{:}) = -phi(b{:});
end; end; end

S = ifftn(abs(F).*exp(1i*phi));
S = S + m
end


function i = findmin(x)
    [~, i] = min(abs(x));
end