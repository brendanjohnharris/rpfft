function S = rpfft2(X, options)
arguments
    X (:,:) {mustBeNumeric}
    options.subsample logical = false
    options.oddfix Integer = 0
end

if any(isnan(X)); error("Input contains NaN's"); end

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

for i = fs{1}; for j = fs{2};
    a = idx(i, j, k);
    b = idx(-i, -j, -k);
    phi(a{:}) = -phi(b{:});
end; end

phi(idx(0, min(fs{2}))) = 0
phi(idx(min(fs{1}), 1)) = 0

S = ifftn(abs(F).*exp(1i*phi));
S = S + m
end


function i = findmin(x)
    [~, i] = min(abs(x));
end