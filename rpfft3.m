function S = rpfft3(X, options)
% RPFFT3 Compute one random fourier surrogate of the input array X by taking the three dimensional Fourier transform, assigning random phases to all frequency components, then inverse fourier transforming the result.
% The symmetry condition is F(i, j, k) = F(-i, -j, -k) for i, j, k > 0 (?)
arguments
    X (:,:,:) {mustBeNumeric}
    options.subsample logical
end

if any(isnan(X)); error("Input contains NaN's"); end

m = mean(X, 'all');
Y = X - m;
F = fftn(Y);
sz = size(Y);
fs = arrayfun(@fftfreqs, sz, 'un', 0);
phi = rand(sz).*2*pi;
% flop = @(phi) -flip(flip(phi, 1), 2);


function out = idx(i, j, k)
    if i > max(fs{1}); i = -i; end % We have the unmatched negative frequency component
    if j > max(fs{2}); j = -j; end
    if k > max(fs{3}); k = -k; end
    out = cellfun(@find, {fs{1} == i, fs{2} == j, fs{3} == k}, 'un', 0);
    out = cellfun(@(x) x(1), out, 'un', 0);
end


% phi(fs{1} < 0,  fs{2} == 0, fs{3} == 0) = flop(phi(fs{1} > 0, fs{2} == 0, fs{3} == 0));
% phi(fs{1} == 0,  fs{2} < 0, fs{3} == 0) = flop(phi(fs{1} == 0, fs{2} > 0, fs{3} == 0));
% phi(fs{1} == 0,  fs{2} == 0, fs{3} < 0) = flop(phi(fs{1} == 0, fs{2} == 0, fs{3} > 0));


% phi(fs{1} < 0,  fs{2} < 0, fs{3} < 0)  = flop(phi(fs{1} > 0, fs{2} > 0, fs{3} > 0));

% phi(fs{1} < 0,  fs{2} > 0)  = flop(phi(fs{1} > 0, fs{2} < 0));
for i = fs{1}; for j = fs{2}; for k = fs{3};
    a = idx(i, j, k);
    b = idx(-i, -j, -k);
    phi(a{:}) = -phi(b{:});
end; end; end


disp(fftshift(phi))


% These ones are specifically for even dimensions
% phi(isnan(fs{1}), fs{2} < 0) = flop(phi(isnan(fs{1}), fs{2} > 0));
% phi(fs{1} < 0, isnan(fs{2})) = flop(phi(fs{1} > 0, isnan(fs{2})));
% phi(isnan(fs{1}), isnan(fs{2})) = 0.0;
% phi(isnan(fs{1}), fs{2} == 0) = 0.0;
% phi(fs{1} == 0, isnan(fs{2})) = 0.0;
% disp(fftshift(phi))


S = ifftn(abs(F).*exp(1i*phi));
end
