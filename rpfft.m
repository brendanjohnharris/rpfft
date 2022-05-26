function S = rpfft(X, options)
% RPFFT Compute one random fourier surrogate of the input vector X.
arguments
    X {mustBeNumeric}
    options.subsample logical
end

if any(isnan(X)); error("Input contains NaN's"); end

m = mean(X, 'all');
Y = X - m;
F = fft(Y);
sz = length(Y);
fs = fftfreqs(sz);
sz = sum(fs > 0);
phi = angle(F);
phi_ = rand(1, sz).*2*pi;
phi(fs > 0) = phi_;
if ~mod(length(Y), 2); phi_ = [phi_, 0]; end
phi(fs < 0) = -flip(phi_); % Maintain symmetry along the temporal dimension
disp(phi)
S = ifft(abs(F).*exp(1i*phi));
S = S + m
end