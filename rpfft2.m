function S = rpfft2(X, options)
    arguments
        X (:,:) {mustBeNumeric}
        options.subsample logical
    end
    if any(isnan(X)); error("Input contains NaN's"); end
    m = mean(X, 'all');
    Y = X - m;
    F = fftn(Y);
    sz = size(Y);
    fs = {fftfreqs(sz(1), [], true), fftfreqs(sz(2), [], true)};
    phi = rand(sz).*2*pi;
    flop = @(phi) -flip(flip(phi, 1), 2);
    phi(fs{1} < 0,  fs{2} == 0) = flop(phi(fs{1} > 0, fs{2} == 0));
    phi(fs{1} == 0, fs{2} < 0)  = flop(phi(fs{1} == 0, fs{2} > 0));
    phi(fs{1} < 0,  fs{2} < 0)  = flop(phi(fs{1} > 0, fs{2} > 0));
    phi(fs{1} < 0,  fs{2} > 0)  = flop(phi(fs{1} > 0, fs{2} < 0));

    % These ones are specifically for even dimensions
    phi(isnan(fs{1}), fs{2} < 0) = flop(phi(isnan(fs{1}), fs{2} > 0));
    phi(fs{1} < 0, isnan(fs{2})) = flop(phi(fs{1} > 0, isnan(fs{2})));
    phi(isnan(fs{1}), isnan(fs{2})) = 0.0;
    phi(isnan(fs{1}), fs{2} == 0) = 0.0;
    phi(fs{1} == 0, isnan(fs{2})) = 0.0;
    
    S = ifftn(abs(F).*exp(1i*phi));
end
