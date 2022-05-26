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
    
    fs_ = arrayfun(@fftfreqs, sz, 'un', 0);
    fs = cellfun(@oddup, fs_, 'un', 0);
    sz = cellfun(@length, fs, 'un', 1);
    phi = rand(sz).*2*pi;
    
    function out = idx(i, j)
        out = cellfun(@find, {fs{1} == i, fs{2} == j}, 'un', 0);
        out = cellfun(@(x) x(1), out, 'un', 0);
    end
    
    for i = fs{1}; for j = fs{2};
        a = idx(i, j);
        b = idx(-i, -j);
        phi(a{:}) = -phi(b{:});
    end; end;

    a = idx(0, min(fs{1}));
    b = flip(a);
    phi(a{:}) = 0;
    phi(b{:}) = pi;
    a = idx(min(fs{1}), min(fs{2}));
    phi(a{:}) = 0;
    
    idxs = arrayfun(@(i) arrayfun(@(x) ismember(x, fs_{i}), fs{i}), 1:length(sz), 'un', 0);
    phi = phi(idxs{:});
    disp(fftshift(phi))

    S = ifftn(abs(F).*exp(1i*phi));
    S = S + m;
end


function y = oddup(x)
    if mod(length(x), 2)
        y = x;
    else
        y = [x, -min(x)]; % Symmetrize the frequencies
    end
end
