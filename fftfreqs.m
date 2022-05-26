
function f=fftfreqs(n, dt, symflag)
% FFTFREQS Return the frequencies for an FFT of a signal with length n and sampling period dt.
% If `symflag = true`, return the SYMMETRIC frequencies. This means for an odd signal, all FFT frequencies are returned, byt for an EVEN signal the lowest nagative frequency IS NOT RETURNED because IT DOES NOT HAVE A POSITIVE FREQUENCY COUNTERPART.
% !
    if nargin < 2; dt = 1; end
    if nargin < 3; symflag = false; end
    if isempty(dt); dt = 1; end
    L = n*dt;
    Nq= 1/(2*dt);

    if mod(n, 2)
        f = [0:1/L:(Nq - 1/(2*L)), (-Nq + 1/(2*L)):1/L:-1/L];
    else
        if symflag
            f = [0:1/L:(Nq - 1/L), NaN, (-Nq+1/L):1/L:-1/L];
        else
            f = [0:1/L:(Nq - 1/L), -Nq:1/L:-1/L];
        end
    end
end