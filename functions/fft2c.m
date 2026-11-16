function y = fft2c(x)
% FFT2C  Centered two-dimensional fast Fourier transform.
%   Y = FFT2C(X) computes the 2-D DFT of X with zero-frequency
%   components shifted to the array center (DC at index N/2+1).
%
%   This convention is consistent throughout the DHD repository:
%   all forward and inverse transforms use fft2c / ifft2c so that
%   spatial and frequency arrays share the same centered coordinate
%   system.
%
%   See also IFFT2C.

y = fftshift(fft2(ifftshift(x)));
end
