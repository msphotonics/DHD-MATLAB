function y = ifft2c(x)
% IFFT2C  Centered two-dimensional inverse fast Fourier transform.
%   Y = IFFT2C(X) computes the 2-D inverse DFT of X with
%   zero-frequency components assumed at the array center.
%
%   See also FFT2C.

y = fftshift(ifft2(ifftshift(x)));
end
