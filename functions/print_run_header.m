function print_run_header(title, N, K, qName, q, D, wvl, zName, z, p, nBits, g_AD, sigma_n2)
% PRINT_RUN_HEADER  Standard console header for the illustrative scripts.
%
%   See also PRINT_OBJECT_REPORT.

bar = repmat('=',1,70);
fprintf('\n%s\n %s\n%s\n', bar, title, bar);
fprintf('  N = %d, K = %d, %s = %.3f\n', N, K, qName, q);
fprintf('  D = %.3f m, lambda = %.2e m, %s = %.4f m, p = %.4e m\n', D, wvl, zName, z, p);
fprintf('  %d-bit A/D, g_AD = %.4f DN/pe, sigma_n^2 = %.2f pe^2\n', nBits, g_AD, sigma_n2);
end
