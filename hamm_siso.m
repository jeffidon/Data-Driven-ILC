function y = hamm_siso(x)

y = x ./ sqrt( 0.1 + 0.9*x.^2 );