function y = hamm_bidistill(x)

y = zeros(2,1);
y(1) = x(1) - 0.0236*x(1)*x(2);
y(2) = x(2) - 0.1823*x(1)*x(2);

return

