function [Qcal, L, Psi, Se, Theta_hat, Sigma_Theta, beta, alpha, Q_nominal, L_nominal, S_nominal, R_nominal]  = robust_dd_ilc_hammwien(u, y, H, G, s, f, Se, q, S, Rn, w)
         
%+----------------------------------------------------------------+%
%                                                                  %
% data driven ILC robust to the errors in est. open loop predictor %
%              for Hammerstein and Wiener systems                  %
% Inputs: u,y - the I/O data for the data driven design            %
%         the number of rows in u and y equals the data length     %
%         H, G - the orders of the polynomial bases in ascending   %
%             order, without zero order, i.e. H, G = [1, 2, 3,...],%
%             respectively for Hammerstein and Wiener systems      %
%         s and f - the past the future horizon, scalar integer    %
%         Se - noise covariance is known, if [] then estimated here%
%         q, S, Rn, w - weighting factors according definitions    %
% Outputs:Qcal, L- gain and weighting matrices of the robust ILC   %
%         solution                                                 % 
%         Se, Theta_hat, Sigma_Theta, beta, alpha - estimated noise%
%         covariance and Markov parameters                         %
%         Q_nominal, L_nominal, S_nominal, R_nominal - gain and    %
%         weighting matrices of the norminal ILC solution          % 
% Author: Jianfei Dong (c), Apr. 26, 2021, jfeidong@hotmail.com    %
%                                                                  %
%+----------------------------------------------------------------+%


%------------------------------------------------------------------%
% error check and read the inputs %

nin = nargin;
if nin < 10
    error('Wrong number of input arguments');
end
if size(u,1) < size(u,2)
    u = u'; % the data length must be much larger than the input dimension
end
if size(y,1) < size(y,2)
    y = y'; % the data length must be much larger than the output dimension
end
if size(y,1) ~= size(u,1)
    warning('I/O data have different length. The longer one is trancated.');
    if size(y,1) > size(u,1) 
        y = y(size(u,1),:);
    else
        u = u(size(y,1),:);
    end
end

if f ~= s+1
    warning('by default, we set f = s+1!!!'); % the future horizon must satisfy: f <= s+1.
    f = s+1;
end

%----------------------------------------------------------------%
% define the dimensions %
N  = size(y,1);
r  = length(H);
nu = size(u,2) * r; % the lifted nu = nu * length of H = rm
g  = length(G);
ny = size(y,2) * g; % the lifted ny = ny * length of G = ql
k  = s+1;           % the first "prensent" time instant
P  = N-1;           % the number of columns in the data matrices

%----------------------------------------------------------------%
% preprocess iputs by the Hammerstein polynomial fit by lifting  %
% preprocess outputs by the Wiener polynomial fit by lifting     %
Ulft = []; Ylft = [];
for ii = 1:N
    v = [];
    for jj = 1:r
        v = [ v, u(ii, :).^H(jj) ];
    end
    Ulft = [Ulft; v];
    z = [];
    for kk = 1:g
        z = [ z, y(ii, :).^G(kk) ];
    end
    Ylft = [Ylft; z]; 
end

% replace the I/Os with the lifted ones
u = Ulft;
y = Ylft;
clear Ulft Ylft

%----------------------------------------------------------------%
% wrap the data in future and past data matrices %
% by Up and Yf defined as follows, we estimate [CAs-1B ... CB], dim: ny x s*nu
% by default, we set f = s+1

% the past inputs
Up = []; upc = []; 
for ii = 1:N-s
    for jj = 0:s-1
        upc = [upc; u(ii+jj,:)'];
    end
    Up = [Up, upc]; % extended I/O data matrix
    upc = [];
end

% the future
Yf  = [];
for ii = s+1:N     % note that there is no direct feedthrough because Yf is one step ahead of Up!!!
    Yf  = [Yf, y(ii,:)'];
end

%----------------------------------------------------------------%
% closed-loop identification %
L = triu(qr([Up; Yf]',0))'; 
L11 = L(1:s*nu,1:s*nu);
L21 = L(s*nu+1:end,1:s*nu);
cdnum = cond(Up); 
if cdnum > 1e3
    Omega = inv(L11*L11'+1e-2*eye(s*nu));
    Theta0 = L21 * L11' * Omega;
else
    Theta0 = L21/L11;
    Omega = inv(L11*L11');
end

clear L L11 L21

% estimate the noise covariance matrix
if isempty(Se)
    Ef = Yf - Theta0*Up;
    Se = cov(Ef');
end

% compute the error covariance matrix
Sigma_Theta = kron(Omega, Se);

%-----------------------------------------------------------------%
% add SVD decomposition of \beta and \Xi, for more than one bases %
if g > 1
    [UT, ST, VT] = svd(Theta0);
    beta = pinv(UT(:, 1:ny/g)); 
    Theta0 = ST(1:ny/g,1:ny/g) * VT(:,1:ny/g)';
    ny = ny/g; 
    Se = UT(:, 1:ny)' * Se * UT(:, 1:ny);
else
    beta = 1;
end

% estimate alpha from the decomposed Theta 
[nrt, nct] = size(Theta0);
mnH = nct/s;
Theta_breve = [];
for ii = 1:s
    Theta_breve = [ Theta_breve; Theta0( :, end-ii*mnH+1:end-(ii-1)*mnH ) ];
end
[UB, SB, VB] = svd(Theta_breve);
m = nu/r;
alpha = VB(:, 1:m )';
Theta_hat = [];
temp_theta = UB(:, 1:m) * SB(1:m, 1:m);
for ii = 1:s
    Theta_hat = [ Theta_hat, temp_theta(end-ii*ny+1:end-(ii-1)*ny, :) ];
end

% compute the error covariance matrix
Sigma_Theta = kron(Omega, Se);

%% call the routine to compute the robust and norminal ILC gain matrices 
[Qcal, L, Psi, Q_nominal, L_nominal, S_nominal, R_nominal] = ilc_gain(Theta0, Omega, Se, s, f, nu, ny, q, S, w, Rn);

%____________________end of the main program_____________________%

