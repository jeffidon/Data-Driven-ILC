%+----------------------------------------------------------------+%
%                                                                  %
% Demo of implementing the robust data driven ILC alogrithm for a
% SISO H-W system                                                  %
% Author: Jianfei Dong (c), Apr. 26, 2021                          %
%                                                                  %
%+----------------------------------------------------------------+%

close all
clear 
clc

%% define the true plant
% define true para
A = [1 -1.5 0.7];
B = [0 0.5 0.25];
Ts = 1;
m = idpoly(A,B,1,1,1,1,Ts);
ssmdl = idss(m);

% the true predictor parameters
A = ssmdl.a; B = ssmdl.b; C = ssmdl.c; D = ssmdl.d; K = ssmdl.k;

% dimensions
nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

%% tuning parameters for robust_dd_ilc
sgm = 0.01;
on  = 0.001;
s = 39; 
f = s+1; 
q = 1; 
G = [1 2 3];                 % i.e. G(y) = a1*y + a2*y^2  + a3*y^3;
nG = length(G);
H = [1 2 3 4];               % define the polynomial approximator of the Hammerstein model
nH = length(H);
S = 1e-3*eye((f-1)*nu*nH);   % weighting on inputs
Rn = 1*eye((f-1)*nu*nH);     % weighting of input changes in nominal case
r_m = 5e-2;                  % training input amplitude  
w = 10;

%% simulate the ILC controller

if ~isempty(find(D~=0))                                                    
    % with feedthrough, Uk is f*nu dim
    ni = f;
else
    % without feedthrough, Uk is (f-1)*nu dim, the iteration duration
    % shall be f-1, instead of f, otherwise, u(f) will never be set,
    % and x(f+1) can never be computed!!!
    ni = f-1;
end


%% tuning parameters for ddo ilc 
T = f;
U0   = 0.5 * ones(T,1);
Phi0 = 0.5 * toeplitz_lower(ones(T,1),1,1,T);
para = ones(6,1);
para(1) = T; 
para(2) = 1e-1; % rho, step size on u 
para(3) = 1;    % eta, step size on Phi 0~2
para(4) = 1e0;  % lambda, penalize u
para(5) = 1e0;  % mu, penalize Phi
para(6) = 1e-3; % epsilon, resetting bound

%% MC parameters
mc_num = 200; 
it_num = 200;   
N = it_num * f; 
kk = (0:f-1)';
Rk = 10 + 2*(1+sin(2*pi*kk/f));  
r = kron( ones(it_num,1), Rk  );
MC_robust  = [];
MC_nominal = [];
MC_ddoilc  = [];
CB  = [];
CAB = [];

%% simulate the model to get I/O data
Nt = 200;        
t = 0:Ts:Nt*Ts; 
x = zeros(nx,Nt);
y = zeros(ny,Nt);
v = zeros(nu,Nt);
z = zeros(ny,Nt);
e = sgm*randn(ny, Nt);
e = filter(1, [1 -0.9], e); 
c = on*randn(ny, Nt);
u = 0.1 + r_m*randn(nu,Nt);     % 1e8*randn(nu,N);

for kk = 2:Nt
    x(:,kk) = A*x(:,kk-1) + B*v(:,kk-1);
    z(:,kk) = C*x(:,kk) + e(:,kk);
    v(:,kk) = hamm_siso(u(:,kk)); 
    y(:,kk) = wien_siso(z(:,kk)) + c(:,kk);
end

% save the I/Os
training_IO = [u;y];

%% estimate the parameters and their errors from data
disp(['<<<<<<<<<<<<<<<<<<<<<<<< running the main algorithm >>>>>>>>>>>>>>>>>>>>>>>>']);
[Qcal, L, Psi, Se, Theta0, Sigma_Theta, beta, alpha, Q_nominal, L_nominal, S_nominal, R_nominal] = robust_dd_ilc_hammwien(u, y, H, G, s, f, [], q, S, Rn, w);
Beta = kron(eye(f), beta);
Alpha = kron(eye(f-1), alpha);
norm_Psi = norm(Psi);
disp( sprintf('2 norm of Psi is %f', norm_Psi ) );
norm_DXi = norm(Sigma_Theta);
disp( sprintf('2 norm of Delta Xi is %f', norm_DXi ) );

%% MC simulations
disp(['<<<<<<<<<<<<<<<<<<<<<<<< running Monte Carlo Simulations >>>>>>>>>>>>>>>>>>>>>>>>']);

for mc = 1:mc_num
    
    if mod(mc, 10) == 0
        disp(['<<<<<<<<<<<<<<<<<<<<<<<< ' num2str(mc) ' Monte Carlo Simulations Completed >>>>>>>>>>>>>>>>>>>>>>>>']);
    end

    %% test the robust ILC
    e = sgm*randn(ny,N);
    e = filter(1, [1 -0.9], e); 
    c = on*randn(ny, N);
    x = zeros(nx,N); x(:,1) = 0.1*ones(nx,1);
    z  = zeros(ny,N); z(:,1) = C*x(:,1);
    y  = zeros(ny,N);
    u  = zeros(nu,N);
    Uk = ones(nu*ni, 1);
    Uk_t = ones(nu*ni*length(H), 1);
    k = 0;
    Yk = [];
    iter_rmse = [];

    for t = 1:N-1

        jj = t - k*ni; % rem(t, ni);

        % simulate the dynamics with the previously learnt Uk
        u(:,t) = Uk( (jj-1)*nu+1 : jj*nu );
        v(:,t) = hamm_siso(u(:,t));
        x(:,t+1) = A*x(:,t) + B*v(:,t);
        z(:,t)   = C*x(:,t) + e(:,t);
        y(:,t) = wien_siso(z(:,t)) + c(:,t);
        Yk = [Yk; y(:,t)];

        % start the next iteration
        if mod(t, ni) == 0
            k = k + 1;
            z(:,t+1) = C*x(:,t+1) + e(:,t+1);
            y(:,t+1) = wien_siso(z(:,t+1)) + c(:,t+1);
            Yk   = [Yk; y(:,t+1)];
            lifted_Yk = wiener_lifted(G, Yk, ny, f);
            lifted_Rk = wiener_lifted(G, Rk, ny, f);
            lifted_Ek = lifted_Rk - lifted_Yk;
            Ek = Beta * (lifted_Rk - lifted_Yk);
            Uk_t = Qcal * Uk_t + L * Ek; 
            ek = Rk - Yk;
            iter_rmse = [iter_rmse; sqrt(mean(ek.*ek))];
            Yk   = [];
        end

        Uk = solve_hamm(Uk_t, H, nu, ni); 

    end
    
    MC_robust = [MC_robust; iter_rmse(end)];    

    %% test the nominal ILC
    e2 = sgm*randn(ny,N); e2 = filter(1, [1 -0.9], e2); 
    c2 = on*randn(ny, N);
    x2 = zeros(nx,N); x2(:,1) = 0.1*ones(nx,1);
    z2 = zeros(ny,N); z2(:,1) = C*x2(:,1);
    y2 = zeros(ny,N);
    u2  = zeros(nu,N);
    Uk = ones(nu*ni, 1);
    Uk_t = ones(nu*ni*length(H), 1);
    k = 0;
    Yk = [];
    iter_rmse2 = [];

    for t = 1:N-1

        jj = t - k*ni; 

        % simulate the dynamics with the previously learnt Uk
        u2(:,t) = Uk( (jj-1)*nu+1 : jj*nu );
        v2(:,t) = hamm_siso(u2(:,t));
        x2(:,t+1) = A*x2(:,t) + B*v2(:,t);
        z2(:,t)   = C*x2(:,t) + e2(:,t);
        y2(:,t) = wien_siso(z2(:,t)) + c2(:,t);
        Yk = [Yk; y2(:,t)];

        % start the next iteration
        if mod(t, ni) == 0
            k = k + 1;
            z2(:,t+1) = C*x2(:,t+1) + e2(:,t+1);
            y2(:,t+1) = wien_siso(z2(:,t+1)) + c2(:,t+1);
            Yk = [Yk; y2(:,t+1)];
            lifted_Yk = wiener_lifted(G, Yk, ny, f);
            lifted_Rk = wiener_lifted(G, Rk, ny, f);
            lifted_Ek   = lifted_Rk - lifted_Yk;
            Ek = Beta * (lifted_Rk - lifted_Yk);
            Uk_t = Q_nominal * Uk_t + L_nominal * Ek; 
            ek = Rk - Yk;
            iter_rmse2 = [iter_rmse2; sqrt(mean(ek.*ek))];
            Yk = [];
        end

        Uk = solve_hamm(Uk_t, H, nu, ni);

    end
    
    MC_nominal = [MC_nominal; iter_rmse2(end)];

end

%% print the results to check the performance
good_lim = 0.5;
bad_lim  = 2;
rob_good_perf_num = length(find(MC_robust < good_lim));
rob_bad_perf_num = length(find(MC_robust >= good_lim & MC_robust < bad_lim));
rob_unstable_num = length(find(MC_robust >= bad_lim ));
disp(['robust ILC, percentage of MC simulations with good performance is ' num2str(rob_good_perf_num/mc_num*100) '%']);
disp(['robust ILC, percentage of MC simulations with bad performance is ' num2str(rob_bad_perf_num/mc_num*100) '%']);
disp(['robust ILC, percentage of MC simulations with instability is ' num2str(rob_unstable_num/mc_num*100) '%']);

nom_good_perf_num = length(find(MC_nominal < good_lim));
nom_bad_perf_num = length(find(MC_nominal >= good_lim & MC_nominal < bad_lim));
nom_unstable_num = length(find(MC_nominal >= bad_lim ));
disp(['nominal ILC, percentage of MC simulations with good performance is ' num2str(nom_good_perf_num/mc_num*100) '%']);
disp(['nominal ILC, percentage of MC simulations with bad performance is ' num2str(nom_bad_perf_num/mc_num*100) '%']);
disp(['nominal ILC, percentage of MC simulations with instability is ' num2str(nom_unstable_num/mc_num*100) '%']);

%% box plots of the MC simulations
figure
boxplot( [log10(MC_robust), log10(MC_nominal)], 'Labels', {'robust ILC','nominal ILC'} )
ylabel('log_{10}(RMSE)')

