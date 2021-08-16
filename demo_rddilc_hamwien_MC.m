%+----------------------------------------------------------------+%
%                                                                  %
% Demo of implementing the robust data driven ILC alogrithm for a  %
% MIMO Hammerstein system                                          %
% Author: Jianfei Dong (c), Apr. 26, 2021                          %
%                                                                  %
%+----------------------------------------------------------------+%

close all
clear all
clc

%% define the true plant
% define true para
Num = { [0 -0.0157], [0 -0.0047]; [0 -0.0201], [0 -0.0302] };
Den = { [1 -0.9522], [1 0.02754]; [1 -0.9060], [1 -0.9991] };
Ts = 1;
m = tf(Num, Den, Ts, 'Variable', 'z^-1');
ssmdl = idss(m);
ulim = [-1, 1];

% the true predictor parameters
A = ssmdl.a; B = ssmdl.b; C = ssmdl.c; D = ssmdl.d; K = ssmdl.k;
pole(ssmdl)

% dimensions
nx = size(A,1);
nu = size(B,2);
ny = size(C,1);


%% tuning parameters for robust_dd_ilc
sgm = 0.02;
s = 15; 
f = s+1; 
q = 1e2; 
G = 1;
nG = length(G);
H = [1 2];                   % define the polynomial approximator of the Hammerstein model
nH = length(H);
S = 1e-2*eye((f-1)*nu*nH);   % weighting on inputs
Rn = 1e0*eye((f-1)*nu*nH);   % weighting of input changes in nominal case
w = 0.1;
ulim = [-1, 1];

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

%% MC parameters

mc_num = 500;
MC_robust  = [];
MC_nominal = [];
THETA = [];
SIGMA = [];
NORM_Psi = []; 
NORM_DXi = [];

% PRBS excitation signals
Nt = 200;
ut = idinput([Nt, nu],'prbs',[0,1/1],[-1 1])';

%% MC simulations
disp(['<<<<<<<<<<<<<<<<<<<<<<<< running Monte Carlo Simulations >>>>>>>>>>>>>>>>>>>>>>>>']);

for mc = 1:mc_num
    
    if mod(mc, 10) == 0
        disp(['<<<<<<<<<<<<<<<<<<<<<<<< ' num2str(mc) ' Monte Carlo Simulations Completed >>>>>>>>>>>>>>>>>>>>>>>>']);
    end

    %% simulate the model to get I/O data
    x = zeros(nx,Nt);
    y = zeros(ny,Nt);
    v = zeros(nu,Nt);
    e  = sgm*randn(ny, Nt);

    for kk = 2:Nt
        x(:,kk) = A*x(:,kk-1) + B*v(:,kk-1);
        y(:,kk) = C*x(:,kk) + e(:,kk);
        v(:,kk) = hamm_bidistill(ut(:,kk)); 
    end

    %% estimate the parameters and their errors from data
    [Qcal, L, Psi, Se, Theta0, Sigma_Theta, beta, alpha, Q_nominal, L_nominal, S_nominal, R_nominal] = robust_dd_ilc_hammwien(ut, y, H, G, s, f, [], q, S, Rn, w);
    Beta = kron(eye(f), beta);
    Alpha = kron(eye(f-1), alpha);
    THETA = [THETA; Theta0];
    SIGMA = [SIGMA; Sigma_Theta];
    NORM_Psi = [NORM_Psi; norm(Psi)];
    NORM_DXi = [NORM_DXi; norm(Sigma_Theta)];

    %% test the robust ILC
    it_num = 100; 
    N = it_num * f; 
    kk = (0:f-1)';
    yr = [0.02; 0.98];
    Rk = kron( ones(f,1), yr );
    r  = kron( ones(1,N), yr );
    e  = sgm*randn(ny,N);
    x  = zeros(nx,N); 
    y  = zeros(ny,N);
    u  = zeros(nu,N);
    Uk = ones(nu*ni, 1);
    Uk_t = ones(nu*ni*nH, 1);
    k = 0;
    Yk = [];
    iter_rmse = [];
    
    for t = 1:N-1

        jj = t - k*ni; 

        % simulate the dynamics with the previously learnt Uk
        u(:,t) = Uk( (jj-1)*nu+1 : jj*nu );
        v(:,t) = hamm_bidistill(u(:,t));
        x(:,t+1) = A*x(:,t) + B*v(:,t);
        y(:,t)   = C*x(:,t) + e(:,t);
        Yk = [Yk; y(:,t)];

        % start the next iteration
        if mod(t, ni) == 0
            k = k + 1;
            y(:,t+1) = C*x(:,t+1) + e(:,t+1);
            Yk   = [Yk; y(:,t+1)];
            Ek = Rk - Yk;
            Uk_t = Qcal * Uk_t + L * Ek; 
            iter_rmse = [iter_rmse; sqrt(mean(Ek.*Ek))];
            Yk   = [];
        end

        Uk = solve_hamm(Uk_t, H, nu, ni);
        
       %% uncomment this for bounded inputs, otherwise, no bound on inputs
        % Uk = satu(Uk, ulim);
    end
    
    MC_robust = [MC_robust; iter_rmse(end)];    

    %% test the nominal ILC
    e2 = sgm*randn(ny,N); 
    x2 = zeros(nx,N); 
    y2 = zeros(ny,N);
    u2  = zeros(nu,N);
    Uk = ones(nu*ni, 1);
    Uk_t = ones(nu*ni*length(H), 1);
    k = 0;
    Yk = [];
    iter_rmse2 = [];

    for t = 1:N-1

        jj = t - k*ni; % rem(t, ni);

        % simulate the dynamics with the previously learnt Uk
        u2(:,t) = Uk( (jj-1)*nu+1 : jj*nu );
        v2(:,t) = hamm_bidistill(u2(:,t));
        x2(:,t+1) = A*x2(:,t) + B*v2(:,t);
        y2(:,t)   = C*x2(:,t) + e2(:,t);
        Yk = [Yk; y2(:,t)];

        % start the next iteration
        if mod(t, ni) == 0
            k = k + 1;
            y2(:,t+1) = C*x2(:,t+1) + e2(:,t+1);
            Yk = [Yk; y2(:,t+1)];
            Ek = Rk - Yk;
            Uk_t = Q_nominal * Uk_t + L_nominal * Ek; 
            iter_rmse2 = [iter_rmse2; sqrt(mean(Ek.*Ek))];
            Yk = [];
        end

        Uk = solve_hamm(Uk_t, H, nu, ni);
       %% uncomment this for bounded inputs, otherwise, no bound on inputs
        % Uk = satu(Uk, ulim);
    end
    
    MC_nominal = [MC_nominal; iter_rmse2(end)];

end

%% treat NaN as 1E3
MC_nominal( find( isnan(MC_nominal)==1 ) ) = 1E3;
MC_robust( find( isnan(MC_robust)==1 ) ) = 1E3;

% % rob_good_perf_num = length(find(MC_robust < 1e-2));
% % rob_bad_perf_num = length(find(MC_robust >= 1e-2 & MC_robust < 1e-1));
% % rob_unstable_num = length(find(MC_robust >= 1e-1 ));
% % disp(['robust ILC, percentage of MC simulations with MSE<1e-2 is ' num2str(rob_good_perf_num/mc_num*100) '%']);
% % disp(['robust ILC, percentage of MC simulations with 1e-2<MSE<1e-1 is ' num2str(rob_bad_perf_num/mc_num*100) '%']);
% % disp(['robust ILC, percentage of MC simulations with MSE>1e-1 is ' num2str(rob_unstable_num/mc_num*100) '%']);
% % 
% % nom_good_perf_num = length(find(MC_nominal < 1e-2));
% % nom_bad_perf_num = length(find(MC_nominal >= 1e-2 & MC_nominal < 1e-1));
% % nom_unstable_num = length(find(MC_nominal >= 1e-1 ));
% % disp(['nominal ILC, percentage of MC simulations with MSE<1e-2 is ' num2str(nom_good_perf_num/mc_num*100) '%']);
% % disp(['nominal ILC, percentage of MC simulations with 1e-2<MSE<1e-1 is ' num2str(nom_bad_perf_num/mc_num*100) '%']);
% % disp(['nominal ILC, percentage of MC simulations with MSE>1e-1 is ' num2str(nom_unstable_num/mc_num*100) '%']);

disp(['Binary Dist: robust+ulim ave = ' num2str(mean(MC_robust)) ', std = ' num2str(std(MC_robust)) ...
      '; nominal+ulim ave = ' num2str(mean(MC_nominal)) ', std = ' num2str(std(MC_nominal))]);
disp(['corr of robust tracking error with Psi = ' num2str(corr(NORM_Psi,MC_robust) )] );
disp(['corr of nominal tracking error with Psi = ' num2str(corr(NORM_Psi,MC_nominal) )] );

%% show results

% box plot
figure
boxplot( [log10(MC_robust), log10(MC_nominal)], 'Labels', {'robust ILC','nominal ILC'} )
ylabel('log_{10}(RMSE)')


% plot RMSE against Psi and fitting to line
result_NORM_Psi = NORM_Psi;
temp = inv([result_NORM_Psi, ones(size(result_NORM_Psi))]'*[result_NORM_Psi, ones(size(result_NORM_Psi))]) * [result_NORM_Psi, ones(size(result_NORM_Psi))]';
result_MC_robust = MC_robust;
a = temp * result_MC_robust;
rob_hat = [result_NORM_Psi, ones(size(result_NORM_Psi))] * a;
figure
plot(result_NORM_Psi, result_MC_robust, 'bo', 'MarkerSize', 3)
hold on
plot(result_NORM_Psi, rob_hat, 'k-')
xlabel('||\Psi||_2')
ylabel('RMSE')



