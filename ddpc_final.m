format long
EPS_SYS = 0.002;
EPS = EPS_SYS/2;
EPS_DATA = EPS/10000;

check_algo_for_random_good_sys (EPS, EPS_SYS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% REPEATED OPEN-LOOP OC %%%%%%%%%%%%%%%%%%%%%%%%%
function [isFeasible, uopt, yopt] = control_over_T (sys, x_init, Hu, Hy, u_init, y_init, GG, gg, umax, umin, EPS)
    isFeasible = true;
    t_init = length(u_init) + 1; %when we start control
    T = length(Hu(:, 1)); %prediction horizon
    up = u_init;
    yp = y_init;
    GGtau = GG;
    ggtau = gg;
    for tau = t_init:T        
        Up = Hu(1:tau-1, :);
        Uf = Hu(tau:end, :);
        Yp = Hy(1:tau-1, :);
        Yf = Hy(tau:end, :);
        
        %-2 inequality constraints each time
        if tau > t_init
            GGtau = GGtau(1:end-2, 1:end - 1);
            ggtau = ggtau(1:end-2);
        end
        
        [isFeasibleatTau, uf, yf] = control_at_tau (sys, x_init, Up, Yp, Uf, Yf, up, yp, EPS, GGtau, ggtau, umax, umin);
        if isFeasibleatTau == false
            isFeasible = false;
            break;
        end
        diary on
        uf
        yf
        up' * up + uf' * uf %cost function
        diary off
        up = vertcat(up, uf(1));
        yp = vertcat(yp, yf(1) + EPS * rand);
    end
    uopt = up;
    [yopt, ~, ~] = lsim(sys, uopt, [], x_init);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPEN-LOOP OC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isFeasible, uf, yf] = control_at_tau (sys, x_init, Up, Yp, Uf, Yf, up, yp, EPS, GG, gg, umax, umin)
   [alphaDet, chiDet, solutionExistsDet] = ...
       estimation_Det (Up, Yp, Uf, Yf, up, yp, EPS, GG);
   %optimal input
   [uf, exitflagDet] = ...
       control_Det (Up, Yp, Uf, Yf, GG, gg, umax, umin, chiDet);

    if exitflagDet == 1 && solutionExistsDet == 1 %Function converged to the solution x.
        isFeasible = true;
    else
    %error('control at tau infeasible! lost feasibility');
        isFeasible = false;
    end
    
    [yf, ~, ~] = lsim(sys, vertcat(up, uf), [], x_init);
    yf = yf(length(yp) + 1 : end);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATION OF A SYSTEM WITH FEASIBLE SOLUTION %%%%%%%%%%%%%%%
function [sys, x_init, u_init, y_init_clear, Hu, Hy, GG, gg, umax, umin] = generate_a_nice_sys (EPS)
for attempt = 1:100
    %%%%%
    %%%LTI sys
    %%%%%%%%
    a = randi([-15, 15], 10, 1)/10;
    A_SYS = [0 1 0; 0 0 1; a(1) a(2) a(3)];
    B_SYS = [0; 0; 1];
    C_SYS = [a(4) a(5) a(6)]; 
    D_SYS = [a(7)];
    X0_SYS = transpose([a(8)/10 a(9)/10 a(10)/10]); % initial condition
    sys = ss(A_SYS, B_SYS, C_SYS, D_SYS, -1);
    
    %%%%%%%%%%%%%%
    %Input and output constraints
    %%%%%%%%%%%%%%
    % -1 <= y <=1
    G = [1; -1];
    g = [1;1]; 
    % -1 <= u <= 1
    umax = 1;
    umin = -1;

    %%%%%%%%%%
    %DATA GENERATION
    %%%%%%%%%%
    %TODO: make it a function, not a global variable
    %global n; 
    n = 3; %system dimension
    %global L ;
    L = n + 8; %at first it was n + 2, prediction horizon
    
    %generating an INPUT signal u, uniformly distributed in (0,1)
    u = rand(2*(n+L),1);
    %creating a Hankel matrix (n + L) x (n + L + 1)of order n+L for u 
    H_u_check = hankel(u(1:n+L), u(n+L:end));
    %check for persistency of excitation 
    rank(H_u_check)
    %TODO: add a ref on lin ind prob for random matrices
    %hooray, full row rank!!

    %generate OUTPUTS y
    [y, ~, ~] = lsim(sys, u, 0:1:length(u)-1, X0_SYS);
    %Hankel matrix for y, u
    Hy = hankel(y(1:L), y(L:end));
    Hu = hankel(u(1:L), u(L:end));

    %we start control at tau >= n, part of trajectory before n is  generated here
    tau = n;%just for now, current position in time
    x_init = transpose([rand rand rand]); % our actual system initial value [0 0 1] worked
    u_init = rand(tau,1); %actual system inputs, uniformly distributed in (0, 1)
    [y_init, ~, ~] = lsim(sys, u_init, [], x_init);

    %add some noise to y_init
    y_init_clear = y_init;
    y_init = y_init + EPS * rand(length(y_init), 1);

    %%%%%%%%%%%
    %%%MATRICES at time TAU
    %%%%%%%%%%%
    %Up, Yp, Uf, Yf
    Up = Hu(1:tau, :);
    Uf = Hu(tau+1:end, :);
    Yp = Hy(1:tau, :);
    Yf = Hy(tau+1:end, :);

    %GG & gg - block matrices formed from G(t), g(t) 
    %TODO: general case
    gg = ones(2*(L-tau),1);
    GG = zeros(2*(L -tau), L-tau);
    for i = 1:L-tau
        GG(2*i-1,i) = 1;
        GG(2*i,i) = -1; 
    end

    %worst-case estimation
    [alphaDet, chiDet, solutionExistsDet] = ...
        estimation_Det (Up, Yp, Uf, Yf, u_init, y_init_clear, EPS, GG);
    %optimal input
    [uOptDet, exitflagDet] = ...
        control_Det (Up, Yp, Uf, Yf, GG, gg, umax, umin, chiDet);

    if exitflagDet == 1 && solutionExistsDet == 1 %Function converged to the solution x.
        break;
    end

end %end of attempt to generate
end



%%%%%%%%%%%%%%%%%%  ESTIMATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alphaDet, chiDet, allexitflagsEqualOne] = estimation_Det (Up, Yp, Uf, Yf, u_init, y_init, EPS, GG)
%DETERMINANT system: pre-collected data is NOISE-FREE
%Linear Program for ESTIMATION
allexitflagsEqualOne = true;

tau = length(u_init);
L = tau + length(Uf(:,1));
alphaDet = cell(2*(L-tau), 1); %correspond to worst case
chiDet = zeros(2*(L-tau),1); %for each row of GG   

%precalculate MATRICES for linprog
A = vertcat(Yp, -Yp);
b = vertcat(y_init + EPS.*ones(length(y_init),1), ...
    -y_init+EPS.*ones(length(y_init),1) );
Aeq = vertcat(Up, Uf);
beq = vertcat(u_init, zeros(L-length(u_init), 1));

options = optimoptions('linprog', 'Display','off', 'MaxIter', 5000);
for i = 1:length(chiDet)
    [x, fval, exitflag, ~] = ...
        linprog(-GG(i,:)*Yf, ... %-GG(i,:) instead of GG(i, :) for a max instead of min
            A, b, Aeq, beq, ...
                [], [], options);
            
    if exitflag == 1
        alphaDet{i} = x;
        chiDet(i) = fval;
    else
        allexitflagsEqualOne = false;
    end

end
chiDet = -chiDet; % minus before chi to contain a max instead of min 
 
end



%%%%%%%%%%%%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uOptDet, exitflag] = control_Det (Up, Yp, Uf, Yf, GG, gg, umax, umin, chiDet) 
%Quadratic Problem for CONTROL
%alphaC stays for control
[alphaCdet,~, exitflag, ~] = quadprog(transpose(Uf)*Uf, [], vertcat(GG*Yf, Uf, -Uf), ... %%for -1< = Uf*alphaCdet <= 1
    vertcat(gg-chiDet, umax * ones(length(Uf(:,1)), 1), -umin * ones(length(Uf(:, 1)), 1)),...
    vertcat(Up, Yp), zeros(size(Up,1)+size(Yp,1),1));

%optimal input
uOptDet = Uf*alphaCdet;
 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO CHECK DIFFERENT STUFF %%%%%%%%%%%%%%%%
function check_algo_for_random_good_sys(EPS, EPS_SYS)
totalFeasible = 0;

[sys, x_init, u_init, y_init_clear, Hu, Hy, GG, gg, umax, umin] = generate_a_nice_sys (EPS_SYS);
for i = 1:100
y_init = y_init_clear + EPS * rand(length(y_init_clear), 1);
[isFeasible, uopt, yopt] = control_over_T (sys, x_init, Hu, Hy, u_init, y_init, GG, gg, umax, umin, EPS);
totalFeasible = totalFeasible + isFeasible;
diary on
sys
x_init
u_init
y_init

isFeasible
uopt
yopt
diary off

if isFeasible
    figure
    stem(length(u_init):length(uopt)-1, uopt(length(u_init) + 1:end), 'DisplayName', 'u');
    legend('u')
    figure
    stem(length(y_init):length(yopt)-1, yopt(length(y_init)+1:end), 'DisplayName', 'y');
    legend('y')
end

end

diary on
totalFeasible
diary off
end



function [totalFeasible, totalFeasibleButNotSucceed] = check_how_noise_affects (EPS, EPS_SYS, EPS_DATA)
[sys, x_init, u_init, y_init_clear, Hu, Hy, GG, gg, umax, umin] = generate_a_nice_sys (EPS_SYS);
totalFeasible = 0;
totalFeasibleButNotSucceed = 0;
y_init = y_init_clear + EPS * rand(length(y_init_clear), 1);

for i = 1:1000
%here we generate different noises
%TODO derive n, L here from Hu
n = 3;
L = n + 8;
ksi = EPS_DATA * rand(2*(n+L),1);
HyNoisy = hankel(ksi(1:L), ksi(L:end))+Hy;

[isFeasible, uopt, yopt] = control_over_T (sys, x_init, Hu, HyNoisy, u_init, y_init, GG, gg, umax, umin, EPS);
totalFeasible = totalFeasible + isFeasible;
if isFeasible & max(abs(yopt(length(y_init) + 1:end))) > 1
    totalFeasibleButNotSucceed = totalFeasibleButNotSucceed + 1;
    diary on
    uopt
    yopt
    diary off 
end
end

diary on
totalFeasible
totalFeasibleButNotSucceed
diary off

end