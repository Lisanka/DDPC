format long

%%%%%%%%%%%%%%%
%LTI SISO system
%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%IN SEARCH FOR CHI. My pain section%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
detSysSolutionExistsCount = 0;
noisyFeasibleTotal = 0;
noisySucceedTotal = 0;
properNoisySucceedTotal = 0;
goodDetsys = 0;
for i=1:100
    diary on
    disp exp#
    i
    diary off
    EPS = 0.005; % upper Bound on noise by default
    EPS_DATA = 0.000000005;
    [detSysSolutionExistsn, noisyFeasible, noisySucceed, chiNoisy] = experiment_session(EPS, EPS_DATA);
    if detSysSolutionExistsn == 10 && noisyFeasible == 10
       % goodDetsys = goodDetsys + 1;
       diary on
       chiNoisy
       diary off      
    end    
end
%}

EPS = 0.005;
[sys, x_init, u_init, y_init, Hu, Hy, GG, gg, umax, umin] = ...
    generate_a_nice_sys (EPS);
uopt = ...
    control(sys, x_init, Hu, Hy, umax, umin, EPS/10, u_init, y_init);
[yopt, ~, ~] = lsim(sys, vertcat(u_init, uopt), 0:1:length(u_init) + length(uopt)-1, x_init);
diary on
yopt
diary off

%all iterations
function uopt = control(sys, X0_SYS, Hu, Hy, umax, umin, EPS, uinit, yinit)
T = length(Hu(:,1)); % total length
tauinit = length(uinit) + 1; % start of control

up = uinit;
yp = yinit;
for tau = tauinit:tauinit %here should stay T, it's just for testing
    Up = Hu(1:tau-1, :);
    Uf = Hu(tau:end, :);
    Yp = Hy(1:tau-1, :);
    Yf = Hy(tau:end, :);
    
    remained = T-tau + 1;
    gg = ones(2 * remained,1);
    GG = zeros(2 * remained, remained);
    for i = 1:remained
    GG(2*i-1,i) = 1;
    GG(2*i,i) = -1; %% changed -1 to 1 just for that riny experiment with noisy data! then returned to initial setup.
    end
    GG
    gg
    [uf, solutionExists, exitflag] = control_at_tau (Up, Uf, Yp, Yf, GG, gg, umax, umin, EPS, up, yp);
    
    %up = vertcat(up, uf(1));
    
    %output
    %[ypp, ~, ~] = lsim(sys, up, 0:1:length(up)-1, X0_SYS);
    %yp = vertcat(yp, ypp(tau) + EPS * rand);
end
uopt = vertcat(up, uf);%up;
solutionExists
exitflag
end



%single iteration
function [uOptDet, solutionExistsDet, exitflagDet] = control_at_tau(Up, Uf, Yp, Yf, GG, gg, umax, umin, EPS, up, yp)
%generates optimal control sequence at time tau
[alphaDet, chiDet, solutionExistsDet] = estimation_Det (Up, Yp, Uf, Yf, up, yp, EPS, GG);
%optimal input
[uOptDet, exitflagDet]= control_Det (Up, Yp, Uf, Yf, GG, gg, umax, umin, chiDet);
end



function [sys, x_init, u_init, y_init, Hu, Hy, GG, gg, umax, umin] = generate_a_nice_sys (EPS)
for attempt = 1:100
        
a = randi([-15, 15], 7, 1)/10; %normrnd(0, 1, 7, 1); % ~N(0, 1)
%a = [-0.7; 0.1; 0.9; -1.4; 1.0; 0; -0.5];

A_SYS = [0 1 0; 0 0 1; a(1) a(2) a(3)];
B_SYS = [0; 0; 1];
C_SYS = [a(4) a(5) a(6)]; 
D_SYS = [a(7)];
X0_SYS = transpose([0 0 0]); % initial condition
sys = ss(A_SYS, B_SYS, C_SYS, D_SYS, -1);

%%%%%%%%%%%%%%
%Input and output constraints
%%%%%%%%%%%%%%
% -1 <= y <=1
G = [1; -1];
%в G сейчас вообще просто мусор, считай и не используем, сразу GG пишем
g = [1;1]; 

% -1 <= u <= 1
umax = 1;
umin = -1;


%TODO: make it a function, not a global variable
global n; 
n = 3; %system dimension
global L ;
L = n + 8; %at first it was n + 2

%%%%%%%%%%
%DATA GENERATION
%%%%%%%%%%
%generating an INPUT signal u, uniformly distributed in (0,1)
u = rand(2*(n+L),1);
%creating a Hankel matrix of order n+L for u 
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

%detSysSolutionExistsn = 0;




%for i = 1:10
%generate NOISE ksi for outputs
ksi = EPS * rand(2*(n+L),1);%normrnd(0, EPS_4 * EPS_4, 2*(n+L), 1); %ksi ~ N(0, EPS^2)
%wgn(2*(n+L),1,0); %rand(2*(n+L),1); it was before, now white noise 0 dB ~ 1 Wt
HyNoisy = hankel(ksi(1:L), ksi(L:end))+Hy;

%we start control at tau >= n, part of trajectory before n is  generated here
tau = n;%just for now, current position in time
x_init = transpose([0 0 1]); % our actual system initial value
u_init = rand(tau,1); %actual system inputs, uniformly distributed in (0, 1)
[y_init, ~, ~] = lsim(sys, u_init, [], x_init);%0:1:length(u_init)-1

%add some noise to y_init
y_init_clear = y_init;
y_init = y_init + EPS * rand(length(y_init), 1);%normrnd(0, EPS*EPS, length(y_init), 1); %~N(0, EPS^2)
%wgn(length(y_init),1,0);%rand(length(y_init), 1); it was before, now white noise
                                          %0 dB ~ 0.001 Wt



%%%%%%%%%%%
%%%MATRICES at time TAU
%%%%%%%%%%%
%Up, Yp, Uf, Yf
Up = Hu(1:tau, :);
Uf = Hu(tau+1:end, :);
Yp = Hy(1:tau, :);
Yf = Hy(tau+1:end, :);
%same but affected with noise
%YpNoisy = HyNoisy(1:tau, :);
%YfNoisy = HyNoisy(tau+1:end, :);

%GG & gg - block matrices formed from G(t), g(t) 
%TODO: general case
gg = ones(2*(L-tau),1);
GG = zeros(2*(L -tau), L-tau);
for i = 1:L-tau
    GG(2*i-1,i) = 1;
    GG(2*i,i) = -1; %% changed -1 to 1 just for that riny experiment with noisy data! then returned to initial setup.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%worst-case estimation
[alphaDet, chiDet, solutionExistsDet] = estimation_Det (Up, Yp, Uf, Yf, u_init, y_init_clear, EPS, GG);
%optimal input
[uOptDet, exitflagDet]= control_Det (Up, Yp, Uf, Yf, GG, gg, umax, umin, chiDet);


if exitflagDet == 1 & solutionExistsDet == 1 %Function converged to the solution x.
    break;
end

end %end of attempt to generate

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT SESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%

function [detSysSolutionExistsn, noisyFeasible, noisySucceed, chiNoisy] = experiment_session(EPS, EPS_DATA)
    %, noisyFeasible, noisySucceed, properNoisySucceed]  = experiment_session(EPS)

%generate a system or set your own
%comment out undesired
a = randi([-15, 15], 7, 1)/10; %normrnd(0, 1, 7, 1); % ~N(0, 1)
%a = [-0.7; 0.1; 0.9; -1.4; 1.0; 0; -0.5];

A_SYS = [0 1 0; 0 0 1; a(1) a(2) a(3)];
B_SYS = [0; 0; 1];
C_SYS = [a(4) a(5) a(6)]; 
D_SYS = [a(7)];
X0_SYS = transpose([0 0 0]); % initial condition
sys = ss(A_SYS, B_SYS, C_SYS, D_SYS, -1);

diary on
sys
diary off

%%%%%%%%%%%%%%
%Input and output constraints
%%%%%%%%%%%%%%
% -1 <= y <=1
G = [1; -1];
%в G сейчас вообще просто мусор, считай и не используем, сразу GG пишем
g = [1;1]; 

% -1 <= u <= 1
umax = 1;
umin = -1;


%TODO: make it a function, not a global variable
global n; 
n = 3; %system dimension
global L ;
L = n + 8; %at first it was n + 2



%%%%%%%%%%
%DATA GENERATION
%%%%%%%%%%
%generating an INPUT signal u, uniformly distributed in (0,1)
u = rand(2*(n+L),1);
%creating a Hankel matrix of order n+L for u 
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

detSysSolutionExistsn = 0;
for i = 1:10
%generate NOISE ksi for outputs
ksi = EPS * rand(2*(n+L),1);%normrnd(0, EPS_4 * EPS_4, 2*(n+L), 1); %ksi ~ N(0, EPS^2)
%wgn(2*(n+L),1,0); %rand(2*(n+L),1); it was before, now white noise 0 dB ~ 1 Wt
HyNoisy = hankel(ksi(1:L), ksi(L:end))+Hy;

%we start control at tau >= n, part of trajectory before n is  generated here
tau = n;%just for now, current position in time
x_init = transpose([0 0 1]); % our actual system initial value
u_init = rand(tau,1); %actual system inputs, uniformly distributed in (0, 1)
[y_init, ~, ~] = lsim(sys, u_init, [], x_init);%0:1:length(u_init)-1

%add some noise to y_init
y_init_clear = y_init;
y_init = y_init + EPS * rand(length(y_init), 1);%normrnd(0, EPS*EPS, length(y_init), 1); %~N(0, EPS^2)
%wgn(length(y_init),1,0);%rand(length(y_init), 1); it was before, now white noise
                                          %0 dB ~ 0.001 Wt



%%%%%%%%%%%
%%%MATRICES at time TAU
%%%%%%%%%%%
%Up, Yp, Uf, Yf
Up = Hu(1:tau, :);
Uf = Hu(tau+1:end, :);
Yp = Hy(1:tau, :);
Yf = Hy(tau+1:end, :);
%same but affected with noise
YpNoisy = HyNoisy(1:tau, :);
YfNoisy = HyNoisy(tau+1:end, :);

%GG & gg - block matrices formed from G(t), g(t) 
%TODO: general case
gg = ones(2*(L-tau),1);
GG = zeros(2*(L -tau), L-tau);
for i = 1:L-tau
    GG(2*i-1,i) = 1;
    GG(2*i,i) = -1; %% changed -1 to 1 just for that riny experiment with noisy data! then returned to initial setup.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%worst-case estimation
[alphaDet, chiDet, solutionExistsDet] = estimation_Det (Up, Yp, Uf, Yf, u_init, y_init, EPS, GG);
%optimal input
[uOptDet, exitflagDet]= control_Det (Up, Yp, Uf, Yf, GG, gg, umax, umin, chiDet);


if exitflagDet == 1 && solutionExistsDet == 1 %Function converged to the solution x.
    detSysSolutionExistsn = detSysSolutionExistsn + 1;
else
    break;
end

    
end %end of experiments for one particular matrix

noisyFeasible = 0;
noisySucceed = 0;
if detSysSolutionExistsn == 10 %we wanna take first good det case
    for j = 1:10
        %generate NOISE ksi for outputs
ksi = EPS_DATA * rand(2*(n+L),1);%normrnd(0, EPS_4 * EPS_4, 2*(n+L), 1); %ksi ~ N(0, EPS^2)
%wgn(2*(n+L),1,0); %rand(2*(n+L),1); it was before, now white noise 0 dB ~ 1 Wt
HyNoisy = hankel(ksi(1:L), ksi(L:end))+Hy;

%we start control at tau >= n, part of trajectory before n is  generated here
tau = n;%just for now, current position in time
x_init = transpose([0 0 1]); % our actual system initial value
u_init = rand(tau,1); %actual system inputs, uniformly distributed in (0, 1)
[y_init, ~, ~] = lsim(sys, u_init, [], x_init);%0:1:length(u_init)-1

%add some noise to y_init
y_init = y_init + EPS * rand(length(y_init), 1);%normrnd(0, EPS*EPS, length(y_init), 1); %~N(0, EPS^2)
%wgn(length(y_init),1,0);%rand(length(y_init), 1); it was before, now white noise
                                          %0 dB ~ 0.001 Wt



%%%%%%%%%%%
%%%MATRICES at time TAU
%%%%%%%%%%%
%Up, Yp, Uf, Yf
Up = Hu(1:tau, :);
Uf = Hu(tau+1:end, :);

%same but affected with noise
YpNoisy = HyNoisy(1:tau, :);
YfNoisy = HyNoisy(tau+1:end, :);

%GG & gg - block matrices formed from G(t), g(t) 
%TODO: general case
gg = ones(2*(L-tau),1);
GG = zeros(2*(L -tau), L-tau);
for i = 1:L-tau
    GG(2*i-1,i) = 1;
    GG(2*i,i) = -1; %% changed -1 to 1 just for that riny experiment with noisy data! then returned to initial setup.
end

        [alphaDetButNoisy, chiDetButNoisy, solutionExistsDetButNoisy] = ...
        estimation_Det (Up, YpNoisy, Uf, YfNoisy, u_init, y_init, EPS_DATA, GG);
        [uOptDetButNoisy, exitflagDetButNoisy]= ...
        control_Det (Up, YpNoisy, Uf, YfNoisy, GG, gg, umax, umin, chiDetButNoisy);
        [yOptDetButNoisy, ~, ~] = lsim(sys, vertcat(u_init, uOptDetButNoisy), [], x_init);
    
        if solutionExistsDetButNoisy == 1 && exitflagDetButNoisy == 1
            noisyFeasible = noisyFeasible + 1;
            if max(abs(yOptDetButNoisy(tau + 1 : end))) < 1
                noisySucceed = noisySucceed  + 1;
            else 
                diary on
                yOptDetButNoisy(tau +1:end)
                diary off
            end
        end
    end
end %end of if detSysSolutionExistsn == 10



if noisyFeasible == 10
    [chiNoisy,~] = estimation_Noisy (Up, Yp, Uf, Yf, u_init, y_init_clear, EPS, EPS_DATA, GG); 
else
    chiNoisy = 0;
end % end of if noisyFeasible == 10


%{
detSysSolutionExistsn = 0;
noisyFeasible = 0;
noisySucceed = 0;
properNoisySucceed = 0;
if exitflagDet == 1 %Function converged to the solution x.
    detSysSolutionExistsn = detSysSolutionExistsn + 1;
    %output for optimal input
    [yOptDet, ~, ~] = lsim(sys, vertcat(u_init, uOptDet), [], x_init);
    %diary on
    %uOptDet
    %yOptDet
    %diary off
    
    
   
    %same for noisy
    [alphaDetButNoisy, chiDetButNoisy, solutionExistsDetButNoisy] = ...
        estimation_Det (Up, YpNoisy, Uf, YfNoisy, u_init, y_init, EPS, GG);
    [uOptDetButNoisy, exitflagDetButNoisy]= ...
        control_Det (Up, YpNoisy, Uf, YfNoisy, GG, gg, umax, umin, chiDetButNoisy);
    [yOptDetButNoisy, ~, ~] = lsim(sys, vertcat(u_init, uOptDetButNoisy), [], x_init);
    
    if solutionExistsDetButNoisy == 1 & exitflagDetButNoisy == 1
        noisyFeasible = noisyFeasible + 1;
        if max(abs(yOptDetButNoisy(tau + 1:end))) < 1
            noisySucceed = noisySucceed + 1;
        else
            %analyze estimate chi, output yOpt
            diary on
            yOptDetButNoisy
            yOptDet
            chiDetButNoisy
            chiDet
            diary off
             %same for noisy but with a proper algo!
           [chiNoisy,~] = estimation_Noisy (Up, YpNoisy, Uf, YfNoisy, u_init, y_init, EPS, GG, alphaDet);
           [chiZeroStateNoisyUp, chiZeroStateNoisyDown] = chi_Zero_State_Noisy (Up, YpNoisy, Uf, YfNoisy, EPS, GG);
           [~, solutionExists] = optC(transpose(chiZeroStateNoisyUp), transpose(chiZeroStateNoisyDown), gg-chiNoisy);
           if solutionExists
           properNoisySucceed = properNoisySucceed + 1;
           end
        end
    end
    

    
    diary on
    ksi
    sys
    a
    EPS
    diary off
    
   
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESULTS ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
[chiNoisy1,~] = estimation_Noisy (Up, YpNoisy, Uf, YfNoisy, u_init, y_init, EPS_1, GG, alphaDet);
[chiNoisy2, ~] = estimation_Noisy (Up, YpNoisy, Uf, YfNoisy, u_init, y_init, EPS_2, GG, alphaDet);
[chiNoisy4, ~] = estimation_Noisy (Up, YpNoisy, Uf, YfNoisy, u_init, y_init, EPS_4, GG, alphaDet);
[chiNoisy, ~] = estimation_Noisy (Up, YpNoisy, Uf, YfNoisy, u_init, y_init, EPS, GG, alphaDet);
figure
stem([0 1 2 3], chiDet, 'DisplayName', 'chi');
hold on
stem([0 1 2 3], chiNoisy, 'DisplayName', 'chiNoisy');
stem([0 1 2 3], chiNoisy1, 'DisplayName', 'chiNoisy1');
stem([0 1 2 3], chiNoisy2, 'DisplayName', 'chiNoisy2');
stem([0 1 2 3], chiNoisy4, 'DisplayName', 'chiNoisy4');
hold off
%}


%{
for i = 1:5
    [chiNoisy, v] = estimation_Noisy (Up, Yp, Uf, Yf, u_init, y_init, EPS, GG, alphaDet);
    diary on
    chiNoisy
    v{1}
    horzcat(Yf, zeros(length(Yf(:,1)), length(ksi)))*v{1}
    diary off
end
%}

%{
[alphaDet, chiDet, solutionExistsDet] = estimation_Det (Up, Yp, Uf, Yf, u_init, y_init, EPS, GG);
[chiNoisy, ~] = estimation_Noisy (Up, YpNoisy, Uf, YfNoisy, u_init, y_init, EPS, GG, alphaDet);
[chiZeroStateNoisyUp, chiZeroStateNoisyDown] = chi_Zero_State_Noisy (Up, YpNoisy, Uf, YfNoisy, EPS, GG);
[~, solutionExists] = optC(transpose(chiZeroStateNoisyUp), transpose(chiZeroStateNoisyDown), gg-chiNoisy);
diary on
solutionExists
diary off
%}



end %end of experiments session





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% FUNCTIONS DECLARATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                [], [], options)
            
    if exitflag == 1
        alphaDet{i} = x;
        chiDet(i) = fval;
    else
        allexitflagsEqualOne = false;
    end

end
chiDet = -chiDet; % minus before chi to contain a max instead of min 
 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
function [alphaCdet, uOptDet, exitflag] = control_Det (Up, Yp, Uf, Yf, GG, gg, umax, umin, chiDet) 
%Quadratic Problem for CONTROL
%alphaC stays for control
[alphaCdet,~, exitflag, ~] = quadprog(transpose(Uf)*Uf, [], vertcat(GG*Yf, Uf, -Uf), ... %%for -1< = Uf*alphaCdet <= 1
    vertcat(gg-chiDet, umax * ones(length(Uf(:,1)), 1), -umin * ones(length(Uf(:, 1)), 1)),...
    vertcat(Up, Yp), zeros(size(Up,1)+size(Yp,1),1));

%optimal input
uOptDet = Uf*alphaCdet;
 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBJECTIVE FCN
function [y,grady] = quadobj(x,Q,f,c)
y = 1/2*x'*Q*x + f'*x + c;
if nargout > 1
    grady = Q*x + f;
end
end

%NONLINEAR CONSTRAINT FCN
function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d)
jj = length(H); % jj is the number of inequality constraints
y = zeros(1,jj);
for i = 1:jj
    y(i) = 1/2*x'*H{i}*x + k{i}'*x + d{i};
end
yeq = [];
    
if nargout > 2
    grady = zeros(length(x),jj);
    for i = 1:jj
        grady(:,i) = H{i}*x + k{i};
    end
end
gradyeq = [];
end

%%same for control: but with equality instead of inequality
function [y,yeq,grady,gradyeq] = quadconstr_for_control(x,J,p) %1/2x'Jx + p'x = 0
jj = length(J); % jj is the number of inequality constraints
yeq = zeros(1,jj);
for i = 1:jj
    yeq(i) = 1/2*x'*J{i}*x + p{i}'*x ;
end
y = [];
    
if nargout > 2
    gradyeq = zeros(length(x),jj);
    for i = 1:jj
        gradyeq(:,i) = J{i}*x + p{i};
    end
end
grady = [];
end

%%%%%%%%%%%%%%%%%%%%% unfortunately, didn't work
%{
function chi_for_alpha = worst_for_pattern(alpha_pattern)
    global L
    global tau
    chi_noisy = 1:2*(L-tau);    
    for i = 1:length(chi_noisy)
        f_noisy = GG(i,:)*Yf + sign(GG(i,:))*transpose(GG(i,:))*eps*alpha_pattern;
        A_noisy = vertcat(Yp-ones(size(Yp,1),1)*alpha_pattern,-Yp-ones(size(Yp,1),1)*alpha_pattern, diag(-alpha_pattern));
        b_noisy = vertcat(y_init+eps.*ones(length(y_init),1),-y_init+eps.*ones(length(y_init),1), zeros(length(alpha_pattern),1));
        [~, chi_noisy(i)] = linprog(transpose(f_noisy), A_noisy, b_noisy); %, vertcat(Up, Uf), vertcat(u_init, zeros(L-length(u_init), 1)));
    end
    chi_for_alpha = chi_noisy;
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%functions for building some useful matrices
%% building Q
function Q = Q_alpha_ksi(rownum, colnum, horstart, horpattern)
%horpattern is a column vector, %horstart = start position for 11s, horlen -- num of 11s per row
      currentstart = horstart;
      Q = zeros(rownum, colnum);
      for i=1:rownum
          for j = currentstart:currentstart+length(horpattern)-1
              Q(i, j) = horpattern(j-currentstart+1);
          end
          currentstart = currentstart+1;
      end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%% specifying Hessian
function hess = quadhess(x,lambda,Q,H)
hess = Q;
jj = length(H); % jj is the number of inequality constraints
for i = 1:jj
    hess = hess + lambda.ineqnonlin(i)*H{i};
end
end



%%%%%%%%%%%%%%%%%%%%%%%%% Hessian when we have equalities instead
function hess = quadhess_for_control(x,lambda,Q,J)
hess = Q;
jj = length(J); % jj is the number of equality constraints
for i = 1:jj
    hess = hess + lambda.eqnonlin(i)*J{i};
end
end



function  [chiZeroStateNoisyUp, chiZeroStateNoisyDown] = chi_Zero_State_Noisy (Up, Yp, Uf, Yf, EPS, GG)
%v = (alpha, ksi)
lenalpha = length(Up(1,:));
Hu = vertcat(Up, Uf);
lenksi = length(Hu(:,1)) + length(Hu(1,:))-1;



%QUADRATIC constraints
%v’Jv + p’v = 0
quadconstrnum = length(Yp(:, 1));
JJ = cell(quadconstrnum, 1);

pp = cell(quadconstrnum, 1); 
p = horzcat(Yp, zeros(length(Yp(:,1)), lenksi));

for i = 1:quadconstrnum
   %J for quad constraint
   Q_constr = Q_alpha_ksi(lenalpha, lenksi, i, [1]);
   Jup = horzcat(zeros(lenalpha), Q_constr);
   Jdown = horzcat(transpose(Q_constr), zeros(lenksi));
   J = vertcat(Jup, Jdown);
   JJ{i} = J; 
   
   %p for quad constr
   pp{i} = transpose(p(i, :));
   %transpose merely for fmincon: argument needs to be a column vector
end



%%%%%%LINEAR constraints
%A, b for LINEAR ineq in ksi
A = vertcat(diag(ones(lenksi,1)), -diag(ones(lenksi,1)));
%populate for alpha missing
A = horzcat(zeros(length(A(:,1)),lenalpha), A);
b = vertcat(EPS*ones(lenksi,1), EPS*ones(lenksi,1));

%Aeq, beq for equality in alpha
%single Aeq
Aeq = vertcat(Up, Uf);
%populate for ksi missing
Aeq = horzcat(Aeq, zeros(length(Aeq(:, 1)), lenksi));

%multiple beq - inside main fcn




%%%%%%OBJECTIVE fcn
Q = cell(length(GG(:,1)), 1);
f = cell(length(GG(:,1)), 1);
for i=1:length(GG(:,1))
    % Q_i for i-th optimization problem
    Q_constr = Q_alpha_ksi(lenalpha, lenksi, length(Yp(:,1)) + 1, transpose(GG(i,:)));
    Qup = horzcat(zeros(lenalpha), Q_constr);
    Qdown = horzcat(transpose(Q_constr), zeros(lenksi));
    Q{i} = vertcat(Qup, Qdown);
    
    % f_i
    f{i} = transpose(horzcat(GG(i,:)*Yf, zeros(1, lenksi))); 
    %transpose merely for fmincon: argument needs to be a column vector
    
    % c_i = 0
end



%QUAD problem with QUAD constraints: we use FMINCON 
%%%%%% v = [alpha, ksi] is our decision variable
nonlconstr = @(v)quadconstr_for_control(v,JJ,pp);
v = cell(length(GG(:, 1)), 1);
uCount = length(Uf(:,1)); % one constraint per Uf row - for one [000 .. 010..0] vector
chiZeroStateNoisyUp = zeros(uCount, length(GG(:, 1)));%upper bound
chiZeroStateNoisyDown = zeros(uCount, length(GG(:, 1)));%lower bound
for un = 1:uCount
  for i=1:length(GG(:, 1))
      options = optimoptions(@fmincon,'Algorithm','interior-point', ...
      'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
      ... %'Display', 'iter', ...
      'HessianFcn',@(v,lambda)quadhess_for_control(v,lambda,-Q{i},JJ)); %-Q{i} for max instead of min

      fun = @(v)quadobj(v,-Q{i},-f{i},0); %-Q{i}, -f{i} for max instead of min
      v0 = vertcat(zeros(lenalpha, 1), zeros(lenksi, 1));
    
      %test for correctness: initial condition is feasible,
      %so our solution is "notlessthan" in that starting point
      %[chinoisynotlessthan(i), ~] = quadobj(v0, -Q{i}, -f{i}, 0); 
     
      u = zeros(uCount, 1);
      u(un) = 1; % all 0, single 1    
      beq = vertcat(zeros(length(Up(:, 1)), 1), u);
      
      %max
      [~,chiZeroStateNoisyUp(un,i),eflag,output,lambda] = ...
          fmincon(fun,v0,A,b, Aeq, beq,[],[],nonlconstr,options);
      chiZeroStateNoisyUp(un, i) = -chiZeroStateNoisyUp(un, i);
      
      %min
      fun = @(v)quadobj(v,Q{i},f{i},0); %min
      options = optimoptions(@fmincon,'Algorithm','interior-point', ...
      'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
      ... %'Display', 'iter', ...
      'HessianFcn',@(v,lambda)quadhess_for_control(v,lambda,Q{i},JJ)); %min
      [~,chiZeroStateNoisyDown(un, i),eflag,output,lambda] = ...
          fmincon(fun,v0,A,b, Aeq, beq,[],[],nonlconstr,options);
                                           %change back for min
  end
end
end



function [chiNoisy,v] = estimation_Noisy (Up, Yp, Uf, Yf, u_init, y_init, EPS, EPS_DATA, GG)
%v = (alpha, ksi)
lenalpha = length(Up(1,:));
Hu = vertcat(Up, Uf);
lenksi = length(Hu(:,1)) + length(Hu(1,:))-1;

%QUADRATIC constraints
%FOR one line, before vertcat
quadconstrnum = 2*length(Yp(:, 1));
HH = cell(quadconstrnum, 1);

kk = cell(quadconstrnum, 1);
%k for nonlinear constr 
kup = horzcat(Yp, zeros(length(Yp(:,1)), lenksi));
%where is one more variable k???

dd = cell(quadconstrnum, 1);
%d for nonlin constr
dup = -y_init-EPS*ones(length(y_init), 1);
ddown = y_init - EPS*ones(length(y_init),1);

for i = 1:quadconstrnum/2
   %H for quad constraint
   Q_constr = Q_alpha_ksi(lenalpha, lenksi, i, [1]);
   Hup = horzcat(zeros(lenalpha), Q_constr);
   Hdown = horzcat(transpose(Q_constr), zeros(lenksi));
   H = vertcat(Hup, Hdown);
   HH{i} = H; %upper bound
   HH{quadconstrnum/2 + i} = -H; %lower bound
   
   %k for quad constr
   kk{i} = transpose(kup(i, :)); %upper bound
   %transpose merely for fmincon: argument needs to be a column vector
   kk{quadconstrnum/2 + i} = -transpose(kup(i, :)); %lower bound
   
   %d for quad constr
   dd{i} = dup(i);%upper bound
   dd{quadconstrnum/2 + i} = ddown(i); %lower bound
end



%%%%%%LINEAR constraints
%A, b for LINEAR ineq in ksi
A = vertcat(diag(ones(lenksi,1)), -diag(ones(lenksi,1)));
%populate for alpha missing
A = horzcat(zeros(length(A(:,1)),lenalpha), A);
b = vertcat(EPS_DATA * ones(lenksi,1), EPS_DATA * ones(lenksi,1));

%Aeq, beq for equality in alpha
Aeq = vertcat(Up, Uf);
%populate for ksi missing
Aeq = horzcat(Aeq, zeros(length(Aeq(:, 1)), lenksi));
beq = vertcat(u_init,  zeros(length(Uf(:, 1)), 1));

%%%%%%OBJECTIVE fcn
Q = cell(length(GG(:,1)), 1);
f = cell(length(GG(:,1)), 1);
for i=1:length(GG(:,1))
    % Q_i for i-th optimization problem
    Q_constr = Q_alpha_ksi(lenalpha, lenksi, length(Yp(:,1)) + 1, transpose(GG(i,:)));
    Qup = horzcat(zeros(lenalpha), Q_constr);
    Qdown = horzcat(transpose(Q_constr), zeros(lenksi));
    Q{i} = vertcat(Qup, Qdown);
    
    % f_i
    f{i} = transpose(horzcat(GG(i,:)*Yf, zeros(1, lenksi))); 
    %transpose merely for fmincon: argument needs to be a column vector
    
    % c_i = 0
end

%QUAD problem with QUAD constraints: we use FMINCON 
%%%%%% v = [alpha_hat, ksi] is our decision variable
nonlconstr = @(v)quadconstr(v,HH,kk,dd);
v = cell(length(GG(:, 1)), 1);
chiNoisy = zeros(length(GG(:, 1)), 1);
%chinoisynotlessthan = zeros(length(GG(:, 1)), 1);
for i=1:length(GG(:, 1))
    options = optimoptions(@fmincon,'Algorithm','interior-point', ...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, ...
    ... %'Display', 'iter', ...
    'HessianFcn',@(v,lambda)quadhess(v,lambda,-Q{i},HH)); %-Q{i} for max instead of min

    fun = @(v)quadobj(v,-Q{i},-f{i},0); %-Q{i}, -f{i} for max instead of min
    
    %initial point corresponds to worst case for noise = 0
    %which we kept in alphaDeterminant
    %v0 = vertcat(alphaDeterminant{i}, zeros(lenksi, 1));
    v0 = rand(lenalpha + lenksi, 1); %different starting points
    %test for correctness: initial condition is feasible,
    %so our solution is "notlessthan" in that starting point
    %[chinoisynotlessthan(i), ~] = quadobj(v0, -Q{i}, -f{i}, 0); 
    disp iter
    [v{i},chiNoisy(i),eflag,output,lambda] = fmincon(fun,v0,A,b, Aeq, beq,[],[],nonlconstr,options);
end

chiNoisy = -chiNoisy;%minus before chinoisy to contain max instead of min
end



function [x_opt, solutionExists] = optC(Aup, Adown, b)
solutionExists = false;
dim = length(Aup(1, :)); %dimension of x_opt
b = vertcat(b, zeros(dim, 1));
Abase = zeros(length(Aup(:,1)), length(Aup(1,:)));

for pattern_num= 1:2.^dim
    x_pattern = de2bi(pattern_num-1, dim);
    x_pattern (x_pattern == 0) = -1; %change 0 to -1
    %x >= 0
    for i=1:dim
        if x_pattern == 1
            Abase(:,i) = Aup(:, i);
        end
        if x_pattern == -1
            Abase(:, i) = -Adown(:, i);
        end
    end
    
    A = vertcat(Abase, -eye(dim));
    %qp for each pattern
    [x, fval, exitflag, ~] = quadprog(eye(dim),[], A, b);
    if exitflag == 1
       soulutionExists = true;
    end
    if pattern_num == 1
        fval_opt = fval;
        x_opt = x;
    elseif fval < fval_opt
            fval_opt = fval;
            x_opt = x;
    end
        
end

end

