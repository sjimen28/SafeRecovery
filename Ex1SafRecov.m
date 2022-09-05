%--------------------------------------------------------------------------
% Matlab M-file Project: Safe Recovery from DoS Attacks @  Hybrid Systems Laboratory (HSL), 
% Filename: Ex1DSafRecov.m
%--------------------------------------------------------------------------
% Project: Example - 2D Linear TI
% Author: Santiago Jimenez Leudo
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Make sure to install HyEQ Toolbox (Beta) v3.0.0.22 from
%   https://www.mathworks.com/matlabcentral/fileexchange/102239-hybrid-equations-toolbox-beta 
%   (View Version History) 
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.4 Date: 09/05/2022 12:29:00

clear all
clc 
% --------------------------------------------------------------
%%% Initialization
% --------------------------------------------------------------
%   Paremeters: 
%   TSPAN, Ts, Ta, Tna
%   A, B, tildeC, barC, P, L, Q, tildeL
%   x0, hx0, barE
%   Modify any parameter in this section to simulate the system/case of interest

%% Simulation Horizon
TSPAN = [0 14];    % Second entry is the maximum amount of seconds allowed
Ts = 0.001;        % Steptime
t = 0:Ts:TSPAN(2);

% Attacks
Ta = 1.2;       % Maximum Length of a DoS attack
Tna = 1;      % Minimum lenght of intervals without attacks 
atck = 0;       % Attack flag (1 attack, 0 no attack)

%%% Continuous Dynamics
% \dot x=Ax+B*u     Flow Map
% y = Cx
A = [0 1 ; 0 0];
B = [0 ; 1];
tildeC = [1 0];     % Components that can't be attacked (Position)
barC = [0 1];       % Components that can be attacked (Velocity)
Cf = @(atck) [tildeC; (1-atck)*barC];
C = Cf(atck);

% Controlability
Co = ctrb(A,B);     % Controlability Matrix
if length(A) - rank(Co)>0
    disp('There are uncontrolable modes')
else
    disp('Controlable System')
end

% Observability
Obs = ctrb(A,C);    % Controlability Matrix
if length(A) - rank(Obs)>0
    disp('There are unobservable modes')
else
    disp('Observable System')
end

%%% Lyapunov Function
P = [1 0 ; 0 1];         % Needs to be symmetric and PD
if prod(eig(P))>0
    disp('P is positive definite')
else
    disp('P is NOT positive definite')
end
V = @(x) P*x^2;

%%% Estimation
% When No attacks:
% \dot \hat x=A \hatx+B*u + L(Cx-C\hat x) 
%
L = [2 -1 ;1 2];%[20 -1 ;1 15];             % Size 2x2
Q = -(A-L*C)'*P-P*(A-L*C);  % Needs to be symmetric and PD
if prod(eig(Q))>0
    disp('Q is positive definite')
else
    disp('Q is NOT positive definite')
end
eig(A-L*C)                  % L is designed so the eig(A-LC) are in the open left-half plane

% Under attacks:
% \dot \hat x=A \hatx+B*u + \tildeL(\tildeCx-\tildeC\hat x) 
%
tildeL = 4*[1.5 ; -1.4];          % Size 2x1
eig(A-tildeL*tildeC)        % Some of the eigenvalues are in the open left-half plane

%%% Initial State 
x0 = [-5; 2];       % State
hx0 = [-5.4; 2.3];      % Estimation
e0 = x0-hx0;

% Time Bounds
barE = norm(e0)+0.05;    % Initial error bound
theta = 0.1;%rand(1);
c1 = min(eig(Q));
c2 = min(eig(P))/max(eig(P));
C1 = sqrt(1/c2)*exp(-(1-theta)*c1*c2*Tna/(2*min(eig(P))))
delta = C1*barE;

[Phi,hA] = jordan(A-tildeL*tildeC);

Phi1 = Phi(:,1);        % Block associated to negative eigenvalues of matrix
Phi2 = Phi(:,2);        % Block associated to potentially nonnegative eigenvalues of matrix
%Phi =[Phi1, Phi2] ;    % Herminian Matrix s.t. Phi'*(A-tildeL*tildeC)*Phi=[A11 0; 0 A22]

hA11 = hA(1,1);
hA22 = hA(2,2);
hP = [1  0; 0 1];            % Needs to be positive definite
if prod(eig(hP))>0
    disp('\hatP is positive definite')
else
    disp('\hatP is NOT positive definite')
end
hQ = -hA11'*hP-hP*hA11;     % Needs to be symmetric and PD
if prod(eig(hQ))>0
    disp('\hatQ is positive definite')
else
    disp('\hatQ is NOT positive definite')
end

gc1 = min(eig(hQ));
gc2 = min(eig(hP))/max(eig(hP));

hc1 = norm(Phi)*norm(Phi1)*sqrt(1/gc2);
hc2 = norm(Phi)*norm(Phi2);

lambda1 = (1-theta)*gc1*gc2/(2*min(eig(hP)));
lambda2 = norm(hA22);

fun = @(t) -( hc1*exp(-lambda1*t)+hc2*lambda2*exp(lambda2*t));
[max_error_time_attacks, C2] = fminbnd(fun,0,Ta);
C2 = -C2;
C2 = -fun(Ta) %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C1*C2
figure(3)
hold on
plot(0:0.01:Ta,-fun(0:0.01:Ta))
%%
% Control Design
syms x1 x2
S = @(x1,x2) x1.^2 + 2*x2.^2+2*x1.*x2-25;   % Set to render safe
h = @(x1,x2) x1.^2 + 2*x2.^2+2*x1.*x2-17;   % Barrier function defining bar X_0
hf = x1.^2 + 2*x2.^2+2*x1.*x2-25;
g = gradient(hf,[x1,x2]);
nabla = inline(g');

%
% --------------------------------------------------------------
%%% System Evolution + Attack Launches
% --------------------------------------------------------------
%
% Define time of attacks

for j = 1:6                             % Number of attacks
    if j>1
        t1(j) = t2(j-1)+Tna;%+rand(1);    % Initial times of attacks
    else
        t1(j) = Tna;%+rand(1);            % Initial times of attacks
    end    
        t2(j) = t1(j)+Ta;%*rand(1);       % End times of attacks
end
    x(:,1) = x0;                % Initial State
    hx(:,1) = hx0;              % Initial Estimate
    e(:,1) = e0;                % Initial State
   
    % Controller
    K = [1 1];
    u(1) = -K*x(:,1);   % Intial Control Action arbitrarily DEFINED!!
    j = 1;              % Number of attacks

    for i = 1:length(t)-1

        x(:,i+1) = x(:,i)+(t(i+1)-t(i))*(A*x(:,i)+B*u(i)); % Evolve in state
        hx1t = hx(1,i); hx2t = hx(2,i);
        
        if j<=6 && (abs(t(i)-t1(j)) <= 0.002|| (atck(i) == 1 && t(i) == 0))        % Check if an attack just begun or the process start with an attack
            
            atck(i+1)=1;

        elseif j<=6 && (abs(t(i)-t2(j)) <= 0.002 || (atck(i) == 0 && t(i) == 0))    % Check if an attack just ended or the process starts with no attack

            atck(i+1) = 0;

            if atck(i) == 1

                j = j+1;                   % Report attack finished and start comparing with next attack time

            end

        else

            atck(i+1) = atck(i);           % No switch     

        end
        
        if atck(i) == 1         % Evolution during Attack

            hx(:,i+1) = hx(:,i)+(t(i+1)-t(i))*(A*hx(:,i)+B*u(i)+tildeL*(tildeC*x(:,i)-tildeC*hx(:,i))); % Evolve in estimation
            e(:,i+1) = e(:,i)+(t(i+1)-t(i))*((A-tildeL*tildeC)*e(:,i));                                 % Evolve in error
            ut = quadprog(eye(2),[0 0],[nabla(hx1t,hx2t)*B h(hx1t,hx2t)],-nabla(hx1t,hx2t)*(A*hx(:,i)+tildeL*(tildeC*x(:,i)-tildeC*hx(:,i))));

        else %if atck(i) == 0     % Evolution when there is no attack

            % C = Cf(atck(i));
            hx(:,i+1) = hx(:,i)+(t(i+1)-t(i))*(A*hx(:,i)+B*u(i)+L*(C*x(:,i)-C*hx(:,i)));    % Evolve in estimation 
            e(:,i+1) = e(:,i)+(t(i+1)-t(i))*((A-L*C)*e(:,i));                               % Evolve in error
            ut = quadprog(eye(2),[0 0],[nabla(hx1t,hx2t)*B h(hx1t,hx2t)],-nabla(hx1t,hx2t)*(A*hx(:,i)+L*(C*x(:,i)-C*hx(:,i))));

        end
        
        u(i+1) = ut(1);
        %u(i+1) = -K*x(:,i+1);  % Update input for next time step  
        
    end
    
%%
% --------------------------------------------------------------
%%% Plot
% --------------------------------------------------------------
%

figure (1)
clf
set(0,'defaultfigurecolor',[1 1 1])                 
set(0,'defaulttextinterpreter','latex')
set(gcf,'color','w');

subplot(5,1,1)
plot(t,x)                               % Plot Continuos Response
hold on
plot(t,hx)  
xlabel('$t$','Interpreter','Latex')
ylabel('$x$','Interpreter','Latex')
legend('$x_1$','$x_2$','$\hat x_1$','$\hat x_2$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

subplot(5,1,2)
plot(t, e)                              % Plot Error
ylabel('Error','Interpreter','Latex')
legend('$e_1$','$e_2$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

subplot(5,1,3)
norme = @(t)(e(1,t).^2+e(2,t).^2).^0.5;
plot(t, norme(':'))     % Plot Error
hold on
plot([t(1), t(end)], barE*[1, 1])
hold on
plot([t(1), t(end)], barE*C1*C2*[1, 1])
xlabel('$t$','Interpreter','Latex')
ylabel('norm(Error)','Interpreter','Latex')
legend('$|e|$','$\bar E$','$C_1 C_2 \bar E$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

subplot(5,1,4)
plot(t, atck)                           % Plot Error
xlabel('$t$','Interpreter','Latex')
ylabel('Attack','Interpreter','Latex')
ylim([-0.2, 1.2])
legend('$1 - Attack$', '0 - No attack','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

subplot(5,1,5)
plot(t, u)                              % Plot Action for ContSolution
xlabel('$t$','Interpreter','Latex')
ylabel('$u$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

%%
figure (2)
clf
set(0,'defaultfigurecolor',[1 1 1])                 
set(0,'defaulttextinterpreter','latex')
set(gcf,'color','w');

plot(x(1,:),x(2,:))                     % Plot Continuos Response
hold on
plot(hx(1,:),hx(2,:))
axis([-8,8,-6,6])
hold on 

H = fimplicit(S);
h1 = gcf;
axObjs = h1.Children;
dataObjs = axObjs.Children;
x1H = dataObjs(1).XData;
x2H = dataObjs(1).YData;
epsilon = 1.1;

if epsilon-barE*(1+C1*C2)>0
    disp('epsilon value is correct')
else
    disp('Modify epsilon to apply Theorem 1')
    epsilon
    barE
    res=sprintf('%.6f',barE*(1+C1*C2));
    disp(strcat('barE*(1+C1*C2)=', res))
end

%%% Defining X0 and \tilde X_0
clear Xz Xt
nnn = diff(hf,x2)/diff(hf,x1);
for k = 1:size(x1H,2)

n(k) = subs(nnn,[x1 x2],[x1H(1,k),x2H(1,k)]);  % Slope of normal vector at x1H(1,k),x2H(1,k)

    if k<266 && k>7

        Xz(:,k) = [x1H(k) ; x2H(k)]+epsilon*sqrt(1/(1+double(n(k))^2))*[1;n(k)];
        Xt(:,k) = Xz(:,k)-barE*sqrt(1/(1+double(n(k))^2))*[1;n(k)];

    else %if  k>=266

        Xz(:,k) = [x1H(k) ; x2H(k)]-epsilon*sqrt(1/(1+double(n(k))^2))*[1;n(k)];
        Xt(:,k) = Xz(:,k)+barE*sqrt(1/(1+double(n(k))^2))*[1;n(k)];

    end

end

plot(Xz(1,1:end),Xz(2,1:end))
hold on
plot(Xt(1,1:end),Xt(2,1:end))
hold on
%%% Definition of \hat X_0
xc = @(t) barE*cos(t)+x0(1);
yc = @(t) barE*sin(t)+x0(2);
fplot(xc,yc)
hold on
barX0 = fimplicit(h);       % barX0 set defined by barrier function

set(gca,'TickLabelInterpreter','latex')
xlabel('$x_1$','Interpreter','Latex')
ylabel('$x_2$','Interpreter','Latex')
legend('$x$','$\hat x$','$S$','$X_0$','$\tilde X$', '$\hat X_0(x_0)$','$\bar X$','Interpreter','Latex')

%% Debugging
% Find the indexes of the attacks in the time vector t

changedIndexes_t1 = diff(atck)>0;
changedIndexes_t2 = diff(atck)<0;
t1_indexes=find(changedIndexes_t1);  % Indexes at which attacks are launched
t2_indexes=find(changedIndexes_t2);  % Indexes at which attacks end

%norme(t1_indexes)
%norme(t2_indexes)

%disp('Error cummulated during each attack:')
%norme(t2_indexes)-norme(t1_indexes)  % Change of error over attacks

disp('Factor of increase of error during attack intervals:')
f1 = norme(t2_indexes)./norme(t1_indexes)

disp('Factor of increase of error during no-attack intervals:')
f2 = norme(t1_indexes)./[norme(1), norme(t2_indexes(1:5))]
% This shows that C_1 is not the actual increase over Tna !!

f1.*f2