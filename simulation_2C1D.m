clear all;close all;clc

%% parameters
a1 = [0; 0.5];  a2 = [0; -0.5];    % cable starting points
b1 = [0; 0.05]; b2 = [0; -0.05];   % cable attachment points on the platform
m = 1;       % mass of the platform
g = [0; -9.81];   % gravitational accleration
dt = 0.1;
tf = 2;
t = 0:dt:tf;      % time
taumin=10;
taumax=100;

%% CVX parameters
Delta = 1e-6;  % slack variable
Lambda = 5;    % decay rate of an upper bound of V(x(t))
Gamma = 1;     % decay rate of a lower bound of B(x(t))
H=eye(2)*10;      % controls the weights of inputs
p=1;           % controls the degree of the slack variable

%% pre-allocating
x = zeros(2,length(t));
xdot = zeros(2,length(t));
xddot = zeros(2,length(t));
tau = zeros(2,length(t));

%% desired trajectory
x_i = [0; 0];     % initial position
x_f = [0; 0.2];   % final position
xdot_i = [0; 0];  % initial velocity
xdot_f = [0; 0];  % final velocity
k = t/tf;
x_d = x_i + (x_f-x_i) * (10*k.^3 - 15*k.^4 + 6*k.^5); % 5th order trajectory
vel_d = [zeros(2,1) diff(x_d,1,2)];
acc_d = [zeros(2,1) diff(vel_d,1,2)];

x(:,1) = x_i;
xdot(:,1) = xdot_i;
xddot = [zeros(2,1) diff(xdot,1,2)];

%% loop
for i=1:length(t)
    %% structure matrix
    l_vec1 = a1 - x(:,i) - b1; l_vec2 = a2 - x(:,i) - b2;
    u1 = l_vec1 / norm(l_vec1); u2 = l_vec2 / norm(l_vec2);
    J = [b1, b2];  
    
    %% SVD
%     [U_,S_,V_]=svd(J); s_=diag(S_);                   % compute SVD of J (see sec 2.1.1 of B2019)
%     for j=1:length(s_); if s_(j)>1e-8; r=j; end, end, % determine rank of J
%     U_underbar=U_(:,1:r); U_overbar=U_(:,r+1:2); S_underbar=S_(1:r,1:r); % partition SVD
%     V_underbar=V_(:,1:r); V_overbar=V_(:,r+1:2); % V_overbar is the null space of the matrix J
%     Jplus=V_underbar*inv(S_underbar)*U_underbar';  % determine pseudoinverse of J
% 
%     tauref = Jplus * m * ( acc_d(:,i) - g );
    
    %% tension optimization part..... this is where we need to write CLF-CBF-QP code
% CLF-CBF-QP
    cvx_begin quiet;
        variable tau_opt(2,1)
        tauref =  mean([taumin,taumax])*ones(2,1);
        err = x_d(:,i)-x(:,i);
%         CLF = [2*err.*xdot(:,i); 2*err.*g] + [zeros(2,1);2*J/m*err.*tau_opt] + Lambda*[err.^2; err.^2];
        CLF = 2*err.*xdot(:,i)+2*J/m*err.*tau_opt + Lambda*err.^2;
%         CBF = xddot(:,i).*xdot(:,i) + J/m.*xddot(:,i)*tau_opt + Gamma * (xdot(:,i)-0.1);
        minimize ( (tau_opt-tauref)'*H*(tau_opt-tauref)+p*Delta^2 );
        subject to;
            CLF <= Delta;
%             CBF >= 0;
            max(tau_opt) <= taumax;
            min(tau_opt) >= taumin;
    cvx_end;
        
    tau(:,i) = tau_opt;
    
%     tau(:,i) = Jplus * m * ( acc_d(:,i) - g ); % no optimization, no control
 
    %% solving dynamics, calculate position, velocity
    [x(:,i+1),xdot(:,i+1)] = rk4_ode([x(:,i); xdot(:,i)],m,J,tau(:,i),g,dt);
    xddot = [zeros(2,1) diff(xdot,1,2)];
    
end

figure; subplot(2,1,1); hold on; grid on; ylim([-0.6 0.6]); title('position');
plot(t,x_d(2,:),'b--'); plot(t,x(2,1:end-1),'b');
subplot(2,1,2); hold on; grid on; ylim([taumin-5 taumax+5]); title('tensions');
plot(t,tau(1,:),'r'); plot(t,tau(2,:),'b'); legend('upper cable','lower cable');

% figure; hold on; grid on; ylim([-0.6 0.6])
% plot([a1(1),a2(1)],[a1(2),a2(2)]); 

%% RK4 method
function [pos,vel] = rk4_ode(X, m, J, tau, g, dt)
    f1=ode(X, m, J, tau, g);
    f2=ode(X+dt/2 *f1, m, J, tau, g);
    f3=ode(X+dt/2 *f2, m, J, tau, g);
    f4=ode(X+dt * f3, m, J, tau, g);
    XX = (X + dt * (f1/6 + f2/3 +f3/3 + f4/6));
    pos = XX(1:2,1); vel = XX(3:4,1);
end

function x_prime = ode(X,m,J,tau,g)
   x_prime = [X(3);X(4); 1/m * (J*tau + m*g)];
end
