%% CLF & PD control - 8 Cables 6DOF(x,y,z) model / update 21. 9.13.

clear all;close all;clc

graphs = 1;       % 0: graphs off,  1: graphs on
simulation = 1;   % 0: simulation off,  1: simulation on
method = 1;       % 1: CLF, 2: Null space(PD), 3: Normal

d = 6; % DOF
n = 8; % Cable number

%% time
dt = 0.1;
tf = 11;
t = 0:dt:tf;

%% Control
kp = -2.5;
kd = -1;

%% parameters
%% CVX parameters
Lambda = 10;    % decay rate of an upper bound of V(x(t))
H=10;             % controls the weights of inputs
p=0.001;           % controls the degree of the slack variable

%% outer structure
edge = 0.4318;     % outer structure edge size: 17in
eer = 0.05;        % moving platform edge size: 10cm
eeh = 0.005;       % moving platform height

a = [   0, edge, edge,    0,    0, edge, edge,    0;
        0,    0, edge, edge,    0,    0, edge, edge;
        0,    0,    0,    0, edge, edge, edge, edge];

%% cable attachment points on the platform
bin = 0.025;   % cross cable
b = [ -eer+bin,  eer-bin,  eer-bin, -eer+bin, -eer,  eer,  eer, -eer;
      -eer, -eer,  eer,  eer, -eer, -eer,  eer,  eer;
       eeh,  eeh,  eeh,  eeh, -eeh, -eeh, -eeh, -eeh];
bb= [ -eer,  eer,  eer, -eer, -eer,  eer,  eer, -eer;
      -eer, -eer,  eer,  eer, -eer, -eer,  eer,  eer;
       eeh,  eeh,  eeh,  eeh, -eeh, -eeh, -eeh, -eeh];  

   
m = 0.176;           % mass of the platform
g = -9.81;           % gravitational accleration (m/s^2)

Ig = zeros(3);       % moment of inertia
Ig(1,1) = (1/12) * m * (0.1^2+0.01^2); 
Ig(2,2) = (1/12) * m * (0.1^2+0.01^2); 
Ig(3,3) = (1/12) * m * (0.1^2+0.1^2); 
M = [m*eye(3), zeros(3); zeros(3), Ig*eye(3)];  % mass matrix

%% Tension bounds & cable stiffness: Kevlar
taumin=5;
taumax=100;
taumean=(taumin+taumax) / 2;
kappa=3000000;

%% pre-allocating
x = zeros(d,length(t));
xdot = zeros(d,length(t));
xddot = zeros(d,length(t));
x_d = zeros(d,length(t));
tau = zeros(n,length(t));

%% desired trajectory
x_i = [edge/2; edge/2; edge/2; 0; 0; 0];     % initial position
x_f = [edge/2; edge/2; edge/2; 0; 0; 0];     % final position
xdot_i = [0; 0; 0; 0; 0; 0];  % initial velocity
xdot_f = [0; 0; 0; 0; 0; 0];  % final velocity

%% Trajectory maker
k = t/tf;
x_d = x_i + (x_f-x_i) * (10*k.^3 - 15*k.^4 + 6*k.^5); % 5th order trajectory / linear

%% Circular trajectory
% theta = -2 * pi * (10*k.^3 - 15*k.^4 + 6*k.^5);
% r=0.07; x_d(1,:) = r*cos(theta)+edge/2; x_d(2,:) = r*sin(theta)+edge/2;

% writematrix(x_d,'check.xls');

% %% Reading Trajectory
% x_d = round(readmatrix('motion_CD86_fastroll.xls'),4);    % dt=0.1, 10 seconds fastroll
% x_d = round(readmatrix('motion_CD86_round.xls'),4);    % dt=0.02, 1 seconds circle
% x_d = round(readmatrix('motion_CD86_spiral.xls'),4);    % dt=0.05, 4 seconds spiral
% x_d = round(readmatrix('motion_CD86_fastyaw.xls'),4); x_d(6,:) = x_d(6,:)/2;   % dt=0.05, 2 seconds yaw

x_d = round(readmatrix('motion_CD86_v2.xls'),4);  % dt=0.1, 11 seconds periodic motions
vel_d = [zeros(d,1) diff(x_d,1,2)]/dt;
acc_d = [zeros(d,1) diff(vel_d,1,2)]/dt;

x(:,1) = x_i;
xdot(:,1) = xdot_i;
xddot = [zeros(d,1) diff(xdot,1,2)]/dt;

%% loop
for i=1:length(t)
    the=x(4,i); phi=x(5,i); psi=x(6,i);
    %% Quaternion
    qx = [cos(the/2); sin(the/2); 0; 0]; qy = [cos(phi/2); 0; sin(phi/2); 0]; qz = [cos(psi/2); 0; 0; sin(psi/2)];
    %% Quaternion multiplication
    q1 = quat_multiply(qz,qy); q = quat_multiply(q1,qx); q = q/norm(q); % unit quaternion

    %% Jacobian matrix
    for j = 1:n
        bR(:,j,i) = quat_rotate(q,b(:,j));  % quat rotation
        bbR(:,j,i) = quat_rotate(q,bb(:,j));  % quat rotation
        l_vec(:,j) = a(:,j) - x(1:3,i) - bR(:,j,i);
        u(:,j) = l_vec(:,j)/norm(l_vec(:,j));
        crossbRu(:,j) = cross(bR(:,j,i),u(:,j));
    end
    
    J = [u;crossbRu];
    J = J;
    
    %% SVD
    [U_,S_,V_]=svd(J); s_=diag(S_);                   % compute SVD of J (see sec 2.1.1 of B2019)
    for j=1:length(s_); if s_(j)>1e-8; r=j; end, end, % determine rank of J
    U_underbar=U_(:,1:r); U_overbar=U_(:,r+1:end); S_underbar=S_(1:r,1:r); % partition SVD
    V_underbar=V_(:,1:r); V_overbar=V_(:,r+1:end); % V_overbar is the null space of the matrix J
    Jplus=V_underbar*inv(S_underbar)*U_underbar';  % determine pseudoinverse of J
    
    %% Apply mesuarment errors
%     x(:,i)=x(:,i)+wgn(6,1,1)*0.0001;
    
    %% Errors
    errp = x(:,i)-x_d(:,i);
    errd = xdot(:,i)-vel_d(:,i);
    
    if method==1
    %% method 1: CLF control
%         cvx_begin quiet;
%             variables tau_opt(8,1) Delta(1,1)
%             Vx = (x(:,i)-x_d(:,i)).^2;
%             LfV = 0;
%             LgV = 2*(x(:,i)-x_d(:,i));
%             uinput = inv(M) * (J*tau_opt - M*[0;0;g;0;0;0]);
%             CLF = LfV + LgV.*uinput + Lambda*Vx;
%             minimize ( (tau_opt-taumean)'*H*(tau_opt-taumean) + p*Delta^2 );
%             subject to
%                 CLF <= Delta;
%                 taumin <= tau_opt <= taumax
%                 abs(J*tau_opt-M*(acc_d(:,i)+[0;0;g;0;0;0])) <= 0;            
%         cvx_end;
        cvx_begin quiet;
            variables tau_opt(8,1) Delta(1,1)
            Vx = [(x(:,i)-x_d(:,i)).^2; xdot(:,i)-vel_d(:,i).^2];
            LfV = 0;
            LgV = 2*[(x(:,i)-x_d(:,i)); (xdot(:,i)-vel_d(:,i))];
            uinput = [zeros(6,1); inv(M) * (J*tau_opt - M*[0;0;g;0;0;0])];
            CLF = LfV + LgV.*uinput + Lambda*Vx;
            minimize ( 1/2*(tau_opt-taumean)'*H*(tau_opt-taumean) + p*Delta^2 );
            subject to
                CLF <= Delta;
                taumin <= tau_opt <= taumax
                abs(J*tau_opt-M*(acc_d(:,i)+[0;0;g;0;0;0])) <= 0;            
        cvx_end;

        tau(:,i) = tau_opt;
        
        uinput_rec(:,i)=uinput;
        
        Vx_rec(:,i) = Vx;
    end
%     
    if method==2
    %% method 2: PD control tension optimization w/ null space optimization(standard deviation)
        cvx_begin quiet;
            variable h(n,1);
            taumean =  mean([taumin,taumax]);
            tau_opt = Jplus * M * (acc_d(:,i) + kp*errp + kd*errd + [0;0;g;0;0;0]) + ( eye(n) - (Jplus*J) ) * h;
            minimize ( sum( (tau_opt-taumean).^2 ) / 8)   % square of standard deviation
            subject to;
                taumin <= tau_opt <= taumax 
                abs(J*tau_opt-M*(acc_d(:,i)+[0;0;g;0;0;0]+ kp*errp + kd*errd)) <= 1e-8;
        cvx_end;
        h(isnan(h))=0;
        tau(:,i) = Jplus * M * (acc_d(:,i) + kp*errp + kd*errd + [0;0;g;0;0;0]) + ( eye(n) - (Jplus*J) ) * h;
        
        uinput_rec(:,i)=[zeros(6,1); inv(M) * (J*tau_opt - M*[0;0;g;0;0;0])];
    end

    if method==3
    %% method 3: no optimization, solving non linear equation
        tau(:,i) = Jplus * M * (acc_d(:,i) + [0;0;g;0;0;0] + kp*errp + kd*errd);
    end
    
    %% Cable length
    U(:,i) = kappa*[norm(l_vec(:,1)); norm(l_vec(:,2)); norm(l_vec(:,3)); norm(l_vec(:,4)); norm(l_vec(:,5)); norm(l_vec(:,6)); norm(l_vec(:,7)); norm(l_vec(:,8))]./(tau(:,i)+kappa);
    
    %% solving dynamics, calculate position, velocity
    [x(:,i+1),xdot(:,i+1)] = rk4_ode([x(:,i); xdot(:,i)],M,J,tau(:,i),g,dt);
    xddot = [zeros(d,1) diff(xdot,1,2)];
    
    err_rec(:,i) = errp;
end

if method == 1, writematrix(tau,'tau_CLF.xls'); writematrix(err_rec, 'err_CLF.xls');, end
if method == 2, writematrix(tau,'tau_nul.xls'); writematrix(err_rec, 'err_nul.xls');, end
if method == 3, writematrix(tau,'tau_nor.xls'); writematrix(err_rec, 'err_nor.xls');, end

if graphs == 1,
%% Graphs
figure; subplot(2,2,1); hold on; grid on;  ylim([-edge*0.1 edge*1.1]); title('Position'); xlabel('Time(s)'); ylabel('Position(m)');
plot(t,x_d(1,:),'r--'); plot(t,x(1,1:end-1),'r'); plot(t,x_d(2,:),'b--'); plot(t,x(2,1:end-1),'b');  plot(t,x_d(3,:),'g--'); plot(t,x(3,1:end-1),'g');
subplot(2,2,2); hold on; grid on; ylim([-90 90]); title('Rotation'); xlabel('Time(s)'); ylabel('Angle(degree)');
plot(t,x_d(4,:)*180/pi,'r--'); plot(t,x(4,1:end-1)*180/pi,'r'); plot(t,x_d(5,:)*180/pi,'b--'); plot(t,x(5,1:end-1)*180/pi,'b');  plot(t,x_d(6,:)*180/pi,'g--'); plot(t,x(6,1:end-1)*180/pi,'g');

subplot(2,2,3); hold on; grid on; ylim([taumin-5 taumax+5]); title('Cable tensions'); xlabel('Time(s)'); ylabel('Tension(N)');
plot(t,tau(1:4,:),'r'); plot(t,tau(5:8,:),'b');
plot(t,taumin*ones(1,size(t,2)),'k--');plot(t,taumax*ones(1,size(t,2)),'k--')
subplot(2,2,4); hold on; grid on; title('Cable lengths'); xlabel('Time(s)'); ylabel('Length(m)');
plot(t,U(1:4,:),'r'); plot(t,U(5:8,:),'b');
figure;grid on; plot(t,err_rec);title('Position Errors');

end

%% simulation
if simulation == 1,
    figure; hold on; grid on; view(-170,20); xlabel('x'); ylabel('y'); zlabel('z'); xlim([-edge*0.1 edge*1.1]); ylim([-edge*0.1 edge*1.1]); zlim([-edge*0.1 edge*1.1])
    pause_toggle=0;
    for i=1:n
        l(i) = plot3([a(1,i),x(1,1)+bR(1,i)],[a(2,i),x(2,1)+bR(2,i)],[a(3,i),x(3,1)+bR(3,i,1)],'k','LineWidth',1.2);
    end
    
    %% platform
    % upper shape
    platform(1) = plot3([x(1,1)+bbR(1,1),x(1,1)+bbR(1,2)],[x(2,1)+bbR(2,1),x(2,1)+bbR(2,2)],[x(3,1)+bbR(3,1),x(3,1)+bbR(3,2)],'b','LineWidth',1.2);
    platform(2) = plot3([x(1,1)+bbR(1,3),x(1,1)+bbR(1,2)],[x(2,1)+bbR(2,3),x(2,1)+bbR(2,2)],[x(3,1)+bbR(3,3),x(3,1)+bbR(3,2)],'b','LineWidth',1.2);
    platform(3) = plot3([x(1,1)+bbR(1,3),x(1,1)+bbR(1,4)],[x(2,1)+bbR(2,3),x(2,1)+bbR(2,4)],[x(3,1)+bbR(3,3),x(3,1)+bbR(3,4)],'b','LineWidth',1.2);
    platform(4) = plot3([x(1,1)+bbR(1,1),x(1,1)+bbR(1,4)],[x(2,1)+bbR(2,1),x(2,1)+bbR(2,4)],[x(3,1)+bbR(3,1),x(3,1)+bbR(3,4)],'b','LineWidth',1.2);
    % side shape
    platform(5) = plot3([x(1,1)+bbR(1,1),x(1,1)+bbR(1,5)],[x(2,1)+bbR(2,1),x(2,1)+bbR(2,5)],[x(3,1)+bbR(3,1),x(3,1)+bbR(3,5)],'b','LineWidth',1.2);
    platform(6) = plot3([x(1,1)+bbR(1,6),x(1,1)+bbR(1,2)],[x(2,1)+bbR(2,6),x(2,1)+bbR(2,2)],[x(3,1)+bbR(3,6),x(3,1)+bbR(3,2)],'b','LineWidth',1.2);
    platform(7) = plot3([x(1,1)+bbR(1,3),x(1,1)+bbR(1,7)],[x(2,1)+bbR(2,3),x(2,1)+bbR(2,7)],[x(3,1)+bbR(3,3),x(3,1)+bbR(3,7)],'b','LineWidth',1.2);
    platform(8) = plot3([x(1,1)+bbR(1,8),x(1,1)+bbR(1,4)],[x(2,1)+bbR(2,8),x(2,1)+bbR(2,4)],[x(3,1)+bbR(3,8),x(3,1)+bbR(3,4)],'b','LineWidth',1.2);
    % bottom shape
    platform(9) = plot3([x(1,1)+bbR(1,6),x(1,1)+bbR(1,5)],[x(2,1)+bbR(2,6),x(2,1)+bbR(2,5)],[x(3,1)+bbR(3,6),x(3,1)+bbR(3,5)],'b','LineWidth',1.2);
    platform(10)= plot3([x(1,1)+bbR(1,6),x(1,1)+bbR(1,7)],[x(2,1)+bbR(2,6),x(2,1)+bbR(2,7)],[x(3,1)+bbR(3,6),x(3,1)+bbR(3,7)],'b','LineWidth',1.2);
    platform(11)= plot3([x(1,1)+bbR(1,8),x(1,1)+bbR(1,7)],[x(2,1)+bbR(2,8),x(2,1)+bbR(2,7)],[x(3,1)+bbR(3,8),x(3,1)+bbR(3,7)],'b','LineWidth',1.2);
    platform(12)= plot3([x(1,1)+bbR(1,8),x(1,1)+bbR(1,5)],[x(2,1)+bbR(2,8),x(2,1)+bbR(2,5)],[x(3,1)+bbR(3,8),x(3,1)+bbR(3,5)],'b','LineWidth',1.2);

    %% outer structure
    plot3([a(1,1),a(1,2)],[a(2,1),a(2,2)],[a(3,1),a(3,2)],'g'); plot3([a(1,3),a(1,2)],[a(2,3),a(2,2)],[a(3,3),a(3,2)],'g')
    plot3([a(1,1),a(1,4)],[a(2,1),a(2,4)],[a(3,1),a(3,4)],'g'); plot3([a(1,3),a(1,4)],[a(2,3),a(2,4)],[a(3,3),a(3,4)],'g')
    plot3([a(1,5),a(1,6)],[a(2,5),a(2,6)],[a(3,5),a(3,6)],'g'); plot3([a(1,7),a(1,6)],[a(2,7),a(2,6)],[a(3,7),a(3,6)],'g')
    plot3([a(1,5),a(1,8)],[a(2,5),a(2,8)],[a(3,5),a(3,8)],'g'); plot3([a(1,7),a(1,8)],[a(2,7),a(2,8)],[a(3,7),a(3,8)],'g')
    plot3([a(1,1),a(1,5)],[a(2,1),a(2,5)],[a(3,1),a(3,5)],'g'); plot3([a(1,3),a(1,7)],[a(2,3),a(2,7)],[a(3,3),a(3,7)],'g')
    plot3([a(1,2),a(1,6)],[a(2,2),a(2,6)],[a(3,2),a(3,6)],'g'); plot3([a(1,8),a(1,4)],[a(2,8),a(2,4)],[a(3,8),a(3,4)],'g')
    p = plot3(x(1,1), x(2,1), x(3,1),'r*');
    pause
    delete(l); delete(p); delete(platform); 
    
    %% Plotting motion
    for i=1:size(t,2)
        for j=1:n        
            if tau(j,i) <= taumin || tau(j,i) >= taumax,
                l(j) = plot3([a(1,j),x(1,i)+bR(1,j,i)],[a(2,j),x(2,i)+bR(2,j,i)],[a(3,j),x(3,i)+bR(3,j,i)],'r','LineWidth',1.5);
                pause_toggle = 1;
            end
            l(j) = plot3([a(1,j),x(1,i)+bR(1,j,i)],[a(2,j),x(2,i)+bR(2,j,i)],[a(3,j),x(3,i)+bR(3,j,i)],'k','LineWidth',1.2);
        end
        %% platform
        % upper shape
        platform(1) = plot3([x(1,i)+bbR(1,1,i),x(1,i)+bbR(1,2,i)],[x(2,i)+bbR(2,1,i),x(2,i)+bbR(2,2,i)],[x(3,i)+bbR(3,1,i),x(3,i)+bbR(3,2,i)],'b','LineWidth',1.2);
        platform(2) = plot3([x(1,i)+bbR(1,3,i),x(1,i)+bbR(1,2,i)],[x(2,i)+bbR(2,3,i),x(2,i)+bbR(2,2,i)],[x(3,i)+bbR(3,3,i),x(3,i)+bbR(3,2,i)],'b','LineWidth',1.2);
        platform(3) = plot3([x(1,i)+bbR(1,3,i),x(1,i)+bbR(1,4,i)],[x(2,i)+bbR(2,3,i),x(2,i)+bbR(2,4,i)],[x(3,i)+bbR(3,3,i),x(3,i)+bbR(3,4,i)],'b','LineWidth',1.2);
        platform(4) = plot3([x(1,i)+bbR(1,1,i),x(1,i)+bbR(1,4,i)],[x(2,i)+bbR(2,1,i),x(2,i)+bbR(2,4,i)],[x(3,i)+bbR(3,1,i),x(3,i)+bbR(3,4,i)],'b','LineWidth',1.2);
        % side shape
        platform(5) = plot3([x(1,i)+bbR(1,1,i),x(1,i)+bbR(1,5,i)],[x(2,i)+bbR(2,1,i),x(2,i)+bbR(2,5,i)],[x(3,i)+bbR(3,1,i),x(3,i)+bbR(3,5,i)],'b','LineWidth',1.2);
        platform(6) = plot3([x(1,i)+bbR(1,6,i),x(1,i)+bbR(1,2,i)],[x(2,i)+bbR(2,6,i),x(2,i)+bbR(2,2,i)],[x(3,i)+bbR(3,6,i),x(3,i)+bbR(3,2,i)],'b','LineWidth',1.2);
        platform(7) = plot3([x(1,i)+bbR(1,3,i),x(1,i)+bbR(1,7,i)],[x(2,i)+bbR(2,3,i),x(2,i)+bbR(2,7,i)],[x(3,i)+bbR(3,3,i),x(3,i)+bbR(3,7,i)],'b','LineWidth',1.2);
        platform(8) = plot3([x(1,i)+bbR(1,8,i),x(1,i)+bbR(1,4,i)],[x(2,i)+bbR(2,8,i),x(2,i)+bbR(2,4,i)],[x(3,i)+bbR(3,8,i),x(3,i)+bbR(3,4,i)],'b','LineWidth',1.2);
        % bottom shape
        platform(9) = plot3([x(1,i)+bbR(1,6,i),x(1,i)+bbR(1,5,i)],[x(2,i)+bbR(2,6,i),x(2,i)+bbR(2,5,i)],[x(3,i)+bbR(3,6,i),x(3,i)+bbR(3,5,i)],'b','LineWidth',1.2);
        platform(10)= plot3([x(1,i)+bbR(1,6,i),x(1,i)+bbR(1,7,i)],[x(2,i)+bbR(2,6,i),x(2,i)+bbR(2,7,i)],[x(3,i)+bbR(3,6,i),x(3,i)+bbR(3,7,i)],'b','LineWidth',1.2);
        platform(11)= plot3([x(1,i)+bbR(1,8,i),x(1,i)+bbR(1,7,i)],[x(2,i)+bbR(2,8,i),x(2,i)+bbR(2,7,i)],[x(3,i)+bbR(3,8,i),x(3,i)+bbR(3,7,i)],'b','LineWidth',1.2);
        platform(12)= plot3([x(1,i)+bbR(1,8,i),x(1,i)+bbR(1,5,i)],[x(2,i)+bbR(2,8,i),x(2,i)+bbR(2,5,i)],[x(3,i)+bbR(3,8,i),x(3,i)+bbR(3,5,i)],'b','LineWidth',1.2);

        p = plot3(x(1,i), x(2,i), x(3,i),'r*');
        drawnow;
        delete(l);
    if pause_toggle==1
        pause;
        pause_toggle==0;
    end
     delete(p); delete(platform); 
    end
            
    for j=1:n
        l(j) = plot3([a(1,j),x(1,i)+bR(1,j)],[a(2,j),x(2,i)+bR(2,j)],[a(3,j),x(3,i)+bR(3,j)],'k','LineWidth',1.2);
    end
    %% platform
    % upper shape
    platform(1) = plot3([x(1,i)+bbR(1,1,i),x(1,i)+bbR(1,2,i)],[x(2,i)+bbR(2,1,i),x(2,i)+bbR(2,2,i)],[x(3,i)+bbR(3,1,i),x(3,i)+bbR(3,2,i)],'b','LineWidth',1.2);
    platform(2) = plot3([x(1,i)+bbR(1,3,i),x(1,i)+bbR(1,2,i)],[x(2,i)+bbR(2,3,i),x(2,i)+bbR(2,2,i)],[x(3,i)+bbR(3,3,i),x(3,i)+bbR(3,2,i)],'b','LineWidth',1.2);
    platform(3) = plot3([x(1,i)+bbR(1,3,i),x(1,i)+bbR(1,4,i)],[x(2,i)+bbR(2,3,i),x(2,i)+bbR(2,4,i)],[x(3,i)+bbR(3,3,i),x(3,i)+bbR(3,4,i)],'b','LineWidth',1.2);
    platform(4) = plot3([x(1,i)+bbR(1,1,i),x(1,i)+bbR(1,4,i)],[x(2,i)+bbR(2,1,i),x(2,i)+bbR(2,4,i)],[x(3,i)+bbR(3,1,i),x(3,i)+bbR(3,4,i)],'b','LineWidth',1.2);
    % side shape
    platform(5) = plot3([x(1,i)+bbR(1,1,i),x(1,i)+bbR(1,5,i)],[x(2,i)+bbR(2,1,i),x(2,i)+bbR(2,5,i)],[x(3,i)+bbR(3,1,i),x(3,i)+bbR(3,5,i)],'b','LineWidth',1.2);
    platform(6) = plot3([x(1,i)+bbR(1,6,i),x(1,i)+bbR(1,2,i)],[x(2,i)+bbR(2,6,i),x(2,i)+bbR(2,2,i)],[x(3,i)+bbR(3,6,i),x(3,i)+bbR(3,2,i)],'b','LineWidth',1.2);
    platform(7) = plot3([x(1,i)+bbR(1,3,i),x(1,i)+bbR(1,7,i)],[x(2,i)+bbR(2,3,i),x(2,i)+bbR(2,7,i)],[x(3,i)+bbR(3,3,i),x(3,i)+bbR(3,7,i)],'b','LineWidth',1.2);
    platform(8) = plot3([x(1,i)+bbR(1,8,i),x(1,i)+bbR(1,4,i)],[x(2,i)+bbR(2,8,i),x(2,i)+bbR(2,4,i)],[x(3,i)+bbR(3,8,i),x(3,i)+bbR(3,4,i)],'b','LineWidth',1.2);
    % bottom shape
    platform(9) = plot3([x(1,i)+bbR(1,6,i),x(1,i)+bbR(1,5,i)],[x(2,i)+bbR(2,6,i),x(2,i)+bbR(2,5,i)],[x(3,i)+bbR(3,6,i),x(3,i)+bbR(3,5,i)],'b','LineWidth',1.2);
    platform(10)= plot3([x(1,i)+bbR(1,6,i),x(1,i)+bbR(1,7,i)],[x(2,i)+bbR(2,6,i),x(2,i)+bbR(2,7,i)],[x(3,i)+bbR(3,6,i),x(3,i)+bbR(3,7,i)],'b','LineWidth',1.2);
    platform(11)= plot3([x(1,i)+bbR(1,8,i),x(1,i)+bbR(1,7,i)],[x(2,i)+bbR(2,8,i),x(2,i)+bbR(2,7,i)],[x(3,i)+bbR(3,8,i),x(3,i)+bbR(3,7,i)],'b','LineWidth',1.2);
    platform(12)= plot3([x(1,i)+bbR(1,8,i),x(1,i)+bbR(1,5,i)],[x(2,i)+bbR(2,8,i),x(2,i)+bbR(2,5,i)],[x(3,i)+bbR(3,8,i),x(3,i)+bbR(3,5,i)],'b','LineWidth',1.2);
    p = plot3(x(1,i), x(2,i), x(3,i),'r*');
end

%% RK4 method
function [pos,vel] = rk4_ode(X, M, J, tau, g, dt)
    f1=ode(X, M, J, tau, g);
    f2=ode(X+dt/2 *f1, M, J, tau, g);
    f3=ode(X+dt/2 *f2, M, J, tau, g);
    f4=ode(X+dt * f3, M, J, tau, g);
    XX = (X + dt * (f1/6 + f2/3 +f3/3 + f4/6));
    pos = XX(1:6,1); vel = XX(7:12,1);
end

function x_prime = ode(X,M,J,tau,g)
   x_prime = [X(7:12); inv(M) * (J*tau - M*[0;0;g;0;0;0])];
end


function [s]=quat_multiply(q,r)
% Quaternion multiplication s=q*r  (use [0;r] on input if r has no real part)
q0=q(1); q1=q(2); q2=q(3); q3=q(4); % This line helpful because Matlab enumerates from 1
r0=r(1); r1=r(2); r2=r(3); r3=r(4);
s=[r0*q0-r1*q1-r2*q2-r3*q3;
   r0*q1+r1*q0-r2*q3+r3*q2;
   r0*q2+r1*q3+r2*q0-r3*q1;
   r0*q3-r1*q2+r2*q1+r3*q0]; 
end % function quat_multiply

function [v]=quat_rotate(q,u)
% Note: q needs to be a unit quaternion for this function to work correctly!
q0=q(1); q1=q(2); q2=q(3); q3=q(4); % This line helpful because Matlab enumerates from 1
v=[(q0^2+q1^2-q2^2-q3^2)   2*(q1*q2-q0*q3)     2*(q1*q3+q0*q2);
    2*(q1*q2+q0*q3)    (q0^2-q1^2+q2^2-q3^2)   2*(q2*q3-q0*q1);
    2*(q1*q3-q0*q2)        2*(q2*q3+q0*q1)  (q0^2-q1^2-q2^2+q3^2)]*u;
end % function quat_rotate
