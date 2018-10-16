clc
clear all
% close all
home

global R2D
global D2R
global l1
global l2
global dt
global tf
global ay
global by
global az
global w
global c
global h
global tau
global error_coupling
global n_dofs
global n_bfs
global goal
global y_0


% Variable for converting RADIAN -> DEGREE or DEGREE -> RADIAN %
R2D = 180/pi;
D2R = pi/180;

%-------------------------------------------%
% 0. Record Data
%-------------------------------------------%

q_d = load("data/walking_joint_angle.txt");
dt = 0.001;
timestep = length(q_d);
tf = dt*(timestep-1);
T = 0:dt:tf;



%-------------------------------------------%
% 1. Set the Parameters
%-------------------------------------------%

% Robot Parameters 
l1 = 0.2;
l2 = 0.2;

% Constant Parameters %
az = 1;
ay = [7, 7];
by = [7, 7];

% Canonical System %
tau = 1;
error = 0.0;
error_coupling = 1.0 / (1.0+error);

% DMP Parameters %
n_dofs = 2; % number of degree of freedom
n_bfs = 100; % number of basis function per DMP



%-------------------------------------------%
% 2. Initialize Canonical System
%   Output : z[timestep]
%-------------------------------------------%

% Canonical System %
z = exp(-az/tau*T);

% Basis Function %
des_c = linspace(0, tf, n_bfs);     
c = exp(-az/tau*des_c);             % mean
h = n_bfs^1.5./c;                   % variance
psi = zeros(timestep, n_bfs);

% Psi %
for b=1:n_bfs
    psi(:,b) = exp(-h(b)*(z-c(b)).^2);
end




%------------------------------------------------------%
% 3. Get desired path from kinesthetic teach-in
%   Output : path[path_len, n_dofs]
%          : y_0, goal[n_dofs]
%------------------------------------------------------%

% Cartesian Space(m) desired path
x_d = FK(q_d);
path = x_d;
path_dt = 0.005;
path_len = length(path);
path_T = 0:path_dt:(path_len-1)*path_dt;



%--------------------------------------------------------%
% 4. Imitate the desired path using Interpolation
%   Output : y_d, dy_d, ddy_d[timestep, n_dofs]
%--------------------------------------------------------%

% Interpolate the desired path
y_T = linspace(0, tf, timestep);

y_d = zeros(timestep, n_dofs);
dy_d = zeros(timestep, n_dofs);
ddy_d = zeros(timestep, n_dofs);

for n=1:n_dofs
    y_d(:,n) = interp1(T, path(:,n), y_T);
end

% Velocity and Accelearation of y_d
for k=1:timestep-1
    dy_d(k+1, :) = (y_d(k+1, :)-y_d(k, :))/dt;
    ddy_d(k+1,:) = (dy_d(k+1, :)-dy_d(k, :))/dt;
end

% Set initial state and goal 
y_0 = y_d(1, :);
goal = y_d(path_len, :);




%-----------------------------------------------------------------------%
% 5. Find the reference force required to move along desired path
%   Output : f_ref
%-----------------------------------------------------------------------%
% Calculate reference force value from desired path
f_ref = zeros(timestep, n_dofs);
for n=1:n_dofs
    f_ref(:, n) = ddy_d(:, n) - ay(n)*(by(n)*(goal(n)-y_d(:,n)) - dy_d(:,n));
end


%-----------------------------------------------------------------------%
% 6. Generate weights to realize reference force
%   Output : weights
%-----------------------------------------------------------------------%
% Weight %
% Compute weights via Locally Weighted Regression
w = zeros(n_dofs, n_bfs);
for n=1:n_dofs
    % spatial scailing term
    k = goal(n)-y_0(n);
    for b=1:n_bfs
        numer = sum((z') .* psi(:,b) .* (f_ref(:,n))); % [timestep]->[1]
        denom = sum((z').^2 .* psi(:, b));
        if (k*denom == 0)
            w(n, b) = 0;
        else
            w(n, b) = numer / (k*denom);
        end        
    end
end

% dlmwrite('w_left_back.txt',w)
% dlmwrite('c_left_back.txt',c)
% dlmwrite('h_left_back.txt',h)

%-----------------------------------------------------------------------%
% 7. Change the goal
%   Output : 
%-----------------------------------------------------------------------%
%goal = [-0.1 ,0.2]; %left hit(test)

%-----------------------------------------------------------------------%
% 8. Generate trajectory to track
%   Output : y_t, dy_t, ddy_t -> Joint space
%-----------------------------------------------------------------------%
% Initialize tracking trajectory
y_t = y_0;  %initial
dy_t = zeros(n_dofs);
ddy_t = zeros(n_dofs);
z_t = 1.0;
f_t = [0, 0];

% Iteration numbers
n = 1;        %Iterator for main loop
n_trj = 1;

% Plot Setting %
figure(1)
title('Animation')
grid
hold on
axis([-0.5 0.5 -0.5 0.5]);
   Ax1 = [0, l1];
   Ay1 = [0, 0];
   Ax2 = [l1, l1+l2];
   Ay2 = [0, 0];
   p1 = line(Ax1,Ay1,'LineWidth',[5],'Color','b');
   p2 = line(Ax2,Ay2,'LineWidth',[5],'Color','c');


% Robot Implementation 
for i = 0 : dt : tf
    
    % Run Canonical System %
    z_t = z_t + (-az*z_t*error_coupling)*tau*dt;  % [1]
    [y_t, dy_t, ddy_t, f_t] = dmp_step(z_t, y_t, dy_t, ddy_t);
    
    % Inverse Kinematics
    x2 = y_t(1);    y2 = y_t(2);
    [q1, q2] = IK(x2,y2);
    x1 = l1*cos(q1);    y1 = l1*sin(q1);        
    
    % Save the results of end-effector position, velocity, acceleration
	x1_save(n) = x2;      % Save the end-effector position
	x2_save(n) = y2;
    dx1_save(n) = dy_t(1);
    dx2_save(n) = dy_t(2);
    ddx1_save(n) = ddy_t(1);
    ddx2_save(n) = ddy_t(2);
    f1_save(n) = f_t(1);
    f2_save(n) = f_t(2);
    q_save(n, 1) = q1*R2D;
    q_save(n, 2) = q2*R2D;
    
    % Calculate the coordinates of robot geometry for animation 
	Ax1 = [0, x1];  
	Ay1 = [0, y1];
   	Ax2 = [x1, x2];   
	Ay2 = [y1, y2];
   
    % Update the animation
	if rem(n,5) == 0
        set(p1,'X', Ax1, 'Y',Ay1);
        set(p2,'X', Ax2, 'Y',Ay2);
        drawnow
    end
    %pause(0.1);    
  
    % Save 1st and 2nd joint's location, (x1, y1) and (x2, y2)
	if rem(n,1) == 0
        x_save(n_trj) = x2;
        y_save(n_trj) = y2;
        n_trj = n_trj + 1;        
    end
    
    % Increase the iteration number
	n=n+1;    

end


% Plot the graph of end-effector
figure(2)
subplot(2,1,1)
plot(T, y_d(:,1), 'g', T, x1_save, 'b', 'linewidth', 0.8)
title('< end-effector x position >')
legend('demonstration', 'reproduction')
xlabel('t (sec)')
ylabel('x (m)')
xlim([0 tf])
subplot(2,1,2)
plot(T, y_d(:,2), 'g',  T, x2_save, 'b', 'linewidth', 0.8)
title('< end-effector y position >')
legend('demonstration', 'reproduction')
xlabel('t (sec)')
ylabel('y (m)')
xlim([0 tf])

figure(3)
subplot(2,1,1)
plot(T, dy_d(:,1), 'g', T, dx1_save, 'b', 'linewidth', 0.8)
title('< end-effector x velocity >')
legend('demonstration', 'reproduction')
xlabel('t (sec)')
ylabel('xdot (m)')
xlim([0 tf])
subplot(2,1,2)
plot(T, dy_d(:,2), 'g', T, dx2_save, 'b', 'linewidth', 0.8)
title('< end-effector y velocity >')
legend('demonstration', 'reproduction')
xlabel('t (sec)')
ylabel('ydot (m)')
xlim([0 tf])

figure(4)
subplot(2,1,1)
plot(T, ddy_d(:,1), 'g', T, ddx1_save, 'b', 'linewidth', 0.8)
title('< end-effector x acceleration >')
legend('demonstration', 'reproduction')
xlabel('t (sec)')
ylabel('xddot (m)')
xlim([0 tf])
subplot(2,1,2)
plot(T, ddy_d(:,2), 'g', T, ddx2_save, 'b', 'linewidth', 0.8)
title('< end-effector y acceleration >')
legend('demonstration', 'reproduction')
xlabel('t (sec)')
ylabel('yddot (m)')
xlim([0 tf])

figure(5)
subplot(2,1,1)
plot(T, f_ref(:,1), 'g', T, f1_save, 'b', 'linewidth', 0.8)
title('< Reference Force >')
legend('x (m/s)', 'y (m/s)')
xlabel('t (sec)')
ylabel('x (m)')
xlim([0 tf])
subplot(2,1,2)
plot(T, f_ref(:,2), 'g', T, f2_save, 'b', 'linewidth', 0.8)
title('< Force >')
legend('x (m/s^2)', 'y (m/s^2)')
xlabel('t (sec)')
ylabel('x (m)')
xlim([0 tf])

figure(6)
out = zeros(timestep, n_bfs, n_dofs);
sum = zeros(timestep);
for n=1:n_dofs
    for b=1:n_bfs
        out(:, b, n) = psi(:, b) * w(n, b).*z(:);
        sum = sum + psi(:, b);
    end
end
subplot(2,1,1)
plot(T, psi);
xlabel('t')
ylabel('\psi_i')
xlim([0 tf])
subplot(2,1,2)
plot(T, out(:, :, 1))
xlabel('t')
ylabel('\psi_i')
xlim([0 tf])

figure(7)
subplot(2,1,1)
plot(T, z);
subplot(2,1,2)
plot(T, sum);

figure(8)
% title('Animation')
grid
hold on
axis([-0.2 0.2 0 0.4]);
% Draw trajectory 
%plot(x_d(:, 1), x_d(:, 2), 'g', 'linewidth', 1)
plot(x_save, y_save, 'b')
xlabel('x (m)')
ylabel('y (m)')
%legend('demonstration', 'reproduction')


