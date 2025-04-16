
%% FallingPaperUserInput.m
% Simulation of a paper sheet falling with flutter & flipping.
% User inputs: drop height H and initial angle theta.

clear; clc; close all;

%% 1. USER INPUTS
H = input('Enter drop height H (m): ');
theta_deg = input('Enter initial tilt angle θ from vertical (deg): ');
theta = theta_deg * pi/180;

%% 2. PARAMETERS & INITIAL CONDITIONS

% Sheet geometry & mass
Lx = 0.15;      % A4 width (m)
Ly = 0.2;       % A4 height (m)
m  = 0.015;     % mass (kg)
rho_air = 1.225; Cd = 1.2;  % air density & drag coeff
g = [0;0;-9.81];            % gravity (m/s^2)
C_tau = 0.2;                % aero torque coeff

% Inertia tensor
Ix = (1/12)*m*Ly^2;
Iy = (1/12)*m*Lx^2;
Iz = (1/12)*m*(Lx^2+Ly^2);
I  = diag([Ix Iy Iz]);

% Time stepping
dt = 0.005;

% Initial state
r     = [0; 0; H];               % start at height H
v     = [0.15; -0.05; 0];        % small lateral push
% Tilt the sheet by θ about the y-axis:
R = [ cos(theta), 0, sin(theta);
           0    , 1,     0    ;
     -sin(theta), 0, cos(theta)];
omega = [2; 1; 0.6];             % initial angular velocity

% Preallocate history (max length)
maxSteps = ceil(5*sqrt(2*H/9.81)/dt);
r_hist         = zeros(3, maxSteps);
v_hist         = zeros(3, maxSteps);
omega_hist     = zeros(3, maxSteps);
R_hist         = zeros(3,3,maxSteps);
time_hist      = zeros(1, maxSteps);

%% 3. EULER INTEGRATION UNTIL GROUND CONTACT

i = 0;
while true
    i = i + 1;
    % Store
    r_hist(:,i)     = r;
    v_hist(:,i)     = v;
    omega_hist(:,i) = omega;
    R_hist(:,:,i)   = R;
    time_hist(i)    = (i-1)*dt;
    
    % Stop if hit ground
    if r(3) <= 0
        t_ground = time_hist(i);
        fprintf('Ground contact at t = %.3f s\n', t_ground);
        break
    end
    
    %% TRANSLATION
    n_world = R * [0;0;1];
    if norm(v)<1e-8, v_hat=[0;0;0]; else v_hat=v/norm(v); end
    A_proj = Lx*Ly * abs(dot(n_world, -v_hat));
    F_drag = -0.5 * rho_air * Cd * A_proj * norm(v) * v;
    a = g + F_drag/m;
    v = v + a*dt;
    r = r + v*dt;
    
    %% ROTATION
    tau_aero = -C_tau * cross(n_world, v) + 0.002*randn(3,1);
    Lb = I*omega;
    domega = I \ (tau_aero - cross(omega, Lb));
    omega = omega + domega*dt;
    Omega_mat = [  0,       -omega(3),  omega(2);
                omega(3),     0,      -omega(1);
               -omega(2),  omega(1),     0     ];
    R = R + Omega_mat * R * dt;
    [U,~,V] = svd(R); R = U*V';  % re-orthonormalize
end

% Trim history arrays
N = i;
r_hist         = r_hist(:,1:N);
v_hist         = v_hist(:,1:N);
omega_hist     = omega_hist(:,1:N);
R_hist         = R_hist(:,:,1:N);
time_hist      = time_hist(1:N);

%% 4. SLOW-MOTION ANIMATION

sheet_corners = [ Lx/2,  Ly/2, 0;
                 -Lx/2,  Ly/2, 0;
                 -Lx/2, -Ly/2, 0;
                  Lx/2, -Ly/2, 0 ]';

figure('Color','w');
axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([0 H]);
view(35,20); hold on;
h = patch('Faces',[1 2 3 4],'Vertices',zeros(4,3), ...
          'FaceColor','b','FaceAlpha',0.6);

for k = 1:N
    Rk = R_hist(:,:,k);
    rk = r_hist(:,k);
    pts = (Rk*sheet_corners)' + rk';
    set(h,'Vertices',pts);
    title(sprintf('t = %.2f s, z = %.2f m', time_hist(k), rk(3)));
    drawnow; pause(0.05);
end

%% 5. COMPUTE & PLOT DERIVED QUANTITIES

% Preallocate
A_proj_hist     = zeros(1,N);
Lmag_hist       = zeros(1,N);
Fdrag_mag_hist  = zeros(1,N);
v_mag_hist      = zeros(1,N);
omega_mag_hist  = zeros(1,N);
z_hist          = zeros(1,N);

for k = 1:N
    Rk = R_hist(:,:,k);
    vk = v_hist(:,k);
    % projected area
    nk = Rk*[0;0;1];
    if norm(vk)<1e-8, v_hat=[0;0;0]; else v_hat=vk/norm(vk); end
    A_proj_hist(k) = Lx*Ly * abs(dot(nk, -v_hat));
    % angular momentum magnitude
    Lb = I * omega_hist(:,k);
    Lw = Rk * Lb;
    Lmag_hist(k) = norm(Lw);
    % drag magnitude
    Fd = -0.5 * rho_air * Cd * A_proj_hist(k) * norm(vk) * vk;
    Fdrag_mag_hist(k) = norm(Fd);
    % velocity & angular velocity magnitudes
    v_mag_hist(k)     = norm(vk);
    omega_mag_hist(k) = norm(omega_hist(:,k));
    % z-position
    z_hist(k) = r_hist(3,k);
end

% Plot Z-Position
figure;
plot(time_hist, z_hist, 'LineWidth',2);
grid on; xlabel('Time (s)'); ylabel('Z (m)');
title('Height (Z) vs. Time');

% Plot Projected Area
figure;
plot(time_hist, A_proj_hist, 'LineWidth',2);
grid on; xlabel('Time (s)'); ylabel('A_{proj} (m^2)');
title('Projected Area vs. Time');

% Plot Angular Momentum
figure;
plot(time_hist, Lmag_hist, 'LineWidth',2);
grid on; xlabel('Time (s)'); ylabel('|L| (kg·m^2/s)');
title('Angular Momentum Magnitude vs. Time');

% Plot Drag Force
figure;
plot(time_hist, Fdrag_mag_hist, 'LineWidth',2);
grid on; xlabel('Time (s)'); ylabel('|F_{drag}| (N)');
title('Drag Force Magnitude vs. Time');

% Plot Velocity Magnitude
figure;
plot(time_hist, v_mag_hist, 'LineWidth',2);
grid on; xlabel('Time (s)'); ylabel('|\bfv| (m/s)');
title('Velocity Magnitude vs. Time');

% Plot Angular‐Velocity Magnitude
figure;
plot(time_hist, omega_mag_hist, 'LineWidth',2);
grid on; xlabel('Time (s)'); ylabel('|\bomega| (rad/s)');
title('Angular‐Velocity Magnitude vs. Time');

%% 6. GENERATE TABLE FOR EACH 0.2 SECONDS

% Desired time intervals
t_samples = (0:0.02:time_hist(end))';
nS = numel(t_samples);

% Preallocate sampled arrays
A_sampled     = zeros(nS,1);
L_sampled     = zeros(nS,1);
Fd_sampled    = zeros(nS,1);
v_sampled     = zeros(nS,1);
omega_sampled = zeros(nS,1);
z_sampled     = zeros(nS,1);

for j = 1:nS
    [~, idx] = min(abs(time_hist - t_samples(j)));
    A_sampled(j)     = A_proj_hist(idx);
    L_sampled(j)     = Lmag_hist(idx);
    Fd_sampled(j)    = Fdrag_mag_hist(idx);
    v_sampled(j)     = v_mag_hist(idx);
    omega_sampled(j) = omega_mag_hist(idx);
    z_sampled(j)     = z_hist(idx);
end

% Create table
T = table( ...
    t_samples, ...
    z_sampled, ...
    A_sampled, ...
    L_sampled, ...
    Fd_sampled, ...
    v_sampled, ...
    omega_sampled, ...
    'VariableNames', { ...
      'Time_s', ...
      'Z_m', ...
      'ProjectedArea_m2', ...
      'AngularMomentum_kgm2s', ...
      'DragForce_N', ...
      'Velocity_m_s', ...
      'AngVel_rad_s'});

% Display
disp('--- Sampled Data Every 0.2 Seconds ---');
disp(T);

% Export to CSV
writetable(T,'falling_sheet_data4sduf.csv');
fprintf('Table saved to falling_sheet_data.csv\n');
