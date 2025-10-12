clc;clear
%% Introduction: double pendulum (approx)
%Physics parametres
L1 = 1; L2 = 0.5; m1 = 0.75; m2 = 0.5; g = 9.81;

%Condition
t_init = [0 10];
initial_condition = [20 * pi / 180; 0; 0; 0];
initial_condition_aprox = [20 * pi / 180; 0; 0; 0];

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t, z] = ode45(@(t, theta)equations(t, theta, L1, L2, g, m1, m2), t_init, initial_condition, odeset);
[t_aprox, z_aprox] = ode45(@(t, theta)equations_aprox(t, theta, L1, L2, g, m1, m2), t_init, initial_condition_aprox, odeset);

%% Animation: with a little tricky

if (length(t) >= length(t_aprox))
    nFrames = length(t_aprox);
else
    nFrames = length(t);
end

% Energy vector
Ek = zeros(length(t),1);
Ep = zeros(length(t),1);
E_total = zeros(length(t),1);

Ek_aprox = zeros(length(t_aprox),1);
Ep_aprox = zeros(length(t_aprox),1);
E_total_aprox = zeros(length(t_aprox),1);

% Energy calculus
for j = 1:length(t)
    Ek(j) = 0.5*(m1 + m2)*L1^2*z(j,2)^2 + 0.5*m2*L2^2*z(j,4)^2 + m2*L1*L2*z(j,2)*z(j,4)*cos(z(j,1)-z(j,3));
    Ep(j) = -(m1 + m2)*g*L1*cos(z(j,1)) - m2*g*L2*cos(z(j,3));
    E_total(j) = Ek(j) + Ep(j);
end

for j = 1:length(t_aprox)
    Ek_aprox(j) = 0.5*(m1 + m2)*L1^2*z_aprox(j,2)^2 + 0.5*m2*L2^2*z_aprox(j,4)^2 + m2*L1*L2*z_aprox(j,2)*z_aprox(j,4)*cos(z_aprox(j,1)-z_aprox(j,3));
    Ep_aprox(j) = -(m1 + m2)*g*L1*cos(z_aprox(j,1)) - m2*g*L2*cos(z_aprox(j,3));
    E_total_aprox(j) = Ek_aprox(j) + Ep_aprox(j);
end

%Animation
figure(1);
for j = 1:nFrames
    clf;
    x = L1 * sin(z(j,1));
    y = - L1 * cos(z(j,1));
    x1 = x + L2 * sin(z(j,3));
    y1 = y - L2 * cos(z(j,3));

    x_aprox = L1 * sin(z_aprox(j,1));
    y_aprox = - L1 * cos(z_aprox(j,1));
    x1_aprox = x_aprox + L2 * sin(z_aprox(j,3));
    y1_aprox = y_aprox - L2 * cos(z_aprox(j,3));

    hold on
    % Plot exact pendulum
    plot([0 x], [0 y], 'k-', 'LineWidth', 2);
    plot([x x1], [y y1], 'k-', 'LineWidth', 2);
    plot(x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
    plot(x1, y1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','g');
   


    % Plot approximated pendulum (shifted right)
    plot([4 4+x_aprox], [0 y_aprox], 'k-', 'LineWidth', 2);
    plot([x_aprox+4 x1_aprox+4], [y_aprox y1_aprox], 'k-', 'LineWidth', 2);
    plot(x_aprox+4, y_aprox, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
    plot(x1_aprox+4, y1_aprox, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','g');

    % Ground line
    plot([-2 7], [0 0], 'k', 'LineWidth', 2);

    axis equal; axis([-2, 7, -3, 1]);
    drawnow;
end
%% Total Energy graphs

figure(2);

subplot(3,1,1);
plot(t_aprox, Ek_aprox, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Kinetic energy (J)');
title('Kinetic Energy');

subplot(3,1,2);
plot(t_aprox, Ep_aprox, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Potential energy (J)');
title('Potential Energy');

subplot(3,1,3);
plot(t_aprox, E_total_aprox, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Total energy (J)');
title('Energy Conservation (Approximation)');
grid on;

figure(3);

subplot(3,1,1)
plot(t, Ek, 'LineWidth', 2)
xlabel('Time (s)');
ylabel('Total kinetic energy');
title('Kinetic Energy');

subplot(3,1,2)
plot(t, Ep, 'LineWidth', 2)
xlabel('Time (s)');
ylabel('Potencial energy');
title('Potencial Energy');

subplot(3,1,3)
plot(t(1:length(t)), E_total, 'LineWidth', 2);
xlabel('Time (s)');ylabel('Total energy (J)');
title('Energy conservation');
grid on;

%% Functions
function dy = equations(t, theta, L1, L2, g, m1, m2)
    
    theta1 = theta(1);
    dtheta1 = theta(2);
    theta2 = theta(3);
    dtheta2 = theta(4);

    A = [(m1 + m2) * L1^2, m2 * L1 * L2 * cos(theta1 - theta2);
        m2 * L1 * L2 * cos(theta1 - theta2), m2 * L2^2];

    B = [-m2 * L1 * L2 * sin(theta1 - theta2) * dtheta2^2 - (m1 + m2) * g * L1 * sin(theta1);
     m2 * L1 * L2 * sin(theta1 - theta2) * dtheta1^2 + m2 * g * L2 * sin(theta2)];

    sol = A \ B;

    dtheta1_sol = sol(1);
    dtheta2_sol = sol(2);
    
    
    dy = zeros(4,1);
    dy(1) = dtheta1;
    dy(2) = dtheta1_sol;
    dy(3) = dtheta2;
    dy(4) = dtheta2_sol;

end

function dy_aprox = equations_aprox(t, theta, L1, L2, g, m1, m2)
    
    theta1 = theta(1);
    dtheta1 = theta(2);
    theta2 = theta(3);
    dtheta2 = theta(4);

    A = [(m1 + m2) * L1, m2 * L2;
         m2 * L1,        m2 * L2];

    B = [- (m1 + m2) * g * theta1;
         - m2 * g * theta2];

    sol = A \ B;

    dtheta1_aprox = sol(1);
    dtheta2_aprox = sol(2);

    dy_aprox = zeros(4,1);
    dy_aprox(1) = dtheta1;
    dy_aprox(2) = dtheta1_aprox;
    dy_aprox(3) = dtheta2;
    dy_aprox(4) = dtheta2_aprox;
    
end