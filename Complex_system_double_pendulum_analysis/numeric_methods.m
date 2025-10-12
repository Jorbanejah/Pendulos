%Probar el ode45, ode15s, y el metodo de Verlet

%Physics parametres
L1 = 1; L2 = 0.5; m1 = 0.75; m2 = 0.5; g = 9.81;

%Condition
angle = linspace(1, 45, 44);
t_init = [0 10];

opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
%% Ode45
E_total_average_ode45 = zeros(length(angle),2);
for i = 1:length(angle)
    initial_condition = [angle(i) * pi / 180; 0; 0; 0];
    [t_ode45, z_ode45] = ode45(@(t, theta)equations(t, theta, L1, L2, g, m1, m2), t_init, initial_condition, odeset);
    
    Ek_ode45 = zeros(length(t_ode45),1);
    Ep_ode45 = zeros(length(t_ode45),1);
    E_total_ode45 = zeros(length(t_ode45),1);
    for j = 1:length(t_ode45)
        Ek_ode45(j) = 0.5*(m1 + m2)*L1^2*z_ode45(j,2)^2 + 0.5*m2*L2^2*z_ode45(j,4)^2 + m2*L1*L2*z_ode45(j,2)*z_ode45(j,4)*cos(z_ode45(j,1)-z_ode45(j,3));
        Ep_ode45(j) = -(m1 + m2)*g*L1*cos(z_ode45(j,1)) - m2*g*L2*cos(z_ode45(j,3));
        E_total_ode45(j) = Ek_ode45(j) + Ep_ode45(j);
    end
    
    E_total_average_ode45(i,1) = max(E_total_ode45);
    E_total_average_ode45(i,2) = min(E_total_ode45);

end

figure(1)
plot(angle, E_total_average_ode45(:,1), 'b-')
hold on
plot(angle,E_total_average_ode45(:,2), 'r-')
xlabel('Initial angle (degrees)'); ylabel('Average Total Energy');
legend('Max Energy', 'Min Energy')
grid on

energy_range = E_total_average_ode45(:,1) - E_total_average_ode45(:,2);
figure(2)
plot(angle, energy_range, 'k-', 'LineWidth', 2)
xlabel('Initial angle (degrees)');
ylabel('Energy Range (J)');
title('Energy Spread vs Initial Angle');
grid on
%% Ode15s

[t_ode15s, z_ode15s] = ode15s(@(t, theta)equations(t, theta, L1, L2, g ,m1, m2), t_init, initial_condition, odeset);
initial_condition = [angle(i) * pi / 180; 0; 0; 0];
%% Simplectic Verlet integrator
% Parámetros de integración Verlet

dt = 1e-3;              % Paso de tiempo fijo
tspan = 0:dt:10;        % Vector de tiempo
N = length(tspan);

z_verlet = zeros(N,4);  % [theta1, dtheta1, theta2, dtheta2]
z_verlet(1,:) = initial_condition';  % condiciones iniciales

% Calcular aceleración inicial
accel = acceleration(z_verlet(1,:), L1, L2, g, m1, m2);

for i = 1:N-1
    % Posición siguiente
    z_verlet(i+1,1) = z_verlet(i,1) + z_verlet(i,2)*dt + 0.5*accel(1)*dt^2;
    z_verlet(i+1,3) = z_verlet(i,3) + z_verlet(i,4)*dt + 0.5*accel(2)*dt^2;

    % Aceleración en la nueva posición
    accel_new = acceleration([z_verlet(i+1,1), z_verlet(i,2), z_verlet(i+1,3), z_verlet(i,4)], L1, L2, g, m1, m2);

    % Velocidad siguiente
    z_verlet(i+1,2) = z_verlet(i,2) + 0.5*(accel(1)+accel_new(1))*dt;
    z_verlet(i+1,4) = z_verlet(i,4) + 0.5*(accel(2)+accel_new(2))*dt;

    % Actualizar aceleración
    accel = accel_new;
end

% Energía con Verlet
Ek_verlet = zeros(N,1);
Ep_verlet = zeros(N,1);
E_total_verlet = zeros(N,1);

for i = 1:N
    Ek_verlet(i) = 0.5*(m1+m2)*L1^2*z_verlet(i,2)^2 + 0.5*m2*L2^2*z_verlet(i,4)^2 + m2*L1*L2*z_verlet(i,2)*z_verlet(i,4)*cos(z_verlet(i,1)-z_verlet(i,3));
    Ep_verlet(i) = -(m1+m2)*g*L1*cos(z_verlet(i,1)) - m2*g*L2*cos(z_verlet(i,3));
    E_total_verlet(i) = Ek_verlet(i) + Ep_verlet(i);
end

figure;
subplot(3,1,1); plot(tspan,Ek_verlet); title('Kinetic Energy (Verlet)');
subplot(3,1,2); plot(tspan,Ep_verlet); title('Potential Energy (Verlet)');
subplot(3,1,3); plot(tspan,E_total_verlet); title('Total Energy (Verlet)'); grid on;

%% Auxiliary function: acceleration
function acc = acceleration(state, L1, L2, g, m1, m2)
    theta1 = state(1); dtheta1 = state(2);
    theta2 = state(3); dtheta2 = state(4);

    A = [(m1 + m2) * L1^2, m2 * L1 * L2 * cos(theta1 - theta2);
         m2 * L1 * L2 * cos(theta1 - theta2), m2 * L2^2];

    B = [-m2 * L1 * L2 * sin(theta1 - theta2) * dtheta2^2 - (m1 + m2) * g * L1 * sin(theta1);
          m2 * L1 * L2 * sin(theta1 - theta2) * dtheta1^2 + m2 * g * L2 * sin(theta2)];

    sol = A \ B;
    acc = sol(:)';  % [ddtheta1, ddtheta2]
end

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
