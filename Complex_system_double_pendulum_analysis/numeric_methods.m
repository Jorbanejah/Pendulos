%Probar el ode45, ode15s, y el metodo de Verlet

%% Ode45
%% Ode15s
%% Integrador simpléctico Verlet
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
    Ek_verlet(i) = 0.5*(m1+m2)*L1^2*z_verlet(i,2)^2 + ...
                   0.5*m2*L2^2*z_verlet(i,4)^2 + ...
                   m2*L1*L2*z_verlet(i,2)*z_verlet(i,4)*cos(z_verlet(i,1)-z_verlet(i,3));
    Ep_verlet(i) = -(m1+m2)*g*L1*cos(z_verlet(i,1)) - m2*g*L2*cos(z_verlet(i,3));
    E_total_verlet(i) = Ek_verlet(i) + Ep_verlet(i);
end

figure;
subplot(3,1,1); plot(tspan,Ek_verlet); title('Kinetic Energy (Verlet)');
subplot(3,1,2); plot(tspan,Ep_verlet); title('Potential Energy (Verlet)');
subplot(3,1,3); plot(tspan,E_total_verlet); title('Total Energy (Verlet)'); grid on;

%% Función auxiliar: aceleraciones
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
