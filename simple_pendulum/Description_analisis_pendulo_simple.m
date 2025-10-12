clc;clear
% A continuación describiremos la mecánica de uno de los tres sistemas que
% pueden los físicos resolver en su totalidad: El oscilador armónico lineal (OAL), cuya
% ecuación de movimiento es: x'' + (g/L) sin(x + x0)  = 0, donde L es la longitug, la g es
% la gravedad. La solución que nosotros
% Nota (importante!!)
%  - Estamos en campos conservativos (donde F = - ∇U).
%  - El pto. de equilibrio lo cogemos como origen de coordenadas y de
%  potencial x0 = 0. 
%  - La resolución de la ecuaciación diferencial se puede hacer de manera
%  análitica: x(t) = x0*cos(omega*t - delta), donde omega = sqrt(g/L) y
%  delta es un desfase. Sin embargo, para hacer más instructivo el código,
%  partiremos de la ecuación diferencial con los siguientes supuestos.

g = 9.8;
L = 2;
x01 = 0; %CC

syms x(t) %Haremos la ecuación con variables simbólicas
derivada_segunda = diff(x, t, 2);

eqn = derivada_segunda + g/L*sin(x - x01) == 0;

%La resolución de esta ecuación la podemos hacer con métodos numéricos o
%con las funciones de Matlab. Haremos un análisis. 
%
%Con las funciones de Matlab: 

% Dsolve: no puede resolver sistemas no lineales por eso linealizamos el
% sistema con ángulo pequeños: 
eqn1 = derivada_segunda + g/L*(x - x01) == 0;
solucion_dsolve = dsolve(eqn1);
disp(solucion_dsolve)

% ode45 no usa variables simbólicas, asi que, reescribiremos la ecuación
%como un sistema de ecuaciones
pendulo = @(t, x) [x(2); -(g/L) * sin(x(1) - x01)];
t_init = [0 10]; % tiempo 
x_init = [pi/4; 0]; %Ángulo inicial; velocidad angular inicial
[t, x]= ode45(pendulo, t_init, x_init);

% Con los métodos numéricos:

h = 0.01; %Paso
N = 1000; %Número de iteraciones
x0 = pi/4;v0 = 0; %Condiciones iniciales
omega = sqrt(g/L);

%Comparación de los cuatro métodos con la solución del sistema original.

[t1, x1, v1] = euler_implicito_osc(0, x0, v0, h, N, omega);
[t2, x2, v2] = heun_osc(0, x0, v0, h, N, omega);
[t3, x3, v3] = RG4(0, x0, v0, h, N, omega);

plot(t, x(:,1), 'r', t1, x1, 'b--', t2, x2, 'g-.', t3, x3, 'k:');
legend('ODE45', 'Euler Implícito', 'Heun', 'Runge-Kutta 4');
xlabel('Tiempo'); ylabel('Desplazamiento');
title('Comparación de Métodos Numéricos');

%Método de Euler implícito: requiere resolver un sistema de ecuaciones en
%cada paso en nuestro caso, la ecuación cambiará: dx/dt = v; dv/dt = -
%g/L*(x-x0)

function [t, x, v] = euler_implicito_osc(t0, x0, v0, h, N, omega)
    x = zeros(1, N+1);
    v = zeros(1, N+1);
    t = t0:h:(t0 + N*h);
    x(1) = x0;
    v(1) = v0;

   for n = 1:N
        % Definimos el sistema no lineal:
        F = @(y) [y(1) - x(n) - h*y(2); y(2) - v(n) + h*(omega^2)*sin(y(1))];
        % Usamos como estimación inicial el valor anterior:
        y0 = [x(n); v(n)];
        % Resolvemos el sistema no lineal:
        y_next = fsolve(F, y0);
        x(n+1) = y_next(1);
        v(n+1) = y_next(2);
    end
end

%Método de Heun: que es el método de Euler mejorado (con una corrección)

function [t, x, v] = heun_osc(t0, x0, v0, h, N, omega)
    x = zeros(1, N+1);
    v = zeros(1, N+1);
    t = t0:h:(t0 + N*h);
    x(1) = x0;
    v(1) = v0;

    for n = 1:N
        % Hacemos una estimación lineal del resultado a través del Euler
        % Explícito. Llamado: predictor
        x_pred = x(n) + h * v(n);
        v_pred = v(n) - h * omega^2 * sin(x(n));

        % Ahora corregimos la información de "dirección" que nos da el
        % predictor. Llamado: corrector
        x(n+1) = x(n) + (h/2) * (v(n) + v_pred);
        v(n+1) = v(n) - (h/2) * omega^2 * (sin(x(n)) + sin(x_pred));
    end
end

%Método de RK4 (es en el que está basado ode45):

function [t, x, v] = RG4(t0, x0, v0, h, N, omega)
    x = zeros(1, N+1);
    v = zeros(1, N+1);
    t = t0:h:(t0 + N*h);
    x(1) = x0;
    v(1) = v0;

    for n = 1:N
        % k1
        k1x = v(n);
        k1v = -omega^2 * sin(x(n));

        % k2
        k2x = v(n) + 0.5*h*k1v;
        k2v = -omega^2 * sin((x(n) + 0.5*h*k1x));

        % k3
        k3x = v(n) + 0.5*h*k2v;
        k3v = -omega^2 * sin((x(n) + 0.5*h*k2x));

        % k4
        k4x = v(n) + h*k3v;
        k4v = -omega^2 * sin((x(n) + h*k3x));

        % Actualización
        x(n+1) = x(n) + (h/6)*(k1x + 2*k2x + 2*k3x + k4x);
        v(n+1) = v(n) + (h/6)*(k1v + 2*k2v + 2*k3v + k4v);
    end
end
% 
% El método implícito de Euler es numéricamente disipativo
% Esto significa que reduce la energía del sistema con cada paso, aunque el modelo físico no tenga ninguna fuerza de fricción.
% En sistemas oscilatorios como el péndulo, esto se traduce en que la amplitud de las oscilaciones disminuye artificialmente con el tiempo.
% Es una propiedad intrínseca del método, no del sistema que estás modelando.
% El método implícito de Euler prioriza la estabilidad sobre la precisión.
% En cada paso, resuelve una ecuación que "mira hacia adelante", lo que tiende a suavizar las oscilaciones.
% Esto es útil en sistemas rígidos o con alta frecuencia, pero en sistemas conservativos (como el péndulo sin fricción), no representa la física real.