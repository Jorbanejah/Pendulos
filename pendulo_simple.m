% A continuación describiremos la mecánica de uno de los tres sistemas que
% pueden los físicos resolver en su totalidad: El oscilador armónico lineal (OAL), cuya
% ecuación de movimiento es: x'' + (k/M) (x - x0)  = 0, donde M es la masa, la K es
% la una constante dependiente de la segunda derivada del potencial (U) y la x es el recorrido que
% seguirá el sistema.
% Nota (importante!!)
%  - Estamos en campos conservativos (donde F = - ∇U).
%  - El pto. de equilibrio lo cogemos como origen de coordenadas y de
%  potencial x0 = 0. 
%  - En el desarrollo de Taylor del potencial entorno al pto. de equilibrio, nos
%  quedamos en el segundo orden. Con lo que se demuestra que U(x) = 1/2*k/M
%  - La resolución de la ecuaciçon diferencial se puede hacer de manera
%  análitica: x(t) = x0*cos(omega*t - delta), donde omega = sqrt(k/M) y
%  delta es un desfase. Sin embargo, para hacer más instructivo el código,
%  partiremos de la ecuación diferencial con los siguientes supuestos.

M = 200;
k = 2;
x0 = 0;

syms x(t) %Haremos la ecuación con variables simbólicas
derivada_primera = diff(x,t);
derivada_segunda = diff(x, t, 2);

eqn = derivada_segunda + k/M*(derivada_primera - x0) == 0;

% Condiciones iniciales:
cond1 = x(0) == 1;


%La resolución de esta ecuacion la podemos hacer con métodos numéricos o
%con las funciones de Matlab como ode45, dsolve... al utilizar el lenguaje
%simbólico nos tenemos que quedar con dsolve()

solucion = dsolve(eqn, cond1);
disp(solucion)

% Para poder ilustrarlo, cogeremos el ejemplo tipicamente escogido en estos casos: Péndulo simple.
