%% Ecuación I en tiempo contínuo
t = 0:0.001:50;
y = 2*exp(-0.1*t).*sin((2*t)/3);
disp("Máximo global de la función: "+max(y));
disp("Mínimo global de la función: "+min(y));
TF = islocalmax(y);
TF2 = islocalmin(y);
figure;
plot(t,y,t(TF),y(TF),'*g',t(TF2),y(TF2),'*r');
axis padded;
xlabel("Tiempo(t)");
ylabel("Señal y");
grid on; grid minor;

%% Ecuacion 2 en tiempo contínuo

t = 0:0.001:50;
y = 2*exp(0.1*t).*sin((2*t)/3);
disp("Máximo global de la función: "+max(y));
disp("Mínimo global de la función: "+min(y));
TF = islocalmax(y);
TF2 = islocalmin(y);
figure;
plot(t,y,t(TF),y(TF),'*g',t(TF2),y(TF2),'*r');
axis padded;
xlabel("Tiempo(segundos)");
ylabel("Señal y");
grid on; grid minor;

%% FFT de la ecuación I en tiempo contínuo

Fs = 10; % Frecuencia de muestreo
Ts = 1/Fs; % Periodo de muestreo
L = 20; %Muestras
t=(0:L-1).*Ts;
y = 2*exp(-0.1*t).*sin((2*t)/3);
fourier = fft(y);
figure;
plot(abs(fourier));
axis padded;
xlabel("Frecuencia (Hz)");
ylabel("Señal y");
grid on; grid minor;


%% FFT de la ecuación II en tiempo contínuo
Fs = 10; % Frecuencia de muestreo
Ts = 1/Fs; % Periodo de muestreo
L = 20;%Muestras
t=(0:L-1).*Ts;
y = 2*exp(0.1*t).*sin((2*t)/3);
fourier = fft(y);
figure;
plot(abs(fourier));
axis padded;
xlabel("Frecuencia (Hz)");
ylabel("Señal y");
grid on; grid minor;

%% Señal en tiempo discreto
n = -4:46;
y = 2*exp(-0.1*t).*sin((2*t)/3);
figure;
stem(n,y,"filled");
axis padded;
xlabel("Numero de muestras(n)");
ylabel("Señal y");
grid on; grid minor;

%% Punto 9. Señal I en tiempo discreto
n=-5:60;
y1 = (0).*(n>=-5 & n<=3) + (n.^2 + 0.5.*n).*(n>=4 & n<=30)  ...
+ (1.5).*(n>=31 & n<=40) + (0).*(n>=41 & n<=60); 
figure;
stem(n,y1,"filled");
axis padded;
title("Señal discreta I: x1 [n]");
xlabel("n");
ylabel("x1 [n]");
grid on; grid minor;

%% DFT de la señal I en tiempo discreto
n = -5:60; 
y = (0).*(n>=-5 & n<=3) + (n.^2 + 0.5.*n).*(n>=4 & n<=30)  ...
+ (1.5).*(n>=31 & n<=40) + (0).*(n>=41 & n<=60); 

y1 = fft(y); % DFT de x
m = abs(y1); % Magnitud
%y1(m<1e-6) = 0;
p = unwrap(angle(y1)); % Fase
f = (0:length(y1)-1)*100/length(y1); % Vector de frecuencia

subplot(2,1,1)
stem(f,m)
title('Magnitud')
axis padded;
grid on; grid minor;

subplot(2,1,2)
stem(f,p*180/pi)
title('Fase')
axis padded;
grid on; grid minor;

%% Punto 9. Señal II en tiempo discreto
n=-15:50;
y2 = (2*cos(2.*n*pi)).*(n>=-15 & n<=8) ...
+ (0.5).*(n>=9 & n<=17) + (exp((n.*pi)/4)...
.*sin((2/3).*n)).*(n>=18 & n<=40) + (0).*(n>=41 & n<=50); 
figure;
stem(n,y2,"filled");
axis padded;
title("Señal discreta II: x2[n]");
xlabel("n");
ylabel("x2 [n]");
grid on; grid minor;

%% DFT de la señal II en tiempo discreto
n = -15:50; 
y = (2*cos(2.*n*pi)).*(n>=-15 & n<=8) ...
 + (0.5).*(n>=9 & n<=17) + (exp((n.*pi)/4)...
 .*sin((2/3).*n)).*(n>=18 & n<=40) + (0).*(n>=41 & n<=50);

y1 = fft(y); % DFT de y
m = abs(y1); % Magnitud
p = unwrap(angle(y1)); % Fase
f = (0:length(y1)-1)*100/length(y1); % Vector de frecuencia

subplot(2,1,1)
stem(f,m)
title('Magnitud')
axis padded;
grid on; grid minor;

subplot(2,1,2)
stem(f,p*180/pi)
title('Fase')
axis padded;
grid on; grid minor;