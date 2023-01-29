clear all 
close all 
clc
%% Signal x(t)
te=1e-4;
fe=1/te;
t=0:te:5-te;
N=length(t);
f1=500;
f2=400;
f3=50;
x=sin(2*pi*f1*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t);
plot(t,x)
title('signal x(t)')
%% Spectre avec Te=0.0001
 fshift = (-N/2:N/2-1)*(fe/N);
 y = fft(x);
 plot(fshift, fftshift(abs(y)/N)*2);
 title('spectre');

 %%  Spectre avec Te=0.0005
 te =5e-4 ;
 y = fft(x);
 plot(fshift, fftshift(abs(y)/N)*2);
 title('le spectre');

 
%% fonction de transmittance 
te = 1e-4 ;
fe = 1/te ;
t = 0:te:5 ;
N = length(t);
f = (0:N-1)*(fe/N);
K = 1 ;
wc1 = 50 ;
w = 2*pi*f ; 

H = (K*1j*w/wc1)./(1+1j*w/wc1) ;
Module_H = abs(H);

   plot(w,Module_H, 'r')



%%
te = 1e-4 ;
fe = 1/te ;
N = length(t);
f = (0:N-1)*(fe/N);
fshift = (-N/2:N/2-1)*fe/N;
k=1;
w = 2*pi*f; 
fc1 = 30; 
fc2 =200;
fc3=300;
H1 = (k*1j*w/fc1) ./ (1 + 1j*w/fc1);
H2 = (k*1j*w/fc2) ./ (1 + 1j*w/fc2);
H3 = (k*1j*w/fc3) ./ (1 + 1j*w/fc3);
G1 = 20*log(abs(H1));
G2 = 20*log(abs(H2));
G3 = 20*log(abs(H3));
semilogx(f,G1,'g',f,G2,'r',f,G3,'b')
legend('fc = 200 rad/s', 'fc = 500 rad/s', 'fc = 1000 rad/s');
grid on 

%%
yt1 = y.*H1 ;
yt2 = y.*H2 ;
yt3 = y.*H3 ;
YT1 = ifft(yt1,"symmetric");
YT2 = ifft(yt2,"symmetric");
YT3 = ifft(yt3,"symmetric");
YT1t = fft(YT1);
YT2t = fft(YT2);
YT3t = fft(YT3);
subplot(2,2,1) 
plot(fshift,2*fftshift(abs(y))/N);
title('spectre d amplitude');
subplot(2,2,2)
plot(fshift,2*fftshift(abs(YT1t))/N,'r');
title('fc=200rad/s');
subplot(2,2,3)
plot(fshift,2*fftshift(abs(YT2t))/N,'bl');
title('fc=500rad/s');
subplot(2,2,4)
plot(fshift,2*fftshift(abs(YT3t))/N,'g');
title('fc=1000rad/s');

%% Choix d'une pulsation de coupure optimale
K = 1 ;
f_op=1000;
wc_op=2*pi*f_op;
w = 2*pi*f ; 
H = (K*1j*w/wc_op)./(1+1j*w/wc_op) ;

yt = y.*H ;
YT = ifft(yt,"symmetric");
YTt = fft(YT);
plot(fshift,2*fftshift(abs(YTt))/N,'r');
title('spectre filtré ');


