clc;
clear all;
close all;
% disp('enter the IIR filter designspecifications');
rp= 1; % input('enter the passband ripple:');
rs= 50; % input('enter the stopband ripple:');
wp= 0.5; % input('enter the passband freq:');
ws= 0.7; % input('enter the stopband freq:');
fs= 2; % input('enter the sampling freq:');
w1=2*wp/fs;
w2=2*ws/fs;
[n,wn]=buttord(w1,w2,rp,rs,'s');
disp('Frequency response of IIR LPF is:');
[b,a]=butter(n,wn,'low','s');
w=0:.01:pi;
[h,om]=freqs(b,a,w);
m=20*log10(abs(h));
an=angle(h);
figure,subplot(2,1,1);plot(om/pi,m);
title('magnitude response of IIR filter is:');
xlabel('(a) Normalized freq. -->');
ylabel('Gain in dB-->');
subplot(2,1,2);
plot(om/pi,an);
title('phase response of IIR filter is:');
xlabel('(b) Normalized freq. -->');
ylabel('Phase in radians-->');