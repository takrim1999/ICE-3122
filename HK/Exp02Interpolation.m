interpolationFactor = 3;
frequency1 = 0.02;
frequency2 = 0.03;
signalLength = 100;
time = 0:signalLength-1;
interpolatedTime = 0:signalLength*interpolationFactor-1;
signal = sin(2*pi*frequency1*time)+sin(2*pi*frequency2*time);
interpolatedSignal = interp(signal(1:signalLength),interpolationFactor);
subplot(2,1,1)
stem(time,signal(1:signalLength));
title('Signal');
subplot(2,1,2)
% length(interpolatedSignal)
% length(interpolatedTime)
stem(interpolatedTime,interpolatedSignal);
title('Interpolated Signal');
