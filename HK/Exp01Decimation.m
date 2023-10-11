% equation: signal = sin(2*pi*frequency1*time)+sin(2*pi*frequency2*time);
decimationFactor = 5; % input('Down Sample Rate: ');
signalLength = 100; % input('Length of the input signal: ');
frequency1 = 0.02; % input('Frequency of First Sinusidal Signal: ');
frequency2 = 0.03; %input('Frequency of Second Sinusidal Signal: ');
time = 0:signalLength-1;
signal = sin(2*pi*frequency1*time)+sin(2*pi*frequency2*time);
decimatedSignal = decimate(signal,decimationFactor,'fir');
decimatedTime = (0:(signalLength/decimationFactor)-1);
subplot(2,1,1)
stem(time,signal(1:signalLength))
xlabel('Time');
ylabel('Amplitude');
title('Signal')
subplot(2,1,2)
stem(decimatedTime,decimatedSignal(1:signalLength/decimationFactor))
xlabel('Time');
ylabel('Amplitude');
title('Decimated Signal')