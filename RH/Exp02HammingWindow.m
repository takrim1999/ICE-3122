signal = load('Data/plain');
% x(10:15);
startSample = 12001;
interval = 512;
stopSample = startSample+interval-1;
sample = signal(startSample:stopSample);
for i=(1:interval)
hammingWindow(i) = 0.54-0.46*cos(2*pi*(i-1)/(interval-1));
end
windowedSample = sample.*hammingWindow';
% length(windowedSample);
% windowedSample(1:5);
subplot(3,1,1)
plot(x)
title('Voice Signal')
xlabel('Time')
ylabel('Amplitude')
subplot(3,1,2)
plot(sample)
title('Sample')
xlabel('Time')
ylabel('Amplitude')
subplot(3,1,3)
plot(windowedSample)
title('Hamming Window')
xlabel('Time')
ylabel('Amplitude')
