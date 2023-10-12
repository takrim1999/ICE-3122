signal = load('Data/plain');
startSegment = 12001;
interval = 512; 
endSegment = startSegment+interval-1;
segment = signal(startSegment:endSegment);
for i=1:interval
    hanningWindow(i) = 0.5*cos(1-(2*pi*(i-1))/(interval-1));
end
for i=(1:interval)
    hammingWindow(i) = 0.54-0.46*cos(2*pi*(i-1)/(interval-1));
end
hanningSample = sample.*hanningWindow';
hammingSample = sample.*hammingWindow';
subplot(2,2,1)
plot(x)
title('Voice Signal')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,3)
plot(sample)
title('Sample')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,2)
plot(hanningSample)
title('Hanning Window')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,4)
plot(hammingSample)
title('Hamming Window')
xlabel('Time')
ylabel('Amplitude')