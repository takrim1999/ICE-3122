signal = load('Data/plain');
startSegment = 17001;
interval = 512;
endSegment = startSegment+interval-1;
segment = signal(startSegment:endSegment);
fftSegment = abs(fft(segment,interval).^2);
ifftSegment = real(ifft(fftSegment));
subplot(2,2,1)
plot(signal)
title('Voice Signal')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,3)
plot(segment)
title('Segment')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,2)
plot(fftSegment)
title('FFT Segment')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,4)
plot(ifftSegment)
title('IFFT Segment')
xlabel('Time')
ylabel('Amplitude')