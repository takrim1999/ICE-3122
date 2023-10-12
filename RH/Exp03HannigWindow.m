signal = load('Data/plain');
startSegment = 12001;
interval = 512; 
endSegment = startSegment+interval-1;
segment = signal(startSegment:endSegment);
for i=1:interval
    hanningWindow(i) = 0.5*cos(1-(2*pi*(i-1))/(interval-1));
end
% hannig
% length(hanningWindow)
hanningSignal = segment.*hanningWindow';
subplot(3,1,1)
plot(signal);
title('signal');
subplot(3,1,2)
plot(segment);
title('segment');
subplot(3,1,3)
plot(hanningSignal);
title('segment');
