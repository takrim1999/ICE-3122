signalLength = 100;
time = 0:signalLength-1;
frequency1=0.02;
frequency2=0.03;
interFactor = 3;
deciFactor = 5;
signal = sin(2*pi*frequency1*time)+sin(2*pi*frequency2*time);
resampledSignal = resample(signal,interFactor,deciFactor);
length(resampledSignal)
resampledTime = 0:((signalLength*interFactor)/deciFactor)-1
subplot(2,1,1)
stem(time,signal);
subplot(2,1,2)
stem(resampledTime,resampledSignal)