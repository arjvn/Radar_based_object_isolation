% Implement21D CFAR using lagging cells on the given noise and target scenario.

% Close and delete all currently open figures
close all;

% Specify the parameters of a signal with a sampling frequency of 1 kHz 
% and a signal duration of 1.5 seconds.

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

% Form a signal containing a 50 Hz sinusoid of amplitude 0.7 and a 120 Hz 
% sinusoid of amplitude 1.

Signal = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

% Corrupt the signal with zero-mean white noise with a variance of 4
Noisy_signal = Signal + 2*randn(size(t));
Noisy_signal = abs(Noisy_signal);

%Targets location. Assigning bin 100, 200, 300, and 700 as Targets
%  with the amplitudes of 16, 18, 27, 22.
Noisy_signal([100 ,200, 300, 700])=[16 18 27 22];


% Compute the 2 sided spectrum P2. Then compute the single-sided spectrum
% P1 based on P2 and the even-valued signal length L.

FFT_signal = fft(Noisy_signal);
FFT_signal = abs(FFT_signal/L);
% FFT_signal = FFT_signal(1:L/2+1);



% Plot the noisy signal in the time domain. It is difficult to identify
% the frequency components by looking at the signal X(t).
figure(1);
tiledlayout(1,2)

% left plot
nexttile
plot(FFT_signal)
title('Signal corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')





% Determine the number of Training cells for each dimension Tr and Td. 
% Similarly, pick the number of guard cells Gr and Gd
% Apply CFAR to detect the targets by filtering the noise.
% TODO: Define the number of Training Cells
Tr = 20;
Td = 10;
% TODO: Define the number of Guard Cells 
Gr = 5;
Gd = 3;
% TODO: Define Offset (Adding room above noise threshold for the desired SNR)
offset = 5;
% Define CUT cell
CUT = 1;
