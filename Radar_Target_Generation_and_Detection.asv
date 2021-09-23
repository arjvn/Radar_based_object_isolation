clc
clear
clear all
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 77e9; % Radar_frequency;
Radar_Range = 200;
Radar_range_resolution = 1;
Radar_max_Velocity = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e8;
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
Range_input = 110;
velocity_input = 70;
 
%% FMCW Waveform Generation
% *%TODO* :
% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

B = c / (2*Radar_range_resolution);
% The sweep time can be computed based on the time needed for the signal to travel 
% the unambiguous maximum range. In general, for an FMCW radar system, the sweep time 
% should be at least 5 to 6 times the round trip time. This example uses a factor of 5.5.
% https://ww2.mathworks.cn/help/radar/ug/automotive-adaptive-cruise-control-using-fmcw-technology.html
TChirp = 5.5 * (2*Radar_Range/c);
slope = B / TChirp;
                                                         
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*TChirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
range_target = zeros(1,length(t));
td = zeros(1,length(t)); % radar RX signal time delay


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
              
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    range_target(i) = Range_input+(velocity_input*t(i));
    td(i) = 2*range_target(i) / c;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+(slope*(t(i))^2)/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i))+(slope*(t(i)-td(i))^2)/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix = Tx.*Rx;
    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.

Mix = reshape(Mix,[Nr, Nd]);

% *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.

Mix_fft1 = fft(Mix, Nr); % treats each column as a vector

 % *%TODO* :
% Take the absolute value of FFT output
Mix_fft1 = abs(Mix_fft1);
Mix_fft1= Mix_fft1./max(Mix_fft1); %normalize max value of each column to 1
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
Mix_fft1 = Mix_fft1(1:Nr/2);

%plotting the range
figure ('Name','Range from First FFT')

 % *%TODO* :
 % plot FFT output 
plot(Mix_fft1)
axis ([0 200 0 1]);
title('Normalized range from First FFT')
xlabel('Range');
ylabel('Normalized Signal strength')

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd); %speed
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400); %range
figure('Name', 'FFT2 plot');
surf(doppler_axis,range_axis,RDM); % surf -> creates a 3D plot
xlabel('Speed');
ylabel('Range');
zlabel('Signal Strength');

%% CFAR implementation
%Slide Window through the complete Range Doppler Map
% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10; %range 
Td = 8; %doppler

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 6;
Gd = 4;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 1.4;

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

% move the kernel across the RDM
RDM = RDM/max(max(RDM)); % Normalizing
for i=Tr+Gr+1:(Nr/2)-Tr-Gr
    for j = Td+Gd+1:(Nd)-Td-Gd
        %Create a vector to store noise_level for each iteration on training cells
        noise_level = zeros(1,1);
        for x = i-Tr-Gr:i+Tr+Gr
            for y = j-Td-Gd:j+Td+Gd
                if (abs(i-x)>Gr || abs(j-y)>Gd)
                    noise_level = noise_level+db2pow(RDM(x,y));
                end
            end
        end
        no_training_cells = (2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1);
        threshold = pow2db(noise_level/no_training_cells); % avg noise lvl
        threshold = threshold + offset;
        CUT = RDM(i,j);
        if (CUT < threshold)
            RDM(i,j) = 0;
        else
            RDM(i,j) = 1;
        end
    end
end


% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
RDM(RDM~=0 & RDM~=1) = 0;

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name', 'CA-CFAR Filtered RDM')
surf(doppler_axis,range_axis,RDM);
colorbar;
title('CA-CFAR Filtered Radar output');
xlabel('Speed');
ylabel('Range');
zlabel('Signal Strength')

figure('Name', 'CA-CFAR Filtered RDM: Signal Strength Vs Speed')
surf(doppler_axis,range_axis,RDM);
colorbar;
title('CA-CFAR Filtered Radar output: Signal Strength Vs speed');
xlabel('Speed');
ylabel('Range');
zlabel('Signal Strength')
view(0,0)

figure('Name', 'CA-CFAR Filtered RDM: Signal Strength Vs Range')
surf(doppler_axis,range_axis,RDM);
colorbar;
title('CA-CFAR Filtered Radar output: Signal Strength Vs Range');
xlabel('Speed');
ylabel('Range');
zlabel('Signal Strength')
view(90,0)