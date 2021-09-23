%Operating frequency (Hz)
fc = 77.0e9;

%Transmitted power (W)
Ps = 3e-3;

%Antenna Gain (linear) (dBi)
G =  10000;

%Minimum Detectable Power
Pe = 1e-10;

%RCS of a car (m^2)
RCS = 100;

%Speed of light (m/s^2)
c = 3*10^8;

%Range resolution

%Radar maximum range (m)
R_max = 300;

% TASK 1
%TODO: Calculate the wavelength
lamda = c/fc;

%TODO : Measure the Maximum Range a Radar can see. 
% Radar range equation
R = nthroot(((Ps*G^2*lamda^2*RCS)/Pe*(4*pi)^3),4);

% TASK 2
% TODO : Find the B_sweep of chirp for 1 m resolution (d_res)
d_res = 1;
B_sweep = (c)/(2*d_res);

% TODO : Calculate the chirp time based on the Radar's Max Range
% The sweep time (T_s) and the sweep bandwidth (B_sweep) are generated by
% the Sythesizer.
% As a rule of thumb, sweep time should be at least 5 to 6 times the round
% trip time (maximum distance)
T_s = 5.5 * (2 * R_max)/c;

% TODO : define the frequency shifts 
f_b = [0, 1.1e6, 13e6, 24e6];
R = (c*T_s*f_b)/(2*B_sweep)
