% Doppler Velocity Calculation
c = 3*10^8;         %speed of light
frequency = 77e9;   %frequency in Hz

% TODO: Calculate the wavelength
lamda = c/frequency

% TODO: Define the doppler shifts in Hz using the information from above 
fd = [3e3, -4.5e3, 11e3, -3e3];

% TODO: Calculate the velocity of the targets  fd = 2*vr/lambda
% frequency shift due to velocity (fd)
% frequecy shift due to range (fr)
object_v = (fd*lamda)/2
% IMPORTANT: this is not relative velocity but absolute velocity of the
% target


% TODO: Display results