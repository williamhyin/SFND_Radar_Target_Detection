% Doppler Velocity Calculation
c = 3*10^8;         %speed of light
frequency = 77e9;   %frequency in Hz

% TODO : Calculate the wavelength

wavelength = c/frequency;

% TODO : Define the doppler shifts in Hz using the information from above 

doppler_shift= [3e3 -4.5e3  11e3 -3e3];

% TODO : Calculate the velocity of the targets  fd = 2*vr/lambda

vr=doppler_shift*wavelength/2;

% TODO: Display results

disp(vr)