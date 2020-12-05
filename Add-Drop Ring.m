% This code uses a matrix representation to simulate the transfer function
% of an add-drop ring resonator; note that optical power can be launched 
% into both the input and add ports simultaneously.
% Author: JJ Hu @ MIT

clear all;
ii = 1;

lambda_step = 0.0001;   % Wavelength step in nm
lambda = 1484:lambda_step:1485;  % Wavelength range in nm

% Launched wave complex amplitudes (wavelength independent)
s_i = 1;    % Input port
s_a = 1;    % Add port

for lamda = lambda; % wavelength in nm

% Basic ring parameters
L_ring = 400*pi;    % Total length of the ring in um
loss_ring = 0.2;  % Equivalent waveguide loss in the ring in dB/cm
n_eff = 2.1175;    % Effective index of waveguide

% Coupler parameters
k1 = 0.3;    % Amplitude coupling coefficient of the ring coupler #1
t1 = sqrt(1-(abs(k1)).^2);  % Note that a real value is taken
k2 = 0.3;    % Amplitude coupling coefficient of the ring coupler #2
t2 = sqrt(1-(abs(k2)).^2);
L1 = 200*pi;    % Length of ring section 1 between the two couplers in um
L2 = L_ring - L1;   % Length of ring section 2 between the two couplers in um
alpha1 = sqrt(10.^(-loss_ring/1E5*L1));   % amplitude loss coefficient in the ring section 1
alpha2 = sqrt(10.^(-loss_ring/1E5*L2));   % amplitude loss coefficient in the ring section 2
psi1 = L1.*n_eff./(lamda/1000)*2*pi;    % Phase delay in the ring section L1
phi1 = exp(-i*psi1) * alpha1;
psi2 = L2.*n_eff./(lamda/1000)*2*pi;    % Phase delay in the ring section L2
phi2 = exp(-i*psi2) * alpha2;

A = [t1         -1          0           0           0           0;
    conj(k1)    0           0           0           1           0;
    0           0           t2          -1          0           0;
    0           0           conj(k2)    0           0           1;
    0           phi2        -1          0           0           0;
    -1          0           0           phi1        0           0];

B = [-k1*s_i;   conj(t1)*s_i;   -k2*s_a;   conj(t2)*s_a;    0;  0];

a = inv(A)*B;

s_t(ii) = a(5); % Through port complex amplitude
s_d(ii) = a(6); % Drop port complex amplitude
ii = ii + 1;

end

hold on
plot(lambda, abs(s_t).^2, 'b')
plot(lambda, abs(s_d).^2, 'r')