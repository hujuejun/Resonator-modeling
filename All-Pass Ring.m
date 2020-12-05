% Calculate transmission spectrum of a ring resonator in an all-pass
% configuration
% Reference: A. Yariv, Electron. Lett. 36, 321-322 (2000)
% By JJ Hu @ MIT

clear all;

psi = 0 * pi; % phase delay of the ring coupler
k = 0.1032*exp(i*psi);    % amplitude coupling coefficient of the ring coupler
t = sqrt(1-(abs(k)).^2)*exp(i*psi);

lamda = 1494.8:0.0001:1495; % wavelength in nm
L_ring = 40*pi;    % length of the ring (except the coupler part) in micron
n_eff = 2.1175-(lamda-min(lamda))*0.0001;    % effective index of waveguide

% load(['loss_ring.txt'],'loss_ring','-ascii');   % load wavelength depedent loss figure saved in the file loss_ring.txt
% loss_ring = transpose(loss_ring);
loss_ring = 3.7;  % Equivalent waveguide loss in the ring in dB/cm

alpha = sqrt(10.^(-loss_ring/1E5*L_ring));   % round trip amplitude loss in the ring
theta = L_ring.*n_eff./(lamda/1000)*2*pi;    % round trip phase delay in the ring

insertion = 0;   % insertion loss

b1 = (-alpha+t.*exp(-i*theta))./(-alpha.*conj(t)+exp(-i*theta));
hold on;
plot(lamda, log10(abs(b1).^2)*10-insertion,'r')     % plot transmitted intensity in dB scale
% plot(lamda, abs(b1).^2,'r')     % plot transmitted intensity in linear scale
figure;
plot(lamda, angle(b1),'b')  % plot phase of the transmitted light

DataOutput = transpose([lamda; abs(b1).^2]);
% save(['Transmission.txt'],'DataOutput','-ascii');