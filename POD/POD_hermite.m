clear all; close all; clc

%Comparison of POD vs known soln to harmonic oscillator

%Define the domain:

L = 30; n = 512;
x2 = linspace(-L/2, L/2, n+1);
x = x2(1:n); V = (x.^2).';
t = 0:0.2:10;

%wavenumbers
k  = 2*pi/L * [0:n/2-1, -n/2:-1].';

%Initial data INSANELY IMPORTANT!
%The data u sample really do impact the modes that dominate!
u = exp(-0.2*(x-2).^2);
ut = fft(u);

%Solve in fourier domain
[t,utsol] = ode45('harmonic_rhs', t, ut, [], k, V);

for j = 1:length(t)
    usol(j,:) = ifft(utsol(j,:));
end
%% plot the dynamics
surf(x,t,abs(usol)), shading interp, title('Real dynamics')

%% SVD for dim reduction

X = usol.';
[u,s,v] = svd(X, 'econ');

figure(2), plot(diag(s)/sum(diag(s)), 'ko', 'Linewidth', [2])
%SVD outputs the percentage of energy within each mode

figure(3), plot(abs(u(:,1:3))), legend('Mode 1', 'Mode 2', 'Mode 3')









