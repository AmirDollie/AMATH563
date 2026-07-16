%Script that gives Burgers equation data

clear all; close all; clc

% Burgers' equation
%  u_t + u*u_x - eps*u_xx =0 

% setup
dt=0.1;
t=0:dt:10; 

%diffusion regularization
eps=0.1;

%domain size = 16, split into 256 pts
L=16; n=256; 
x2=linspace(-L/2,L/2,n+1); x=x2(1:n); %257 pts, then remove the last as same as 1st
k=(2*pi/L)*[0:n/2-1 -n/2:-1].'; %wave numbers normalized to size 2pi
k2=fftshift(k); %unshifted domain

% initial data
u=exp(-(x+2).^2).'; %gaussian centred at -2
ut=fft(u); 

%evolve the Burger's equation
[t,utsol]=ode45('burgers_rhs',t,ut,[],k,eps);
  
for j=1:length(t)
    %pull the time snapshots from the Fourier domain to time domain
    usol(:,j)=ifft( utsol(j,1:n).' );
end

%plot the slices of burgers equ
figure(1), 
waterfall(x,t,real(usol.')); colormap([0 0 0]);
view(15,35), set(gca,'Fontsize',[12])
set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])

figure(2), surfl(x,t,real(usol.')); colormap(gray); shading interp
view(15,35)


%saves the data for use later
save('burgers.mat','x','t','usol')








%
% figure(2)
% waterfall(x,t,real(usol.')); colormap([0 0 0]);
% view(15,35), axis off
%set(gca,'Fontsize',[12])
%set(gca,'Fontsize',[12],'Xlim',[-L/2 L/2],'Xtick',[-L/2 0 L/2],'Ylim',[0 10],'Ytick',[0 10],'Zlim',[0 1],'Ztick',[0 0.5 1])