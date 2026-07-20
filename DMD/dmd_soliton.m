clear all; close all; clc

%spatial points
xi = linspace(-10, 10, 400);

%time points
t = linspace(0, 4*pi, 200);
dt = t(2) - t(1);

%Make a grid:
[Xgrid, T] = meshgrid(xi,t);

%create two spatio-temporal patterns and mix them. Goal will be how to pull
%them apart

f1 = sech(Xgrid + 3).*(1*exp(i*2.3*T));
f2 = (sech(Xgrid).*tanh(Xgrid)).*(2*exp(i*2.8*T));
f = f1 + f2;

%Here we are sort of putting it into the form:
% x = f = sum bj psij exp(iwjt) with wj = [2.3, 2.8]

figure;
subplot(2,2,1), surfl(Xgrid,T,real(f1)); shading interp, colormap(gray), title('Real(f1)')
subplot(2,2,2), surfl(Xgrid,T,real(f2)); shading interp, colormap(gray), title('Real(f2)')
subplot(2,2,3), surfl(Xgrid,T,real(f)); shading interp, colormap(gray), title('Real(f)')

%How do we get f1 and f2 from f?

%SVD time:
[u,s,v] = svd(f.');
figure(2)
plot(diag(s)/(sum(diag(s))), 'ro')

%can see that almost all of the variance is taken by the 1st 2 modes (as
%expected)

figure(3) %Let's see the space and time dynamics:
subplot(2,1,1), plot(real(u(:,1:3))), legend('Mode 1', 'Mode 2', 'Mode 3') %Spatial
subplot(2,1,2), plot(real(v(:,1:3))), legend('Mode 1', 'Mode 2', 'Mode 3') %Temporal

%Clearly the 3rd mode is USELESS so u can truncate it

figure(4) %Let's see the space and time dynamics:
subplot(2,1,1), plot(real(u(:,1:2))), legend('Mode 1', 'Mode 2') %Spatial
subplot(2,1,2), plot(real(v(:,1:2))), legend('Mode 1', 'Mode 2') %Temporal

%Here we can see that using just the SVD, the 2 modes do not recover the
%actual functions used to create f (namely sech and tanh), because it does
%not enforce linearity

%% DMD 

X = f.'; %rows = space, col = time as usual for DMD
X1 = X(:,1:end-1);
X2 = X(:,2:end);

%STEP 1: Rank truncation
r = 2;
[U,S,V] = svd(X1,'econ');
Ur = U(:,1:r);
Sr = S(1:r, 1:r);
Vr = V(:,1:r);

%STEP 2: Atilde

Atilde = Ur' * X2 * Vr / Sr;

[W,D] = eig(Atilde);
Phi = X2*Vr/Sr*W; %DMD modes

% STEP 3: Extract DMD eigenvalues and compute temporal dynamics
lambda = diag(D);
omega = log(lambda)/dt; 

%Let's look at these modes overlaid on the previous figure
figure(5) %Let's see the space and time dynamics:
subplot(2,1,1), plot(real(u(:,1:2))) %Spatial
%Now for DMD overlaid
subplot(2,1,1), hold on, plot(real(Phi)), legend('SVD Mode 1', 'SVD Mode 2', 'DMD Mode 1', 'DMD Mode 2')
%Clearly recovers the precisely right modes!!!




subplot(2,1,2), plot(real(v(:,1:2))), legend('Mode 1', 'Mode 2') %Temporal

%% Time reconstruction 
%Need to check if f is recovered given the 2 modes

%IC:
x1 = X(:,1);
b = Phi\x1; %See notes
time_dynamics = zeros(r,length(t));

%Build time dynamics
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter))); %outer product
end

X_dmd = Phi * time_dynamics;


figure(6)
subplot(1,2,1), surfl(Xgrid,T,real(f)); shading interp, colormap(gray), title('Real(f)')
subplot(1,2,2), surfl(Xgrid,T,real(X_dmd).'); shading interp, colormap(gray), title('Real(DMD Reconstructed)')

%Perfect! 
%You can use this for forecasting, allowing t to exceed the original set ie
%say t = 0 -> 20 or whatever.


%But again, due to the exponential nature of this, if Re{eig} is large,
%then you will blow up at large t (unstable system)








