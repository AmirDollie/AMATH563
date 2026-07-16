clear all; close all; clc

t = 0:0.01:100;
x0 = [5,5,5];
sigma = 10; b = 8/3; r = 28;
[t, Y] = ode45('lor_rhs', t, x0, [], sigma, b, r);
%solve Lorentz
x = Y(:,1); y = Y(:,2); z = Y(:,3);

plot3(x,y,z, 'Linewidth', [2]), grid on

%Here we have x, y and z. But what if we only had x?
%How do we even know y and z exist?
%In other words, what is the dimension of the data??

%Let's see what each variable is doing:
figure(2)
subplot(3,1,1), plot(t,x, 'Linewidth', [2])
subplot(3,1,2), plot(t,y, 'Linewidth', [2])
subplot(3,1,3), plot(t,z, 'Linewidth', [2])

%Let's now somehow use only x time series data, to infer stuff about y and
%z.
%% Time delaying embedding   

%construct the Hankel matrix:
H = [x(1:5000).'
     x(2:5001).'
     x(3:5002).'
     x(4:5003).'
     x(5:5004).'
     x(6:5005).'
     x(7:5006).'
     x(8:5007).'];

%Now let's take the SVD of this matrix:

[u,s,v] = svd(H,'econ');
%Plot the singular values
figure(3), plot(diag(s)/sum(diag(s)), 'ro', 'LineWidth', [3]), title("Singular values of H")

 %From figure u can see the first ~3 modes dominate and thus there are approximately 3 DOF/terms.dynamical rank.

 %Now let's see what the structures look like:
%EIGENSPACE:
 figure(4), plot(u(:,1:5), 'LineWidth',[2])
 legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5')

 %Note that these modes are polynomials!
 % Thus the time delaying of the single variable produces a polynomial
 % basis set

 %Now for the time series: let's plot the 3 dominant time series ie what
 %are the time dynamics like? (gotten from the right singular matrix)

figure(5), plot3(v(:,1), v(:,2), v(:,3), 'Linewidth', [2])
grid on;

%Warped representation of actual dynamics!!!!!
%It inherits its structure from the full dynamics

%Gives a "shadow" of the true manifold


%The 4th mode (which we truncated) can be seen to SPIKE when the x switches
%lobes!!
figure(6), plot(300*v(:,4), 'Linewidth', [2])
hold on;
plot(x(1:5000))

%The spikes actually happen BEFORE the swap. 
%It acts as a CONTROL signal for the 3x3 embedding!!!




