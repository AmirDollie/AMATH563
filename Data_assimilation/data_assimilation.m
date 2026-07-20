clear all; close all; clc

t=0:0.01:20;
sigma=10; b=8/3; r=28;

x0=[5 5 5];
[t,xsol]=ode45('lor_rhs',t,x0,[],sigma,b,r);

x_true=xsol(:,1);
y_true=xsol(:,2);
z_true=xsol(:,3);
%figure(1), plot3(x_true,y_true,z_true,'Linewidth',[2]);
%grid on

sigma2=1; % strength of noise
for j=1:8
   xic=x0+sigma2*randn(1,3); 
[t,xsol]=ode45('lor_rhs',t,xic,[],sigma,b,r);
subplot(4,2,j)
plot(t,x_true,'m',t,xsol(:,1),'g','Linewidth',[2])

end

close all
% noisy observations
tdata=t(1:50:end);
n=length(tdata);
xn=randn(n,1); yn=randn(n,1); zn=randn(n,1);
sigma3=1; % error strength
xdata=x_true(1:50:end)+sigma3*xn;
ydata=y_true(1:50:end)+sigma3*yn;
zdata=z_true(1:50:end)+sigma3*zn;

plot(tdata,xdata,'ro',t,x_true,'k','Linewidth',[2])
close all
x_da=[];
for j=1:length(tdata)-1
    tspan=0:0.01:0.5;
    [tspan,xsol]=ode45('lor_rhs',tspan,xic,[],sigma,b,r);
    xic0=[xsol(end,1); xsol(end,2); xsol(end,3)];
    xdat=[xdata(j+1); ydata(j+1); zdata(j+1)];
    K=sigma2/(sigma2+sigma3);
    xic=xic0+(K*(xdat-xic0));
    
    x_da=[x_da; xsol(1:end-1,:)];
    
end
x_da=[x_da; xsol(end,:)];
plot(t,x_true,'m',t,x_da(:,1),'g','Linewidth',[2])
legend('True Data', 'Assimilated Estimate')
















































