dt=t(2)-t(1);
X=real(usol); %data
[m,n]=size(X); %[num snapshots (101), num spatial points (256)]

%Construct b in Ax = b, recall b here is u_t
Xdot=zeros(m,n-2);
%terrible evil centre differences, recall ignore 1st and last pts
for jj=1:m  % walk through rows (space)
for j=2:n-1  % walk through time
   Xdot(jj,j-1)=( X(jj,j+1)-X(jj,j-1) )/(2*dt); %Xdot =b
end
end

% derv matrices. Differentiation can be written as matrix multiplication
dx=x(2)-x(1);

D=zeros(m,m); %1st derivative
D2=zeros(m,m); %2nd

for j=1:m-1 %structure of 1st derivative matrix. 1 on upper diagonal, -1 on lower diagonal, zeroes elsewhere
  D(j,j+1)=1;
  D(j+1,j)=-1;
%
  D2(j,j+1)=1; %2nd deriv has 1 on off diags and -2 on the main diagonal
  D2(j+1,j)=1;
  D2(j,j)=-2;
end

%periodic boundary conditions
D(m,1)=1;
D(1,m)=-1;
D=(1/(2*dx))*D; %D is the differentiation matrix enforcing periodic boundary conds
%
D2(m,m)=-2;
D2(m,1)=1;
D2(1,m)=1;
D2=D2/(dx^2); %2nd deriv


%library time:
%reshape u such that u is a stack of the different time points, ie time
%goes down the rows
u=reshape( X(:,2:end-1).',(n-2)*m ,1 ); %also take off 1st and last pts


%construct elements of the library
for jj=2:n-1
   ux(jj-1,:)=((D*X(:,jj)).');  % u_x
   uxx(jj-1,:)=((D2*X(:,jj)).');  % u_xx
   u2x(jj-1,:)=((D* (X(:,jj).^2) ).');  % (u^2)_x
end

Ux=reshape(ux,(n-2)*m,1);
Uxx=reshape(uxx,(n-2)*m,1);
U2x=reshape(u2x,(n-2)*m,1);

A=[u u.^2 u.^3 Ux Uxx U2x Ux.*u Ux.*Ux Ux.*Uxx];

Udot=reshape((Xdot.'),(n-2)*m,1);

%2 options to solve, use the QR algorithm inside the '\' to promote
%sparsity:
%xi=A\Udot;

%OR  to explicitly promote sparsity through Lasso (l1) regularization:
xi=lasso(A,Udot,'Lambda',0.002);
bar(xi)




 