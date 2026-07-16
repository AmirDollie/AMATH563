clear all; close all; clc
% input data
x = 1:24;
y = [75, 77, 76,73, 69, 68, 63, 59, 57, 55, 54, 52, 50, 50, 49, 49, 49, 50, 54, 56, 59, 63, 67, 72];
%plot(x,y)

% Genetic alg
m = 200; %generations
n = 50; %trials (guesses)
n2 = 10; % kept trials

%random guesses [50] (with centres of our first estimates):
A = 12 + randn(n,1); 
B = pi/12 + randn(n,1); 
C = 60 + randn(n,1);

for jgen = 1:m
    for j = 1:n
        E(j) = sum((A(j)*cos(B(j)*x) + C(j)-y).^2); %LSE
    end
    plot(E), pause(0.3); %Error convergence

    [Es, Ej] = sort(E); %sort from small to large
    %Es is value of the function smallest to largest
    %Ej are the vector element
    
    %extract best n2 guesses for A,B and C
    Ak1 = A(Ej(1:n2));
    Ak1 = Ak1(:);
    Bk1 = B(Ej(1:n2));
    Bk1 = Bk1(:);
    Ck1 = C(Ej(1:n2));
    Ck1 = Ck1(:);
    
    %Mutate these to get new guesses
    Ak2 = Ak1 + randn(n2,1)/jgen;
    Bk2 = Bk1 + randn(n2,1)/jgen;
    Ck2 = Ck1 + randn(n2,1)/jgen;
    
    Ak3 = Ak1 + randn(n2,1)/jgen;
    Bk3 = Bk1 + randn(n2,1)/jgen;
    Ck3 = Ck1 + randn(n2,1)/jgen;

    Ak4 = Ak1 + randn(n2,1)/jgen;
    Bk4 = Bk1 + randn(n2,1)/jgen;
    Ck4 = Ck1 + randn(n2,1)/jgen;

    Ak5 = Ak1 + randn(n2,1)/jgen;
    Bk5 = Bk1 + randn(n2,1)/jgen;
    Ck5 = Ck1 + randn(n2,1)/jgen;


    A = [Ak1, Ak2, Ak3, Ak4, Ak5];
    B = [Bk1, Bk2, Bk3, Bk4, Bk5];
    C = [Ck1, Ck2, Ck3, Ck4, Ck5];

    %ie 50 best values for A,B,C
end

