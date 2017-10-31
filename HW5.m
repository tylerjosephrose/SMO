% SMO algorithm for SVM
% Homework Assignment 5
% Seth Dippold and Tyler Rose

clear all; close all;

% Read in the data here
x = randi(15,10,3);%training data x values
y = mod(randi(2,10,1),2);%training data y values
L = length(y);
C = 0.5;
tol = 10e-5;

% Initialize alpha constrained to the sum of y*alpha = 0
alpha = zeros(L,1);
for i=1:L
    num = rand(1)-.5;
    alpha(i) = num;
end
replacement = find(y,1);
totalSum = sum(y.*alpha);
alpha(replacement) = alpha(replacement) + -totalSum;

% Initialize b to 0
b = 0;

% Calculate the weight vector
w = (alpha.*y).'*x;

% Calculate KKT Conditions
KKT = zeros(L,1);
for i=1:L
    KKT(i) = alpha(i)*(y(i)*(dot(w,x(i,:) + b) - 1));
end
% Pick x1 and x2
[~,i1] = max(KKT);
x1 = x(i1,:);
%Ei = sum(alpha.*y.*(dot(x(,
E = zeros(L,1);
for i=1:L
    for j=1:L
        E(i) = E(i) + alpha(j)*y(j)*(dot(x(j),x(1)) - dot(x(j),x(i))) + y(i)- y(i)*dot(x(i),x(j));
    end
end
[~,i2] = max(E);
x2 = x(i2,:);

% Calculate k
k  = dot(x1,x1) + dot(x2,x2) - 2*dot(x1,x2);
