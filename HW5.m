% SMO algorithm for SVM
% Homework Assignment 5
% Seth Dippold and Tyler Rose

clear all; close all;

% Read in the data here
%x = randi(15,10,3);%training data x values
%y = mod(randi(2,10,1),2);%training data y values
load fisheriris;
x = meas;
L = size(meas,1);
y = zeros(L,1);
for i=1:L
   if strcmp(species(i),'setosa');
       y(i) = 1;
   else
       y(i) = -1;
   end
end
C = 0.5;
tol = 10e-5;

% Initialize alpha constrained to the sum of y*alpha = 0 & alpha >= 0
alpha = zeros(L,1);
% for i=1:L-1
%     num = rand(1);
%     alpha(i) = num;
% end
% totalSum = sum(y.*alpha);
% if y(L) == 1
%     alpha(L) = -totalSum;
% else
%     alpha(L) = totalSum;
% end

% Initialize b to 0
b = 0;

% loop
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

if y(i1)~=y(i2)
    L=max(0,alpha(i2)-alpha(i1));
    H=min(C,C+alpha(i2)-alpha(i1));
else 
    L=max(0,alpha(i1)+alpha(i2)-C);
    H=min(C,alpha(i1)+alpha(i2));
end


% Calculate k
k  = dot(x1,x1) + dot(x2,x2) - 2*dot(x1,x2);
% if k<=0
%     continue
% end
% update alphas
alpha1_old = alpha(i1);
alpha2_old = alpha(i2);
alpha(i2) = alpha2_old + (y(i2)*E(i2))/k
alpha(i1) = alpha1_old +y(i1)*y(i2)*(alpha2_old - alpha(i2))

alpha(alpha < L)=L;
alpha(alpha > H)=H;

b1 = b - E(i1) - y(i1)*(alpha(i1)-alpha1_old)*x(i1,:)*x(i1,:)' - y(i2)*(alpha(i2)-alpha2_old)*x(i1,:)*x(i2,:)';
b2 = b - E(i2) - y(i1)*(alpha(i1)-alpha1_old)*x(i1,:)*x(i2,:)' - y(i2)*(alpha(i2)-alpha2_old)*x(i2,:)*x(i2,:)';

if 0<alpha(i1)<C
    bias=b1;
elseif 0<alpha(i2)<C
    bias=b2;
else
    bias=(b1+b2)/2;
end

