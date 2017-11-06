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
C = 1;
tol = 10e-5;

% Initialize alpha constrained to the sum of y*alpha = 0 & alpha >= 0
alpha = zeros(L,1);
for i=1:L
    num = rand(1);
    alpha(i) = num;
end
totalSum = sum(y.*alpha);
if totalSum > 0
    place = find(-1);
    alpha(place) = totalSum + alpha(place);
elseif totalSum < 0
    place = find(1);
    alpha(place) = -totalSum + alpha(place);
end
bias = 0;
count = 0;
while count < 100
    L = size(meas,1);
    % Calculate the weight vector
    w = (alpha.*y).'*x;

    % Calculate KKT Conditions
    KKT = zeros(L,1);
    for i=1:L
        KKT(i) = alpha(i)*(y(i)*(dot(w,x(i,:)) + bias) - 1);
    end
    % Pick x1 and x2
    [~,i1] = max(abs(KKT));
    x1 = x(i1,:);
    %Ei = sum(alpha.*y.*(dot(x(,
    E = zeros(L,1);
    for i=1:L
        for j=1:L
            %E(i) = E(i) + alpha(j)*y(j)*(dot(x(j),x(1)) - dot(x(j),x(i))) + y(i)- y(1)*dot(x(i),x(j));
            tempE(i,j) = alpha(j)*y(j)*dot(x(i),x(j));
        end
    end
    h = sum(tempE');
    E = h-y';
    
    [~,i2] = max(abs(E));
    x2 = x(i2,:);

    if y(i1)~=y(i2)
        L=max(0,alpha(i2)-alpha(i1));
        H=min(C,C+alpha(i2)-alpha(i1));
    else 
        L=max(0,alpha(i1)+alpha(i2)-C);
        H=min(C,alpha(i1)+alpha(i2));
    end


    % Calculate k
    %k  = dot(x1,x1) + dot(x2,x2) - 2*dot(x1,x2);
    k = dot(x(i1),x(i1)) + dot(x(i2),x(i2)) - 2*dot(x(i1),x(i2));
%     if k<=0
%         continue
%     end
    % update alphas
    alpha_old = alpha;
    % not sure which of these to use...one was in class, the other is on HW
    %alpha(i2) = alpha_old(i2) + (y(i2)*E(i2))/k
    alpha(i2) = alpha_old(i2) - y(i2)*(E(i1) - E(i2))/k;
    %alpha(alpha < L)=L;
    %alpha(alpha > H)=H;
    if alpha(i2) > L
        alpha(i2) = L;
    elseif alpha(i2) > H
        alpha(i2) = H;
    end
    alpha(i1) = alpha_old(i1) +y(i1)*y(i2)*(alpha_old(i2) - alpha(i2));


    b1 = bias - E(i1) - y(i1)*(alpha(i1)-alpha_old(i1))*x(i1,:)*x(i1,:)' - y(i2)*(alpha(i2)-alpha_old(i2))*x(i1,:)*x(i2,:)';
    b2 = bias - E(i2) - y(i1)*(alpha(i1)-alpha_old(i1))*x(i1,:)*x(i2,:)' - y(i2)*(alpha(i2)-alpha_old(i2))*x(i2,:)*x(i2,:)';

    if 0<alpha(i1)<C
        bias=b1;
    elseif 0<alpha(i2)<C
        bias=b2;
    else
        bias=(b1+b2)/2;
    end
    
%     if(alpha(i2) == alpha2_old && alpha(i1) == alpha1_old)
%         break;
%     end
    
    fprintf('count, %i\n', count)
    count=count+1;
end

positives = x(1:50,:);
negatives = x(51:150,:);
pos = 0;
neg = 0;
for i=1:50
    pos = pos + sum(positives(i).*w);
end
for i = 1:100
    neg = neg + sum(negatives(i).*w);
end
bias = -(pos/50+neg/100)/2;

countDelicious = 0;
countNasty = 0;
for k = 1:150
    result = bias + sum(x(k).*w);
    if(result >0 && y(k) > 0)
        countDelicious = countDelicious + 1;
    elseif(result<0 && y(k) <0)
        countDelicious = countDelicious + 1;
    else
        countNasty = countNasty + 1;
    end
end



total= sum(y.*alpha)
fprintf('Accuracy of %0.1i percent\n', (countDelicious/150)*100);
