function CV

a=-1:0.25:2.75;
alpha = [10.^a];

cv = zeros(size(alpha));
Mean = zeros(size(alpha));
Std = zeros(size(alpha));
Skewness = zeros(size(alpha));
Kurtosis = zeros(size(alpha));
StdPos = zeros(size(alpha));
StdNeg = zeros(size(alpha));
Max = zeros(size(alpha));
Min = zeros(size(alpha));

for idx1 = 1:length(alpha)

load(sprintf('N=60_alpha_%i.mat',alpha(idx1)));

%CV:
cv(idx1)=std(J(1:end-1,2))/mean(J(1:end-1,2));
Skewness(idx1) = skewness(J(1:end-1,2));
Kurtosis(idx1) = kurtosis(J(1:end-1,2)); 

%mean of flux
Mean(idx1) = mean(J(1:end-1,2));


%standard deviation
Std(idx1) = std(J(1:end-1,2));

Max(idx1) = max(J(1:end-1,2));
Min(idx1) = min(J(1:end-1,2));


n1 = 0;
n2 = 0;
x = 0;
y = 0;

for idx2 = 1:length(J(:,1))
   
    if J(idx2,2)>Mean(idx1)
        
       n1 = n1+1; 
       x = x + (J(idx2,2)-Mean(idx1))^2; 
        
    else
        
       n2 = n2+1; 
       y = y + (J(idx2,2)-Mean(idx1))^2; 
    
    end
    
end

StdPos(idx1) = sqrt(x/n1);
StdNeg(idx1) = sqrt(y/n2);


%alpha(idx)
%num2str(cv(idx))

end

CVPos = StdPos./Mean;
CVNeg = StdNeg./Mean;

loglog(alpha,Max,alpha,Min)

figure
subplot(3,1,1)     
semilogx(alpha,cv)
title('CV')

subplot(3,1,2)       
semilogx(alpha,Skewness)      
title('Skewness')

subplot(3,1,3)       
semilogx(alpha,Kurtosis)      
title('Kurtosis')

%subplot(3,1,3)       
%semilogx(alpha,cv)      
%title('Coeff. of Variation')

%xlabel('\alpha')
