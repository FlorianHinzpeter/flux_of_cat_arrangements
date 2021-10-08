function Dist = SpatialDist(N,n,arrangement_path)

load(arrangement_path);

a = -1:0.25:2.75;
alpha = 10.^a;
DistE1E2 = zeros(2,length(alpha));
Dist = zeros(2,length(alpha));
RadE1 = zeros(2,length(alpha));
RadE2 = zeros(2,length(alpha));
ClosestDist = zeros(2,length(alpha));
cv = zeros(1,length(alpha));


X1Max = zeros(2*n,N/2);
X2Max = zeros(2*n,N/2);

X1Min = zeros(2*n,N/2);
X2Min = zeros(2*n,N/2);

XMax = zeros(2*n,N);
XMin = zeros(2*n,N);

for idx1=1:length(alpha)

load(sprintf('N=%i_alpha_%i.mat',N,alpha(idx1)));

cv(idx1)=std(J(1:end-1,2))/mean(J(1:end-1,2));

[Mmax,imax] = sort(J(1:end-1,2),'descend');

[Mmin,imin] = sort(J(1:end-1,2),'ascend');

for idx11 = 1:n
    
X1Max(2*idx11-1:2*idx11,:) = X(2*imax(idx11)-1:2*imax(idx11),1:2:N);

X2Max(2*idx11-1:2*idx11,:) = X(2*imax(idx11)-1:2*imax(idx11),2:2:N);

XMax(2*idx11-1:2*idx11,:) = X(2*imax(idx11)-1:2*imax(idx11),:);

X1Min(2*idx11-1:2*idx11,:) = X(2*imin(idx11)-1:2*imin(idx11),1:2:N);

X2Min(2*idx11-1:2*idx11,:) = X(2*imin(idx11)-1:2*imin(idx11),2:2:N);

XMin(2*idx11-1:2*idx11,:) = X(2*imin(idx11)-1:2*imin(idx11),:);

end

% average distance btw. E1 and E2

KumDistMaxE1E2 = 0;
KumDistMinE1E2 = 0;
KumDistMax = 0;
KumDistMin = 0;
KumRadE1Max = 0;
KumRadE2Max = 0;
KumRadE1Min = 0;
KumRadE2Min = 0;
KumClosestDistMax = 0;
KumClosestDistMin = 0;
DistListMax = zeros(1,N/2);
DistListMin = zeros(1,N/2);

count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;

for idx2 = 1:n

for idx3 = 1:N/2    
    
    count4 = count4+1;
    
for idx4 = 1:N/2
    
     
   count1 = count1+1;
       
   KumDistMaxE1E2 = KumDistMaxE1E2 + sqrt((X1Max(2*idx2-1,idx3)-X2Max(2*idx2-1,idx4))^2+(X1Max(2*idx2,idx3)-X2Max(2*idx2,idx4))^2);
   KumDistMinE1E2 = KumDistMinE1E2 + sqrt((X1Min(2*idx2-1,idx3)-X2Min(2*idx2-1,idx4))^2+(X1Min(2*idx2,idx3)-X2Min(2*idx2,idx4))^2);  
   
   DistListMax(idx4) = sqrt((X1Max(2*idx2-1,idx3)-X2Max(2*idx2-1,idx4))^2+(X1Max(2*idx2,idx3)-X2Max(2*idx2,idx4))^2);
   
   DistListMin(idx4) = sqrt((X1Min(2*idx2-1,idx3)-X2Min(2*idx2-1,idx4))^2+(X1Min(2*idx2,idx3)-X2Min(2*idx2,idx4))^2);
 
end

KumClosestDistMax = KumClosestDistMax + min(DistListMax);
KumClosestDistMin = KumClosestDistMin + min(DistListMin);

end

for idx5 = 1:N-1
    
    for idx6 = idx5+1:N
      
        count2 = count2+1;
        KumDistMax = KumDistMax +sqrt((XMax(2*idx2-1,idx5)-XMax(2*idx2-1,idx6))^2+(XMax(2*idx2,idx5)-XMax(2*idx2,idx6))^2);
        KumDistMin = KumDistMin +sqrt((XMin(2*idx2-1,idx5)-XMin(2*idx2-1,idx6))^2+(XMin(2*idx2,idx5)-XMin(2*idx2,idx6))^2);
        
        
    end
    
end

for idx7 = 1:N/2
   
    count3 = count3+1;
    
     KumRadE1Max = KumRadE1Max + sqrt(X1Max(2*idx2-1,idx7)^2+X1Max(2*idx2,idx7)^2);
     KumRadE1Min = KumRadE1Min + sqrt(X1Min(2*idx2-1,idx7)^2+X1Min(2*idx2,idx7)^2);
     
     KumRadE2Max = KumRadE2Max + sqrt(X2Max(2*idx2-1,idx7)^2+X2Max(2*idx2,idx7)^2);
     KumRadE2Min = KumRadE2Min + sqrt(X2Min(2*idx2-1,idx7)^2+X2Min(2*idx2,idx7)^2);
end

end

DistE1E2(1,idx1) = KumDistMaxE1E2/count1;
DistE1E2(2,idx1) = KumDistMinE1E2/count1;

Dist(1,idx1) = KumDistMax/count2;
Dist(2,idx1) = KumDistMin/count2;

RadE1(1,idx1) = KumRadE1Max/count3;
RadE1(2,idx1) = KumRadE1Min/count3;

RadE2(1,idx1) = KumRadE2Max/count3;
RadE2(2,idx1) = KumRadE2Min/count3;

ClosestDist(1,idx1) = KumClosestDistMax/count4;
ClosestDist(2,idx1) = KumClosestDistMin/count4;




end

figure
subplot(4,1,1)     
semilogx(alpha,ClosestDist(1,:),'r',alpha,ClosestDist(2,:),'b')
title('Closest Distance btw E_1 and E_2')

subplot(4,1,2)       
semilogx(alpha,RadE1(1,:),'-r',alpha,RadE1(2,:),'-b',alpha,RadE2(1,:),'--r',alpha,RadE2(2,:),'--b')     
title('Radial position')

subplot(4,1,3)
semilogx(alpha,Dist(1,:),'-r',alpha,DistE1E2(1,:),'--r',alpha,Dist(2,:),'-b',alpha,DistE1E2(2,:),'--b')
title('Overall Distance between all E (solid) between E1 and E2 (dashed)')

subplot(4,1,4)     
semilogx(alpha,cv)
title('Coefficient of Variation')

end




