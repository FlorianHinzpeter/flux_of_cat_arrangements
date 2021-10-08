function x = Randomarrangement_clusters(r0,rE,N,mu)
%---------------------------------------------------------------------------------------------

%This function generates a uniform arrangement of catalysts in a cluster with radius r0.
%The cluster center is randomly positioned in the system.
%the catalyst radius is rE, the number of catalysts is N.
%When particles overlap they are relocated apart from each other with the distance mu,
%until no partciles overlap.

%---------------------------------------------------------------------------------------------

R = (r0-(rE+0.0001));

% generating unifor distribution of points in circle as initial arrangement

n = 0;

X = zeros(1,N);
Y = zeros(1,N);

for k = 1:N

r = rand;
phi = 2*pi*rand;
    
X(k) = R*sqrt(r)*cos(phi);
Y(k) = R*sqrt(r)*sin(phi);

end

% calculation of distance of points

while n == 0

n = 1;    
    
for l = 1:N
  
 vecbd = [X(l),Y(l)]; 
 bdist = sqrt(vecbd*vecbd');
 
 if R < bdist
     
     X(l) = vecbd(1) - (vecbd(1)/bdist)*(bdist-R);
     Y(l) = vecbd(2) - (vecbd(2)/bdist)*(bdist-R);
     
     n = 0;
     
 end
    
end


for i=1:N-1
   
    for j = i+1:N
        
        vecd = [X(i)-X(j),Y(i)-Y(j)];
        
        dist = sqrt(vecd*vecd');
        
        if dist < 2*rE
           
            
           vec1 = [X(i),Y(i)];
           vec2 = [X(j),Y(j)];
           
           X(i) = vec1(1) + (vecd(1)/dist)*((1/2*((2*rE+0.0001)-dist))*(mu-mu*rand));
           Y(i) = vec1(2) + (vecd(2)/dist)*((1/2*((2*rE+0.0001)-dist))*(mu-mu*rand));
           
           X(j) = vec2(1) - (vecd(1)/dist)*((1/2*((2*rE+0.0001)-dist))*(mu-mu*rand));
           Y(j) = vec2(2) - (vecd(2)/dist)*((1/2*((2*rE+0.0001)-dist))*(mu-mu*rand));
            
           n = 0;
           
        end
        
    end
    
end

end

r = rand;
phi = 2*pi*rand;

x_cluster = (1-r0)*sqrt(r)*cos(phi);
y_cluster = (1-r0)*sqrt(r)*sin(phi);

X_shift = X+x_cluster;
Y_shift = Y+y_cluster;

x = [X_shift;Y_shift];

end





