function x = Randomarrangement_dynamic(r0,rE,N)
%---------------------------------------------------------------------------------------------

%This function generates a uniform arrangement of catalysts in a system with radius r0.
%the catalyst radius is rE, the number of catalysts is N.
%When particles overlap they are relocated apart from each other
%until no partciles overlap.

%---------------------------------------------------------------------------------------------

R = (r0-(rE+0.0001));

% generating unifor distribution of points in circle

n = 0;

X = zeros(1,N);
Y = zeros(1,N);
Z = zeros(1,N);

for k = 1:N

r = rand;
phi = 2*pi*rand;
theta = acos(2*rand-1);
    
X(k) = R*r^(1/3)*cos(phi)*sin(theta);
Y(k) = R*r^(1/3)*sin(phi)*sin(theta);
Z(k) = R*r^(1/3)*cos(theta);

end

% calculation of distance of points

while n == 0

n = 1;    
    
for l = 1:N
  
 vecbd = [X(l),Y(l),Z(l)]; 
 bdist = sqrt(vecbd*vecbd');
 
 if R < bdist
     
     X(l) = vecbd(1) - (vecbd(1)/bdist)*(bdist-R);
     Y(l) = vecbd(2) - (vecbd(2)/bdist)*(bdist-R);
     Z(l) = vecbd(3) - (vecbd(3)/bdist)*(bdist-R);
     
     n = 0;
     
 end
    
end


for i=1:N-1
   
    for j = i+1:N
        
        vecd = [X(i)-X(j),Y(i)-Y(j),Z(i)-Z(j)];
        
        dist = sqrt(vecd*vecd');
        
        if dist < 2*rE
           
            
           vec1 = [X(i),Y(i),Z(i)];
           vec2 = [X(j),Y(j),Z(j)];
           
           X(i) = vec1(1) + (vecd(1)/dist)*((1/2*((2*rE+0.0001)-dist))*(3-3*rand));
           Y(i) = vec1(2) + (vecd(2)/dist)*((1/2*((2*rE+0.0001)-dist))*(3-3*rand));
           Z(i) = vec1(3) + (vecd(3)/dist)*((1/2*((2*rE+0.0001)-dist))*(3-3*rand));
           
           X(j) = vec2(1) - (vecd(1)/dist)*((1/2*((2*rE+0.0001)-dist))*(3-3*rand));
           Y(j) = vec2(2) - (vecd(2)/dist)*((1/2*((2*rE+0.0001)-dist))*(3-3*rand));
           Z(j) = vec2(3) - (vecd(3)/dist)*((1/2*((2*rE+0.0001)-dist))*(3-3*rand));
            
           n = 0;
           
        end
        
    end
    
end

end

x = [X;Y;Z];

end





