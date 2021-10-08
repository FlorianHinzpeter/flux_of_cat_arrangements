function x = Randomarrangement_pairs(rE,N,fd)
%---------------------------------------------------------------------------------------------

%This function generates a uniform arrangement of catalyst pairs a distance fd apart from each other.
%the catalyst radius is rE, the number of catalysts is N.

%---------------------------------------------------------------------------------------------

R = (1-(rE+0.001));

% generating unifor distribution of points in circle as initial arrangement

n = 0;
count = 0;

X = zeros(1,N);
Y = zeros(1,N);

XX = zeros(1,N/2);
YY = zeros(1,N/2);

theta = zeros(1,N/2);

%%% from 1->N/2 E1 enzymes, from N/2+1->N E2 enzymes

%%% random distribution of E1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:N/2

r = rand;
phi = 2*pi*rand;
    
XX(k) = R*sqrt(r)*cos(phi);
YY(k) = R*sqrt(r)*sin(phi);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while n == 0

X(1,1:N/2) = XX;
Y(1,1:N/2) = YY;
    
for k = 1:N/2

theta(k) = 2*pi*rand;

X(N/2+k) = X(k)+fd*cos(theta(k));
Y(N/2+k) = Y(k)+fd*sin(theta(k));


end


% calculation of distance of points

while n == 0

n = 1; 
count = count+1;

   %%%% distance of enzymes to boundary %%%
   
for l = 1:N/2
  
 vecbd1 = [X(l),Y(l)]; 
 bdist1 = sqrt(vecbd1*vecbd1');
 
 vecbd2 = [X(l+N/2),Y(l+N/2)]; 
 bdist2 = sqrt(vecbd2*vecbd2');
 
 if R < bdist1
     
     
     vecbd_E1 = [X(l),Y(l)];

     
     X(l) = vecbd_E1(1) - (vecbd1(1)/bdist1)*(bdist1-R);
     Y(l) = vecbd_E1(2) - (vecbd1(2)/bdist1)*(bdist1-R);
     
     X(l+N/2) = X(l) + fd*cos(theta(l));
     Y(l+N/2) = Y(l) + fd*sin(theta(l));
     
     n = 0;
     
 elseif R < bdist2
         
     vecbd_E1 = [X(l),Y(l)];
     
     X(l) = vecbd_E1(1) - (vecbd2(1)/bdist2)*(bdist2-R);
     Y(l) = vecbd_E1(2) - (vecbd2(2)/bdist2)*(bdist2-R);
     
     X(l+N/2) = X(l) + fd*cos(theta(l));
     Y(l+N/2) = Y(l) + fd*sin(theta(l));
     
     n = 0;
     
 end
    
end

%%% enzyme-enzyme distance %%%%%%

for i=1:N/2-1
   
    for j = i+1:N/2
        
        vecd1 = [X(i)-X(j),Y(i)-Y(j)];
        dist1 = sqrt(vecd1*vecd1');
        
        vecd2 = [X(i+N/2)-X(j),Y(i+N/2)-Y(j)];
        dist2 = sqrt(vecd2*vecd2');
        
        vecd3 = [X(i)-X(j+N/2),Y(i)-Y(j+N/2)];
        dist3 = sqrt(vecd3*vecd3');
        
        vecd4 = [X(i+N/2)-X(j+N/2),Y(i+N/2)-Y(j+N/2)];
        dist4 = sqrt(vecd4*vecd4');
        
        if dist1 < 2*rE
           
           vec_E11 = [X(i),Y(i)];
           vec_E12 = [X(j),Y(j)];
           
           
           X(i) = vec_E11(1) + (vecd1(1)/dist1)*((1/2*((2*rE+0.0001)-dist1))*(2-2*rand));
           Y(i) = vec_E11(2) + (vecd1(2)/dist1)*((1/2*((2*rE+0.0001)-dist1))*(2-2*rand));
           
           X(j) = vec_E12(1) - (vecd1(1)/dist1)*((1/2*((2*rE+0.0001)-dist1))*(2-2*rand));
           Y(j) = vec_E12(2) - (vecd1(2)/dist1)*((1/2*((2*rE+0.0001)-dist1))*(2-2*rand));
           
           
           X(i+N/2) = X(i) + fd*cos(theta(i));
           Y(i+N/2) = Y(i) + fd*sin(theta(i));
           
           X(j+N/2) = X(j) + fd*cos(theta(j));
           Y(j+N/2) = Y(j) + fd*sin(theta(j));
           
            
           n = 0;
           
        elseif dist2 <2*rE
            
           vec_E11 = [X(i),Y(i)];
           vec_E12 = [X(j),Y(j)];
           
           X(i) = vec_E11(1) + (vecd2(1)/dist2)*((1/2*((2*rE+0.0001)-dist2))*(2-2*rand));
           Y(i) = vec_E11(2) + (vecd2(2)/dist2)*((1/2*((2*rE+0.0001)-dist2))*(2-2*rand));
           
           X(j) = vec_E12(1) - (vecd2(1)/dist2)*((1/2*((2*rE+0.0001)-dist2))*(2-2*rand));
           Y(j) = vec_E12(2) - (vecd2(2)/dist2)*((1/2*((2*rE+0.0001)-dist2))*(2-2*rand));
           
           
           X(i+N/2) = X(i) + fd*cos(theta(i));
           Y(i+N/2) = Y(i) + fd*sin(theta(i));
           
           X(j+N/2) = X(j) + fd*cos(theta(j));
           Y(j+N/2) = Y(j) + fd*sin(theta(j));
            
           n = 0;
           
           elseif dist3 <2*rE
            
           vec_E11 = [X(i),Y(i)];
           vec_E12 = [X(j),Y(j)];
           
           
           X(i) = vec_E11(1) + (vecd3(1)/dist3)*((1/2*((2*rE+0.0001)-dist3))*(2-2*rand));
           Y(i) = vec_E11(2) + (vecd3(2)/dist3)*((1/2*((2*rE+0.0001)-dist3))*(2-2*rand));
           
           X(j) = vec_E12(1) - (vecd3(1)/dist3)*((1/2*((2*rE+0.0001)-dist3))*(2-2*rand));
           Y(j) = vec_E12(2) - (vecd3(2)/dist3)*((1/2*((2*rE+0.0001)-dist3))*(2-2*rand));
           
           
           X(i+N/2) = X(i) + fd*cos(theta(i));
           Y(i+N/2) = Y(i) + fd*sin(theta(i));
           
           X(j+N/2) = X(j) + fd*cos(theta(j));
           Y(j+N/2) = Y(j) + fd*sin(theta(j));
            
           n = 0;
           
           elseif dist4 < 2*rE
            
           vec_E11 = [X(i),Y(i)];
           vec_E12 = [X(j),Y(j)];
           
           
           X(i) = vec_E11(1) + (vecd4(1)/dist4)*((1/2*((2*rE+0.0001)-dist4))*(2-2*rand));
           Y(i) = vec_E11(2) + (vecd4(2)/dist4)*((1/2*((2*rE+0.0001)-dist4))*(2-2*rand));
           
           X(j) = vec_E12(1) - (vecd4(1)/dist4)*((1/2*((2*rE+0.0001)-dist4))*(2-2*rand));
           Y(j) = vec_E12(2) - (vecd4(2)/dist4)*((1/2*((2*rE+0.0001)-dist4))*(2-2*rand));
           
           
           X(i+N/2) = X(i) + fd*cos(theta(i));
           Y(i+N/2) = Y(i) + fd*sin(theta(i));
           
           X(j+N/2) = X(j) + fd*cos(theta(j));
           Y(j+N/2) = Y(j) + fd*sin(theta(j));
            
           n = 0;
           
        end
        
    end
    
end

%%%% Stop iteration when number of interations >1000

if count > 1000
    
    n = 0;
    count = 0;
    
    break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

end


x = [X;Y];

end





