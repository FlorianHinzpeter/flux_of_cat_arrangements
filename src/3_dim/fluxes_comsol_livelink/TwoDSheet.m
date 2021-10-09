function X = TwoDSheet(N1,N2,savepath)

N = N1+N2;

Coord = zeros(3,51*51);

r0 = 0.03;
y = -r0;

for i = 1:51
    
    x = 0;
    y = y+r0;
    
    for j = 1:51
        
    Coord(:,(i-1)*51 + j) =   [x;y;0];  
        
    x = x+r0;    
        
    end
    
end

CenterOfMass(1) = sum(Coord(1,:))/(51*51);
CenterOfMass(2) = sum(Coord(2,:))/(51*51);
CenterOfMass(3) = 0;

for l= 1:51*51

ShiftedCoord(:,l) = Coord(:,l) - CenterOfMass';

end

ShiftedCoord1 = ShiftedCoord(:,1:2:end);
ShiftedCoord2 = ShiftedCoord(:,2:2:end);

d1 = zeros(1,length(ShiftedCoord1));

d2 = zeros(1,length(ShiftedCoord2));

for l = 1:length(ShiftedCoord1)
    
   d1(l) = sqrt(ShiftedCoord1(1,l)^2+ShiftedCoord1(2,l)^2); 
   
end

for l = 1:length(ShiftedCoord2)
   
   d2(l) = sqrt(ShiftedCoord2(1,l)^2+ShiftedCoord2(2,l)^2); 
    
end



[D1,I1] = sort(d1);

[D2,I2] = sort(d2);

SortedCoord1 = ShiftedCoord1(:,I1);
SortedCoord2 = ShiftedCoord2(:,I2);

X = [SortedCoord1(:,1:N1) SortedCoord2(:,1:N2)];


save(savepath,'X')

end