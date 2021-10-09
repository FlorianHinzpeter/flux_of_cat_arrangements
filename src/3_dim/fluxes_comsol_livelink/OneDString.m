function X = OneDString(N1,N2,NStrings,savepath)

N = N1+N2;
X = zeros(3,N);

n = N/NStrings;

r0 = 0.03;
x = -r0/2;

for i = 1:n/2
    
    x = x + r0;
    
    Coord1(:,i) = [x;0;0];
    
end

CoordString1 = [Coord1 -Coord1];

X(:,1:length(CoordString1)) = CoordString1;

y0 = 1.5/NStrings;

for j = 1:NStrings/2

   X(1,(j-1)*length(CoordString1)+1:j*length(CoordString1)) = CoordString1(1,:);
   X(2,(j-1)*length(CoordString1)+1:j*length(CoordString1)) = -ones(1,length(CoordString1))*y0/2 + ones(1,length(CoordString1))*j*y0;
   X(3,(j-1)*length(CoordString1)+1:j*length(CoordString1)) = 0;
    
end


X(1,(NStrings/2)*length(CoordString1)+1:NStrings*length(CoordString1)) = X(1,1:(NStrings/2)*length(CoordString1));
X(2,(NStrings/2)*length(CoordString1)+1:NStrings*length(CoordString1)) = -X(2,1:(NStrings/2)*length(CoordString1));
X(3,(NStrings/2)*length(CoordString1)+1:NStrings*length(CoordString1)) = 0;

save(savepath,'X')

end