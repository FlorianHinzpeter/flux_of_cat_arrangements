function x = FiveFusedArrangement(trials,rE,N1,N2,fd,savepath)
%---------------------------------------------------------------------------------------------

%This function generates arrangements where 5 catalysts are fused togetehr in a symmetric manner.
%The distance between the catalysts is fd.

%---------------------------------------------------------------------------------------------

x = zeros(2*trials,N1+N2);

for h = 1:trials

R = (1-(rE+0.005));

X = zeros(1,N1+N2);
Y = zeros(1,N1+N2);

n = 1;

while n == 1

for k = 1:N2

r = rand;
phi = 2*pi*rand;
    
X(k) = R*sqrt(r)*cos(phi);
Y(k) = R*sqrt(r)*sin(phi);

rot = 2*pi/rand;


X(N2+4*k-3) = X(k)+fd*cos(2*pi/4+rot);
Y(N2+4*k-3) = Y(k)+fd*sin(2*pi/4+rot);

X(N2+4*k-2) = X(k)+fd*cos(4*pi/4+rot);
Y(N2+4*k-2) = Y(k)+fd*sin(4*pi/4+rot);

X(N2+4*k-1) = X(k)+fd*cos(6*pi/4+rot);
Y(N2+4*k-1) = Y(k)+fd*sin(6*pi/4+rot);

X(N2+4*k) = X(k)+fd*cos(8*pi/4+rot);
Y(N2+4*k) = Y(k)+fd*sin(8*pi/4+rot);

end

n = constraint(X,Y,1,rE,N1+N2);

end

x(2*h-1:2*h,:) = [X;Y];

end

save(savepath,'x')

end