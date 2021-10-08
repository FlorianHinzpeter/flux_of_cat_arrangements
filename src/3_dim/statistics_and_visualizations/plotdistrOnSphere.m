function plotdistrOnSphere(x1,x2,N1,N2,rE,R)

[x,y,z] = sphere;

figure

axis([-1 1 -1 1 -1 1])

mesh(R*x,R*y,R*z,'edgecolor','k','facecolor','none');

hold on

for i = 1:N1
    
    surf(rE*x+x1(1,i),rE*y+x1(2,i),rE*z+x1(3,i),'edgecolor','b','facecolor','none')
    
end

for i = 1:N2
    
    surf(rE*x+x2(1,i),rE*y+x2(2,i),rE*z+x2(3,i),'edgecolor','r','facecolor','none')
    
end

pbaspect([1 1 1])

end