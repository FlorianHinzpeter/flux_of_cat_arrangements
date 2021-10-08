function plotdistrNoSys(x1,x2,N1,N2,rE)

[x,y,z] = sphere;

figure

axis([-0.1 0.1 -0.1 0.1 -0.1 0.1])

%mesh(x,y,z,'edgecolor','k','facecolor','none');

%hold on

for i = 1:N1
    
    surf(rE*x+x1(1,i),rE*y+x1(2,i),rE*z+x1(3,i),'edgecolor','b','facecolor','none')
    hold on
end

for i = 1:N2
    
    surf(rE*x+x2(1,i),rE*y+x2(2,i),rE*z+x2(3,i),'edgecolor','r','facecolor','none')
    
end

pbaspect([0.1 0.1 0.1])

end