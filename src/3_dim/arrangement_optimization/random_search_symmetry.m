function optcoordiante = random_search_symmetry(rs_iterations,N,rE,alpha1,alpha2)
%---------------------------------------------------------------------------------------------

%This function computes the optimal arrangement of N second catalysts around a single first catalyst
%at the center. The starting configuration is the Thomson config on a shell with optimal radius
%alpha1 and alpha2 are the dimension less reaction diffusion parameters,
%rE is the catalyst radius, and rs_iterations is the number of random search iterations

%---------------------------------------------------------------------------------------------

X1 = [0;0;0];

load('Coordiantes_snub_cube.mat')

cc = 1;

X_0 = zeros(size(Coord));
X_test = zeros(size(Coord));

while cc == 1

    r = 0.4*rand(1,length(Coord(1,:)));
    
  for j = 1:length(Coord(1,:))
     
      X_0(:,j) = r(j)*Coord(:,j);
      
  end
  
  cc = constraint([X1(1) X_0(1,:)],[X1(2) X_0(2,:)],[X1(3) X_0(3,:)],1,rE+rE/20,N+1);
  
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0;

eff = Flux_bd_loss(1,N,X1(1),X1(2),X1(3),X_0(1,1:N),X_0(2,1:N),X_0(3,1:N),rE,alpha1,alpha2)

eff_0=eff(2);

%%% Random search %%%%%%%%%%%%%%


k = 0;
no_improve = 0;
d = 0.1;

while n < rs_iterations

r_test=r;    
    
n = n + 1;     
    
%%% generate test config %%%%%%%    


if no_improve == 100
    
    d = d/10;
    no_improve = 0;
    
end
    


if k == 0    

enzyme_number = randi(N);    
Delta = sign(2*rand-1)*d;

end

r_test(:,enzyme_number) = r(:,enzyme_number) + Delta;

  for j = 1:length(Coord(1,:))
     
      X_test(:,j) = r_test(j)*Coord(:,j);
      
  end

c = constraint([X1(1) X_test(1,:)],[X1(2) X_test(2,:)],[X1(3) X_test(3,:)],1,rE+rE/20,N+1);

if c == 0

eff = Flux_bd_loss(1,N,X1(1),X1(2),X1(3),X_test(1,1:N),X_test(2,1:N),X_test(3,1:N),rE,alpha1,alpha2)

eff_test = eff(2);

else
    
    eff_test = 0;
end

if eff_test > eff_0
        
        eff_0 = eff_test
        
        r = r_test
        
        
        k = 1;
        
        no_improve = 0;
        
else
    
    k = 0;
    no_improve = no_improve+1
    
end


optcoordiante = X_0;
    

end