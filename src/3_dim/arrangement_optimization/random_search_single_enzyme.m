function optcoordiante = random_search_single_enzyme(rs_iterations,N,rE,alpha1,alpha2)
%---------------------------------------------------------------------------------------------

%This function computes the optimal arrangement of N second catalysts around a single first catalyst
%at the center. alpha1 and alpha2 are the dimension less reaction diffusion parameters,
%rE is the catalyst radius, and rs_iterations is the number of random search iterations.

%---------------------------------------------------------------------------------------------

% choose random initial E2 concentration
X1 = [0;0;0];


%X_00 = Randomarrangement_dynamic(0.2,rE+rE/10,N+1);

%cc = constraint(X_00(1,:),X_00(2,:),X_00(3,:),1,rE+rE/10,N+1);

%while cc == 1

%X_00 = Randomarrangement_dynamic(0.2,rE+rE/10,N+1);

%cc = constraint(X_00(1,:),X_00(2,:),X_00(3,:),1,rE+rE/10,N+1);    
    
%end

%X_0 = X_00(:,2:N+1);

load('Coordiantes_snub_cube.mat')

Coord = coord;

X_0 = 0.03*Coord;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0;

eff = Flux_bd_loss(1,N,X1(1),X1(2),X1(3),X_0(1,1:N),X_0(2,1:N),X_0(3,1:N),0.01,alpha1,alpha2)

eff_0=eff(2);

%%% Random search %%%%%%%%%%%%%%


k = 0;
no_improve = 0;
d = 0.01;

while n < rs_iterations

X_test=X_0;    
    
n = n + 1;     
    
%%% generate test config %%%%%%%    


if no_improve == 100
    
    d = d/10;
    no_improve = 0;
    
end
    

c = 1;

enzyme_number = randi(N);



if k == 0    
   
Delta = sign(2*rand(3,1)-ones(3,1)).*rand(3,1)*d;

end

X_test(:,enzyme_number) = X_0(:,enzyme_number) + Delta;

c = constraint([X1(1) X_test(1,:)],[X1(2) X_test(2,:)],[X1(3) X_test(3,:)],1,rE+rE/20,N+1)

if c==0

eff = Flux_bd_loss(1,N,X1(1),X1(2),X1(3),X_test(1,1:N),X_test(2,1:N),X_test(3,1:N),0.01,alpha1,alpha2)

eff_test = eff(2);

else
    
    eff_test = 0;
end

if eff_test > eff_0
        
        eff_0 = eff_test
        
        X_0 = X_test
        
        
        k = 1;
        
        no_improve = 0;
        
else
    
    k = 0;
    no_improve = no_improve+1
    
end


optcoordiante = X_0;
    

end