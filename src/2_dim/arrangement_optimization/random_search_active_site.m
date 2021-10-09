function optcoordiante = random_search_active_site(rs_iterations,N_random_starts,N,rE,alpha1,alpha2,save_path)
%---------------------------------------------------------------------------------------------

%This function optimizes the arrangement of N second catalysts and the orientation of their active site
%around a single first catalyst at the system center using a random search optimization.
%rs_iterations is the number of random steps,
%N_random_starts is the number of individual optimization runs
%rE is the catalyst radius in units of the system size.catalyst
%alpha1, and alpha2 are the dimensionless reaction-diffusion parameters
%save_path is the path and file name where the results are saved

%---------------------------------------------------------------------------------------------

rng('shuffle')

J = zeros(3*N_random_starts,N+1);
patch = 1/6;

for l=1:N_random_starts
% choose random initial E2 concentration

x1 = 0;
y1 = 0;

cc = 1;

while cc == 1

X_00 = Randomarrangement_dynamic(0.2,rE,N,3);

cc = constraint([x1 X_00(1,:)],[y1 X_00(2,:)],1,rE,N+1);    
    
end

X_0 = X_00;
O_0 = rand(1,N)*2*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0;

eff = Flux_active_site(1,N,x1,y1,X_0(1,:),X_0(2,:),[-2*pi*patch/2 O_0],patch,rE,alpha1/patch,alpha2/patch);

eff_0=eff(2);

%%% Random search %%%%%%%%%%%%%%


k = 0;
no_improve = 0;
d = 0.005;

while n < rs_iterations

X_test = X_0;
O_test = O_0;
    
n = n + 1;     
    
%%% generate test config %%%%%%%    

if no_improve == 250
    
    d = d/4;
    no_improve = 0;
    
end
 
if d<0.0001
    break
end

c = 1;


while c == 1

if k == 0    

    dis = rand;
    
    if dis>0
        
    enzyme_number_X = randi(N);   
    delta_X = (sign(2*rand(2,1)-ones(2,1)).*rand(2,1));
    Delta_X = delta_X./(sqrt(delta_X(1)^2+delta_X(2)^2)).*d;    
    
    else
      
    enzyme_number_O = randi(N);  
    Delta_O = sign(2*rand-1)*rand*2*pi*d;
        
    end

end

    if dis>0

    X_test(:,enzyme_number_X) = X_0(:,enzyme_number_X) + Delta_X;

    else
        
    O_test(:,enzyme_number_O) = O_0(:,enzyme_number_O) + Delta_O;
    
    end
    
    

c = constraint([x1 X_test(1,:)],[y1 X_test(2,:)],1,rE,N+1);

if c == 1 && k == 1 
    k = 0;
end

end


eff = Flux_active_site(1,N,x1,y1,X_test(1,:),X_test(2,:),[-2*pi*patch/2 O_test],patch,rE,alpha1/patch,alpha2/patch);

eff_test = eff(2);


if eff_test > eff_0
        
        eff_0 = eff_test
        
        X_0 = X_test
        O_0 = O_test
        
        
        k = 1;
        
        no_improve = 0;
        
else
    
    k = 0;
    no_improve = no_improve+1
    
end

end


J(3*l-2:3*l-1,1:N) = X_0;
J(3*l,1:N) = O_0;
J(3*l-2:3*l,N+1) = eff_0;

save(save_path,'J')

end

optcoordiante = J;

end