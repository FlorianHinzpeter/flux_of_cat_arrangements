function optcoordiante = random_search_single_enzyme(rs_iterations,N_random_starts,N,rE,alpha1,alpha2,savepath)
%---------------------------------------------------------------------------------------------

%This function optimizes the arrangement of N second catalysts
%around a single first catalyst at the system center using a random search optimization.
%rs_iterations is the number of random steps,
%N_random_starts is the number of individual optimization runs
%rE is the catalyst radius in units of the system size.catalyst
%alpha1, and alpha2 are the dimensionless reaction-diffusion parameters
%savepath is the path and file name where the results are saved

%---------------------------------------------------------------------------------------------

rng('shuffle')

J = zeros(2*N_random_starts,N+1);

for l=1:N_random_starts
% choose random initial E2 concentration

x1 = 0;
y1 = 0;

X_00 = Randomarrangement_dynamic(0.2,rE,N-1,3);

cc = constraint([x1 0 X_00(1,:)],[y1 0.3 X_00(2,:)],1,rE,N+1);

while cc == 1

X_00 = Randomarrangement_dynamic(0.2,rE,N-1,3);

cc = constraint([x1 0 X_00(1,:)],[y1 0.3 X_00(2,:)],1,rE,N+1);    
    
end

X_0 = [[0;0.3],X_00];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0;

eff = Flux(1,N,x1,y1,X_0(1,:),X_0(2,:),rE,alpha1,alpha2);

eff_0=eff(2);

%%% Random search %%%%%%%%%%%%%%


k = 0;
no_improve = 0;
d = 0.05;

while n < rs_iterations

X_test=X_0;    
    
n = n + 1;     
    
%%% generate test config %%%%%%%    

if no_improve == 250
    
    d = d/2;
    no_improve = 0;
    
end
 
if d<0.0005
    break
end

c = 1;


while c == 1

if k == 0    

   enzyme_number = randi(N); 
    
if enzyme_number == 1
 
Delta = [0;sign(2*rand(1,1)-1).*d];    
    
else
    
delta = sign(2*rand(2,1)-ones(2,1)).*rand(2,1);
Delta = delta./(sqrt(delta(1)^2+delta(2)^2)).*d;
end

end

X_test(:,enzyme_number) = X_0(:,enzyme_number) + Delta;

c = constraint([x1 X_test(1,:)],[y1 X_test(2,:)],1,rE,N+1);

if c == 1 && k == 1 
    k = 0;
end

end

eff = Flux(1,N,x1,y1,X_test(1,:),X_test(2,:),rE,alpha1,alpha2);

eff_test = eff(2);


if eff_test > eff_0
        
        eff_0 = eff_test
        
        X_0 = X_test
        
        
        k = 1;
        
        no_improve = 0;
        
else
    
    k = 0;
    no_improve = no_improve+1
    
end

end


J(2*l-1:2*l,1:N) = X_0;
J(2*l-1:2*l,N+1) = Flux(1,N,x1,y1,X_0(1,:),X_0(2,:),rE,alpha1,alpha2);

save(savepath,'J')

end

optcoordiante = J;

end