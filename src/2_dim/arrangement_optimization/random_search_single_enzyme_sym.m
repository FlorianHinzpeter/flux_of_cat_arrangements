function optcoordiante = random_search_single_enzyme_sym(rs_iterations,N_random_starts,N,rE,alpha1,alpha2,savepath)

%---------------------------------------------------------------------------------------------

%This function optimizes the arrangement of N second catalysts
%around a single first catalyst at the system center using a random search optimization.
%To reduce the search space the catalysts arrangement is constrained to be mirror symmetric.
%rs_iterations is the number of random steps,
%N_random_starts is the number of individual optimization runs
%rE is the catalyst radius in units of the system size
%alpha1, and alpha2 are the dimensionless reaction-diffusion parameter of the first catalyst.
%savepath is the path and file name where the results are saved

%---------------------------------------------------------------------------------------------

rng('shuffle')

J = zeros(2*N_random_starts,N+1);

for l=1:N_random_starts
% choose random initial E2 concentration

x1 = 0;
y1 = 0;

if mod(N,2) == 0
    
    x=[[0;0.1],sign(2*rand(2,(N-2)/2)-ones(2,(N-2)/2)).*rand(2,(N-2)/2).*0.5];
    X_00 = [x(1,:),x(1,:);x(2,:),-x(2,:)];
    
else
    
   x=[[0;0.1],sign(2*rand(2,(N-1)/2)-ones(2,(N-1)/2)).*rand(2,(N-1)/2).*0.5];
   X_00 = [[x(1,1);x(2,1)],[x(1,2:(N+1)/2),-x(1,2:(N+1)/2);x(2,2:(N+1)/2),x(2,2:(N+1)/2)]];
    
end


cc = constraint([x1 X_00(1,:)],[y1 X_00(2,:)],1,rE,N+1);

while cc == 1

if mod(N,2) == 0
    
    x=[[0;0.1],sign(2*rand(2,(N-2)/2)-ones(2,(N-2)/2)).*rand(2,(N-2)/2).*0.5];
    X_00 = [x(1,:),x(1,:);x(2,:),-x(2,:)];
    
else
    
   x=[[0;0.1],sign(2*rand(2,(N-1)/2)-ones(2,(N-1)/2)).*rand(2,(N-1)/2).*0.5];
   X_00 = [[x(1,1);x(2,1)],[x(1,2:(N+1)/2),-x(1,2:(N+1)/2);x(2,2:(N+1)/2),x(2,2:(N+1)/2)]];
    
end

cc = constraint([x1 X_00(1,:)],[y1 X_00(2,:)],1,rE,N+1); 
    
end

X_sym = x;

if mod(N,2) == 0
    
   X_0 = [X_sym(1,:),X_sym(1,:);X_sym(2,:),-X_sym(2,:)];
    
else
    
   X_0 = [[X_sym(1,1);X_sym(2,1)],[X_sym(1,2:(N+1)/2),-X_sym(1,2:(N+1)/2);X_sym(2,2:(N+1)/2),X_sym(2,2:(N+1)/2)]];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0;

eff = Flux_bd_loss(1,N,x1,y1,X_0(1,:),X_0(2,:),rE,alpha1,alpha2);

eff_0=eff(2);

%%% Random search %%%%%%%%%%%%%%


k = 0;
no_improve = 0;
d = 0.05;

while n < rs_iterations

X_test=X_sym;    
    
n = n + 1;     
    
%%% generate test config %%%%%%%    


if no_improve == 150
    
    d = d/2;
    no_improve = 0;
    
end
    
if d<0.0005
    break
end

c = 1;


while c == 1

if k == 0    

    if mod(N,2)==0
        
   enzyme_number = randi(N/2); 
    
    else
        
    enzyme_number = randi((N+1)/2);
    
    end
        
if enzyme_number == 1
 
     Delta = [0;sign(2*rand(1,1)-1).*d];
 
else
    
delta = sign(2*rand(2,1)-ones(2,1)).*rand(2,1);
Delta = delta./(sqrt(delta(1)^2+delta(2)^2)).*d;
end

end

X_test(:,enzyme_number) = X_sym(:,enzyme_number) + Delta;

if mod(N,2) == 0
    
   X_0 = [X_test(1,:),X_test(1,:);X_test(2,:),-X_test(2,:)];
    
else
    
   X_0 = [[X_test(1,1);X_test(2,1)],[X_test(1,2:(N+1)/2),-X_test(1,2:(N+1)/2);X_test(2,2:(N+1)/2),X_test(2,2:(N+1)/2)]];
    
end

c = constraint([x1 X_0(1,:)],[y1 X_0(2,:)],1,rE,N+1);

if c == 1 && k == 1 
    k = 0;
end

end

eff = Flux_bd_loss(1,N,x1,y1,X_0(1,:),X_0(2,:),rE,alpha1,alpha2);

eff_test = eff(2);


if eff_test > eff_0
        
        eff_0 = eff_test
        
        X_sym = X_test;
        
if mod(N,2) == 0
    
   X_opt = [X_sym(1,:),X_sym(1,:);X_sym(2,:),-X_sym(2,:)]
    
else
    
   X_opt = [[X_sym(1,1);X_sym(2,1)],[X_sym(1,2:(N+1)/2),-X_sym(1,2:(N+1)/2);X_sym(2,2:(N+1)/2),X_sym(2,2:(N+1)/2)]]
    
end
        
        k = 1;
        
        no_improve = 0;
else
    
    k = 0;
    no_improve = no_improve+1
end

end


J(2*l-1:2*l,1:N) = X_opt;
J(2*l-1:2*l,N+1) = Flux(1,N,x1,y1,X_opt(1,:),X_opt(2,:),rE,alpha1,alpha2);

save(savepath,'J')

end

optcoordiante = J;

end