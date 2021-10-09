function optcoordiante = random_search_single_enzyme_ring(rs_iterations,N_random_starts,N,rE,alpha1,alpha2,savepath)
%---------------------------------------------------------------------------------------------

%This function optimizes the radius of a ring with N second catalysts uniformly arranged
%around a single first catalyst at the system center.
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
X_0 = zeros(2,N);
X_opt = zeros(2,N);

    
    X_sym=0.4;
    
    for t=1:N
    
    X_0(:,t) = X_sym*[cos(t*2*pi/N);sin(t*2*pi/N)];
    
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0;

eff = Flux(1,N,x1,y1,X_0(1,:),X_0(2,:),rE,alpha1,alpha2);

eff_0=eff(2);

%%% Random search %%%%%%%%%%%%%%


k = 0;
no_improve = 0;
d = 0.1;
decreasestep = 10;

while n < rs_iterations

X_test=X_sym;    
    
n = n + 1;     
    
%%% generate test config %%%%%%%    



if no_improve == decreasestep
    
    d = d/2;
    no_improve = 0;
    
end
    
if d<0.0001
    break
end

c = 1;


while c == 1

if k == 0
    
Delta = sign(2*rand-1)*d;

end

X_test = X_sym + Delta;

    for t=1:N
    
    X_0(:,t) = X_test*[cos(t*2*pi/N);sin(t*2*pi/N)];
    
    end

c = constraint([x1 X_0(1,:)],[y1 X_0(2,:)],1,rE,N+1);

if c == 1
    
    k = 0;
    
end

end

eff = Flux_bd_loss(1,N,x1,y1,X_0(1,:),X_0(2,:),rE,alpha1,alpha2);

eff_test = eff(2);


if eff_test > eff_0
        
        eff_0 = eff_test
        
        X_sym = X_test
        
    for t=1:N
    
    X_opt(:,t) = X_test*[cos(t*2*pi/N);sin(t*2*pi/N)];
    
    end
        
        k = 1;
        
        no_improve = 0;
else
    
    k = 0;
    no_improve = no_improve+1
end

end


J(2*l-1:2*l,1:N) = X_opt;
J(2*l-1:2*l,N+1) = Flux_bd_loss(1,N,x1,y1,X_opt(1,:),X_opt(2,:),rE,alpha1,alpha2);

save(savepath,'J')

end

optcoordiante = J;

end