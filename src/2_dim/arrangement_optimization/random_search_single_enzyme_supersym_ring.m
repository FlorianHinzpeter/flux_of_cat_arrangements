function optcoordiante = random_search_single_enzyme_supersym_ring(rs_iterations,N_random_starts,N_rings,N,rE,alpha1,alpha2,savepath)
%---------------------------------------------------------------------------------------------

%This function optimizes the arrangement of N second catalysts
%around a single first catalyst at the system center using a random search optimization.
%To reduce the search space the catalysts are unifmorly arrenged on N_rings.
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
X_00 = zeros(2,N);
X_opt = zeros(2,N);

    
    r = rand(1,N_rings);
    phi = 2*pi*rand(1,(N_rings-1));
    
    for t=1:N/N_rings
    
    X_00(:,t) = r(1,1)*[cos(t*2*pi/(N/N_rings));sin(t*2*pi/(N/N_rings))];
    X_00(:,t+N/N_rings) = r(1,2)*[cos(t*2*pi/(N/N_rings)+phi(1));sin(t*2*pi/(N/N_rings)+phi(1))];
    X_00(:,t+2*N/N_rings) = r(1,3)*[cos(t*2*pi/(N/N_rings)+phi(2));sin(t*2*pi/(N/N_rings)+phi(2))];
    
    end
    

cc = constraint([x1 X_00(1,:)],[y1 X_00(2,:)],1,rE,N+1);

while cc == 1

    r=rand(1,N_rings);
    phi = 2*pi*rand(1,(N_rings-1));
    
    for t=1:N/N_rings
    
    X_00(:,t) = r(1,1)*[cos(t*2*pi/(N/N_rings));sin(t*2*pi/(N/N_rings))];
    X_00(:,t+N/N_rings) = r(1,2)*[cos(t*2*pi/(N/N_rings)+phi(1));sin(t*2*pi/(N/N_rings)+phi(1))];
    X_00(:,t+2*N/N_rings) = r(1,3)*[cos(t*2*pi/(N/N_rings)+phi(2));sin(t*2*pi/(N/N_rings)+phi(2))];
    
    end

cc = constraint([x1 X_00(1,:)],[y1 X_00(2,:)],1,rE,N+1); 
    
end

X_0 = X_00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0;

eff = Flux(1,N,x1,y1,X_0(1,:),X_0(2,:),rE,alpha1,alpha2);

eff_0=eff(2);

%%% Random search %%%%%%%%%%%%%%


k = 0;
no_improve = 0;
d = 0.1;
decreasestep=150;

while n < rs_iterations

r_test = r;   
phi_test = phi;
    
n = n + 1;     
    
%%% generate test config %%%%%%%    


if no_improve == decreasestep
    
    d = d/2;
    
    no_improve = 0;
    
end
    
if d<0.0005
    break
end

c = 1;


while c == 1

if k == 0    

    Delta_r = sign(2*rand(1,N_rings)-ones(1,N_rings)).*rand(1,N_rings)*d;
    Delta_phi = sign(2*rand(1,(N_rings-1))-ones(1,(N_rings-1))).*rand(1,(N_rings-1))*2*pi*d;
       
end
    
   
    r_test = r + Delta_r;
    phi_test = phi + Delta_phi;
    

    for t=1:N/N_rings
    
    X_0(:,t) = r_test(1,1)*[cos(t*2*pi/(N/N_rings));sin(t*2*pi/(N/N_rings))];
    X_0(:,t+N/N_rings) = r_test(1,2)*[cos(t*2*pi/(N/N_rings)+phi_test(1));sin(t*2*pi/(N/N_rings)+phi_test(1))];
    X_0(:,t+2*N/N_rings) = r_test(1,3)*[cos(t*2*pi/(N/N_rings)+phi_test(2));sin(t*2*pi/(N/N_rings)+phi_test(2))];
    
    end


c = constraint([x1 X_0(1,:)],[y1 X_0(2,:)],1,rE,N+1);

if c == 1 && k == 1 
    k = 0;
end

end

eff = Flux(1,N,x1,y1,X_0(1,:),X_0(2,:),rE,alpha1,alpha2);

eff_test = eff(2);


if eff_test > eff_0
        
        eff_0 = eff_test
        X_opt = X_0
        
        r = r_test;
        phi = phi_test;
        
        
        
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