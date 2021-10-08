function X = RandomSearchOnShell(rs_iterations,N,rE,alpha1,alpha2,R)

%---------------------------------------------------------------------------------------------

%This function computes the optimal arrangement of N second catalysts around a single first catalyst
%at the center. The arrangements are constrained onto a spherical shell with radius R.
%alpha1 and alpha2 are the dimension less reaction diffusion parameters,
%rE is the catalyst radius, and rs_iterations is the number of random search iterations

%---------------------------------------------------------------------------------------------

X1 = [0;0;0];
cc = 1;
X_0 = zeros(3,N);
X_test = zeros(3,N);

while cc == 1

    PHI = 2*pi*rand(1,N);
    THETA = pi*rand(1,N);
    
    PHI(N) = 0;
    PHI(N-1) = 0;
    THETA(N) = 0;
    
    for i = 1:N
        
       X_0(1,i) = R*sin(THETA(i))*cos(PHI(i));
       X_0(2,i) = R*sin(THETA(i))*sin(PHI(i));
       X_0(3,i) = R*cos(THETA(i));
        
    end

    cc = constraint([X1(1) X_0(1,:)],[X1(2) X_0(2,:)],[X1(3) X_0(3,:)],1,rE+rE/20,N+1);

end

n = 0;

eff = Flux_bd_loss(1,N,X1(1),X1(2),X1(3),X_0(1,1:N),X_0(2,1:N),X_0(3,1:N),rE,alpha1,alpha2)

eff_0=eff(2);

%%% Random search %%%%%%%%%%%%%%


k = 0;
no_improve = 0;
d = 0.1;

while n < rs_iterations
    
    n = n + 1; 
    
    PHI_test = PHI;
    THETA_test = THETA;
        
    
%%% generate test config %%%%%%%
    
    
    if no_improve == 100
    
        d = d/2;
        no_improve = 0;
    
    end
    

        if k == 0    
            
            enzyme_number = randi(N-1);
            
            if enzyme_number == N-1
                    
                Delta_PHI = 0;
                Delta_THETA = pi*sign(2*rand-1)*d;
                    
            else    
                
            Delta_PHI = 2*pi*sign(2*rand-1)*d;
            Delta_THETA = pi*sign(2*rand-1)*d;
            
            end

        end
        
    PHI_test(enzyme_number) = PHI(enzyme_number)+Delta_PHI;
    THETA_test(enzyme_number) = THETA(enzyme_number)+Delta_THETA;
    
    for i = 1:N
        
       X_test(1,i) = R*sin(THETA_test(i))*cos(PHI_test(i));
       X_test(2,i) = R*sin(THETA_test(i))*sin(PHI_test(i));
       X_test(3,i) = R*cos(THETA_test(i));
        
    end
    
    c = constraint([X1(1) X_test(1,:)],[X1(2) X_test(2,:)],[X1(3) X_test(3,:)],1,rE+rE/20,N+1);
    
    if c == 0

        eff = Flux_bd_loss(1,N,X1(1),X1(2),X1(3),X_test(1,1:N),X_test(2,1:N),X_test(3,1:N),rE,alpha1,alpha2);

        eff_test = eff(2);

    else
    
        eff_test = 0;
        
    end
    
    
    if eff_test > eff_0
        
        eff_0 = eff_test
        
        PHI = PHI_test;
        THETA = THETA_test;
        
        k = 1;
        
        no_improve = 0;
        
        X = X_test;
        
        save(sprintf('Optimization_On_Sphere_Surf_N=%d_alpha1=%d_alpha2=%d.mat',N,alpha1,alpha2),'X')
        
    else
    
        k = 0;
        no_improve = no_improve+1
    
    end   
    
end