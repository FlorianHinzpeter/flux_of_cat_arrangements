function J = DelocalizedFlux(trials,N1,N2,alpha,load_arrang_path,flux_save_path)
%---------------------------------------------------------------------------------------------

%This function loads the delocalized arrangements from the path, load_arrang_path,
%and computes for each arrangement the flux using Flux_bd_loss.
%Subsequently, the arrangements are saved back using the path and file name
%specified in flux_save_path.

%---------------------------------------------------------------------------------------------

load(load_arrang_path);
 J = zeros(trials,2);
 
 for i = 1:trials
     
     try
     
    i
    J(i,:) = Flux_bd_loss(N1,N2,X(3*i-2,1:2:(N1+N2)),X(3*i-1,1:2:(N1+N2)),X(3*i,1:2:(N1+N2)),X(3*i-2,2:2:(N1+N2)),X(3*i-1,2:2:(N1+N2)),X(3*i,2:2:(N1+N2)),0.01,alpha,alpha);
    
     catch
         
       J(i,:)=0;
       
     end
    
    save(flux_save_path,'J');
     
 end

 
end