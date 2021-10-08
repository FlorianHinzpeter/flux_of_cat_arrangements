function J = DelocalizedFlux(flux_fun,trials,N1,N2,alpha,load_arrang_path,flux_save_path)
%---------------------------------------------------------------------------------------------

%This function loads the Delocalized arrangements from the path and file, load_arrang_path,
%and computes for each arrangement the flux using flux_fun, the function that computes the fluxes.
%Subsequently, the arrangements are saved back using the path and file name
%specified in flux_save_path.

%---------------------------------------------------------------------------------------------

load(load_arrang_path);
 J = zeros(trials+1,2);
 j = zeros(1,2);
 
 for i = 1:trials

     J(i,:) = @flux_fun(N1,N2,X(2*i-1,1:2:(N1+N2)),X(2*i,1:2:(N1+N2)),X(2*i-1,2:2:(N1+N2)),X(2*i,2:2:(N1+N2)),0.01,alpha,alpha);
     
     j = j+J(i,:);
     i
     J(i,:)
 end

 J(trials+1,:) = j/trials;
 
 save(flux_save_path,'J')
 
end