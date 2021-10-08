function J = ClusterFlux(flux_fun,trials,N1,N2,alpha,load_arrang_path,flux_save_path)
%---------------------------------------------------------------------------------------------

%This function loads the clustered arrangements from the path, load_arrang_path,
%and computes for each arrangement the flux using the flux_fun, the function that computes the flux.
%Subsequently, the arrangements are saved back using the path and file name
%specified in flux_save_path.

%---------------------------------------------------------------------------------------------

%%load the arrangement from %%
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