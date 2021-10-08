function J = FusedFlux(number,trials,N1,N2,load_arrang_path,flux_save_path)
%---------------------------------------------------------------------------------------------

%This function loads the pair arrangements from the path, load_arrang_path,
%and computes for each arrangement the flux using Flux_bd_loss.
%Subsequently, the fluxes are saved back using the path and file name
%specified in flux_save_path which is an array of several paths for each instance in Alpha.
%this is iterated for the alpha parameters indicated in Alpha.Alpha.

%---------------------------------------------------------------------------------------------

load(load_arrang_path);
Alpha = [0.1 1 10 100 1000];

 J = zeros(trials,2*length(Alpha));
 Npair = N1;
 
 for i = 1:trials

     for j = 1:length(Alpha)
     
          J(i,2*j-1:2*j) = Flux_bd_loss(N1,N2,X(3*i-2,1:Npair),X(3*i-1,1:Npair),X(3*i,1:Npair),X(3*i-2,Npair+1:2*Npair),X(3*i-1,Npair+1:2*Npair),X(3*i,Npair+1:2*Npair),0.01,Alpha(j),Alpha(j));
     
     save(flux_save_path(j),'J')
     
     end
 
     i
 
 end
 
 J
 
 save(flux_save_path(length(Alpha)+1),'J')
 
end