function J = PairFlux(flux_fun,trials,N1,N2,load_arrang_path,flux_save_path)
%---------------------------------------------------------------------------------------------

%This function loads the Pair arrangements from the path and file, load_arrang_path,
%and computes for each arrangement the flux using flux_fun, the function that computes the flux.
%Subsequently, the arrangements are saved back using the path and file name
%specified in flux_save_path.

%---------------------------------------------------------------------------------------------

N = N1+N2;
load(load_arrang_path);
 J = zeros(trials,2);

 Alpha = [0.1 100];
 
 for j = 1:2
     
     alpha = Alpha(j);
 
 for i = 1:trials

     J(i,:) = @flux_fun(N1,N2,X(2*i-1,1:(N1+N2)/2),X(2*i,1:(N1+N2)/2),X(2*i-1,(N1+N2)/2+1:(N1+N2)),X(2*i,(N1+N2)/2+1:(N1+N2)),0.01,alpha,alpha);
     
    i
 end

 save(flux_save_path,'J')
 
 end
 
end