function X = Cluster_Arrangements(trials,N,r0,savepath)
%---------------------------------------------------------------------------------------------

%This function generates an ensemble of cluster arrangements with cluster radius r0.
%The number of samples in the ensemble is specified with the parameter trials
%The number of particles in the arrangemets is N.
%The arrangemets are saved back under the path and file name indicated by savepath.

%---------------------------------------------------------------------------------------------

X = zeros(2*trials,N);

for i=1:trials
    
   X(2*i-1:2*i,:) = Randomarrangement_clusters(r0,0.01,N,3);
    
end

save(savepath,'X')

end

