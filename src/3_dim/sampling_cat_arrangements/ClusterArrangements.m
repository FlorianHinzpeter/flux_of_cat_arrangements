function X = ClusterArrangements(number,trials,N,savepath)
%---------------------------------------------------------------------------------------------

%This function generates an ensemble of cluster arrangements with cluster radius r0.
%The number of samples in the ensemble is specified with the parameter trials
%The number of particles in the arrangemets is N.
%The arrangemets are saved back under the path and file name indicated by savepath.

%---------------------------------------------------------------------------------------------

rng('shuffle')

X = zeros(3*trials,N);

r0 = 0.0737;

for i=1:trials
    
   X(3*i-2:3*i,:) = Randomarrangement_clusters(r0,0.0101,N);
    
end

save(savepath,'X')

end

