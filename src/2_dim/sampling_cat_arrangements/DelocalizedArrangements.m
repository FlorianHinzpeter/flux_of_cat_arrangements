function X = ArrangementDelocalized(N,trials,savepath)
%---------------------------------------------------------------------------------------------

%This function generates an ensemble of Delocalized arrangements.
%The number of samples in the ensemble is specified with the parameter trials
%The number of particles in the arrangemets is N.
%The arrangemets are saved back under the path and file name indicated by savepath.

%---------------------------------------------------------------------------------------------

rng('shuffle') 

X = zeros(2*trials,N);

for i = 1:trials
  
   X(2*i-1:2*i,:) = Randomarrangement_dynamic(1,0.01,N,3); 
    
end

save(savepath,'X');

end