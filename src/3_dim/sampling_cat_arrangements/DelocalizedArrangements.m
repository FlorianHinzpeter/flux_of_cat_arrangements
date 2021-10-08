function X = DelocalizedArrangements(number,trials,N,savepath)
%---------------------------------------------------------------------------------------------

%This function generates an ensemble of Delocalized arrangements.
%The number of samples in the ensemble is specified with the parameter trials
%The number of particles in the arrangemets is N.
%The arrangemets are saved back under the path and file name indicated by savepath.

%---------------------------------------------------------------------------------------------

rng('shuffle')

X = zeros(3*trials,N);

r0 = 1;

for i=1:trials
    
   X(3*i-2:3*i,:) = Randomarrangement_dynamic(r0,0.0105,N);
    
end

save(savepath,'X')

end

