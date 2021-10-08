function X = PairFusedArrangements(N,trials,savepath)
%---------------------------------------------------------------------------------------------

%This function generates an ensemble of Pair arrangements.
%The number of samples in the ensemble is specified with the parameter trials
%The number of particles in the arrangemets is N.
%The arrangemets are saved back under the path and file name indicated by savepath.

%---------------------------------------------------------------------------------------------

rng('shuffle');

X = zeros(2*trials,N);

for i = 1:trials
  
    constr = 1;
    
    while constr == 1
    
    x = Randomarrangement_pairs(0.01,N,2*0.01+0.01);
    
    constr = constraint(x(1,:),x(2,:),1,0.01,N);
    
    end
    
    X(2*i-1:2*i,:) = x;
    i
end

save(savepath,'X')

end