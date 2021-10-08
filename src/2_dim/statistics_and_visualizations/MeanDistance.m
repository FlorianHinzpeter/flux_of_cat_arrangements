function MeanDist = MeanDistance(arrangement_path)

load(arrangement_path);

KD = zeros(size(X(1,1:end-1)));
MinDist = 0;
n = 0;

MeanDist = zeros(size(X(1,:)));

for idx0 = 1:length(X(:,1))

for idx1 = 1:length(X(1,:))
    
    n = n+1;
    
    for idx2 = [1:idx1-1,idx1+1:length(X(1,:))]
        
        KD(idx2) = sqrt((X(1,idx1)-X(1,idx2))^2+(X(2,idx1)-X(2,idx2))^2);
        
    end
    
     MinDist = MinDist+min(KD);
    
end

MeanDist(idx0) = MinDist/n;
n
MinDist = 0;
n= 0;

end



