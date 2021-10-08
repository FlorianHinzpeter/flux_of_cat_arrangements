function Histogramm

a = -1:0.25:2.75;

Alpha = 10.^a;
N = 60;

j = zeros(5000,length(Alpha));

for idx0 = 1:length(Alpha)

load(sprintf('N=%i_alpha_%i.mat',N,Alpha(idx0)));

j(:,idx0) = J(1:end-1,2)./median(J(1:end-1,2));

end



[count,bins] = hist(j,30)

Count = zeros(length(count(:,1)),length(count(1,:))+3*(length(count(1,:))-1));

for idx1 = 1:length(count(1,:))
    
   Count(:,4*idx1-3) = count(:,idx1);
    
end

bar3(bins,Count/3000,'hist')
set(gca,'XTickLabel',{'10^{-1}';[];[];[] ;'10^0';[];[];[] ;'10^1';[];[];[] ;'10^2';[];[];[] })

end

    
    
    
    