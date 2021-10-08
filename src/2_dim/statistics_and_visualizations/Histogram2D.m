function Histogram2D(flux_path)

load(flux_path);

j = J(1:end-1,2)./median(J(1:end-1,2));
 
histogram(j,'BinWidth',0.02,'FaceColor',[0 0.7843 1],'Normalization','pdf')


axis([0.5 1.5 0 7]);

set(gca,'FontSize',50,'FontName','Arial','LineWidth',8.5)


