clear;

load('N=120_alpha_1.000000e-01.mat','J');
Jarray(:,:,1)=J;
load('N=120_alpha_1.778279e-01.mat','J');
Jarray(:,:,2)=J;
load('N=120_alpha_3.162278e-01.mat','J');
Jarray(:,:,3)=J;
load('N=120_alpha_5.623413e-01.mat','J');
Jarray(:,:,4)=J;
load('N=120_alpha_1.mat','J');
Jarray(:,:,5)=J;
load('N=120_alpha_1.778279e+00.mat','J');
Jarray(:,:,6)=J;
load('N=120_alpha_3.162278e+00.mat','J');
Jarray(:,:,7)=J;
load('N=120_alpha_5.623413e+00.mat','J');
Jarray(:,:,8)=J;
load('N=120_alpha_10.mat','J');
Jarray(:,:,9)=J;
load('N=120_alpha_1.778279e+01.mat','J');
Jarray(:,:,10)=J;
load('N=120_alpha_3.162278e+01.mat','J');
Jarray(:,:,11)=J;
load('N=120_alpha_5.623413e+01.mat','J');
Jarray(:,:,12)=J;
load('N=120_alpha_100.mat','J');
Jarray(:,:,13)=J;
load('N=120_alpha_1.778279e+02.mat','J');
Jarray(:,:,14)=J;
load('N=120_alpha_3.162278e+02.mat','J');
Jarray(:,:,15)=J;
load('N=120_alpha_5.623413e+02.mat','J');
Jarray(:,:,16)=J;

%%
load('Delocalized_N=120.mat');
%radial distance1
E1r(1:(size(X,1)/2),1:(size(X,2)/2))=sqrt(X(1:2:size(X,1),1:2:size(X,2)).^2 ...
                                         +X(2:2:size(X,1),1:2:size(X,2)).^2);
E2r(1:(size(X,1)/2),1:(size(X,2)/2))=sqrt(X(1:2:size(X,1),2:2:size(X,2)).^2 ...
                                         +X(2:2:size(X,1),2:2:size(X,2)).^2);
%minimum distance
for i=1:(size(X,1)/2)
    for n=1:2
        for j=1:(size(X,2)/2)
            for k=1:size(X,2)
                d(i,n,j,k)=sqrt((X(2*(i-1)+1,(j-1)*2+n)-X(2*(i-1)+1,k)).^2 ...
                               +(X(2*(i-1)+2,(j-1)*2+n)-X(2*(i-1)+2,k)).^2);
            end
        end
    end
end
d(~d) = nan;
%%
for i=1:2
    for j=1:2
        mnE(1:size(d,1),1:size(d,3),i,j)=mean(d(:,i,:,j:2:size(d,4)),4,'omitnan');
        minE(1:size(d,1),1:size(d,3),i,j)=min(d(:,i,:,j:2:size(d,4)),[],4);
    end
end
%%
M=zeros(16,37);
N=M;
z=zeros(1,size(J,1)-1);
for i=1:16
    M(i,1)=0.1*(10^0.25).^(i-1);
    M(i,2)=mean(Jarray(1:(size(J,1)-1),2,i),1);
    M(i,3)=std(Jarray(1:(size(J,1)-1),2,i),1);
    M(i,4)= M(i,3)./M(i,2);
    [M(i,5),M(i,6),M(i,7)]=regression(mean(E1r,2)',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,8),M(i,9),M(i,10)]=regression(mean(E2r,2)',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,11),M(i,12),M(i,13)]=regression((mean(E2r,2)-mean(E1r,2))',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,14),M(i,15),M(i,16)]=regression(mean(minE(:,:,1,1),2)',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,17),M(i,18),M(i,19)]=regression(mean(minE(:,:,1,2),2)',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,20),M(i,21),M(i,22)]=regression(mean(minE(:,:,2,1),2)',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,23),M(i,24),M(i,25)]=regression(mean(minE(:,:,2,2),2)',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,26),M(i,27),M(i,28)]=regression(mean(mnE(:,:,1,1),2)',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,29),M(i,30),M(i,31)]=regression(mean(mnE(:,:,1,2),2)',Jarray(1:(size(J,1)-1),2,i)');
    [M(i,32),M(i,33),M(i,34)]=regression(mean(mnE(:,:,2,2),2)',Jarray(1:(size(J,1)-1),2,i)');
    
    N(i,1)=0.1*(10^0.25).^(i-1);
    N(i,2)=mean(Jarray(1:(size(J,1)-1),2,i),1);
    N(i,3)=std(Jarray(1:(size(J,1)-1),2,i),1);
    N(i,4)=N(i,3)./N(i,2);
    lm=fitlm(mean(E1r,2),Jarray(1:(size(J,1)-1),2,i));
    N(i,5:7)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(E2r,2),Jarray(1:(size(J,1)-1),2,i));
    N(i,8:10)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(E2r,2)-mean(E1r,2),Jarray(1:(size(J,1)-1),2,i));
    N(i,11:13)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(minE(:,:,1,1),2),Jarray(1:(size(J,1)-1),2,i));
    N(i,14:16)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(minE(:,:,1,2),2),Jarray(1:(size(J,1)-1),2,i));
    N(i,17:19)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(minE(:,:,2,1),2),Jarray(1:(size(J,1)-1),2,i));
    N(i,20:22)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(minE(:,:,2,2),2),Jarray(1:(size(J,1)-1),2,i));
    N(i,23:25)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(mnE(:,:,1,1),2),Jarray(1:(size(J,1)-1),2,i));
    N(i,26:28)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(mnE(:,:,1,2),2),Jarray(1:(size(J,1)-1),2,i));
    N(i,29:31)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
    lm=fitlm(mean(mnE(:,:,2,2),2),Jarray(1:(size(J,1)-1),2,i));
    N(i,32:34)=[lm.Rsquared.Ordinary,lm.Coefficients{2,1},lm.Coefficients{1,1}];
end
%%
close all
%{
figure
for i=1:16
    subplot(4,4,i);
    plot(mean(E1r,2),Jarray(1:(size(Jarray,1)-1),2,i),'.b',...
         mean(E2r,2),Jarray(1:(size(Jarray,1)-1),2,i),'.r',...
         0.5:0.1:0.8,N(i,6)*(0.5:0.1:0.8)+N(i,7),'--c',...
         0.5:0.1:0.8,M(i,6)*(0.5:0.1:0.8)+M(i,7),'--k',...
         0.5:0.1:0.8,N(i,9)*(0.5:0.1:0.8)+N(i,10),':c',...
         0.5:0.1:0.8,M(i,9)*(0.5:0.1:0.8)+M(i,10),':k');
end

figure
for i=1:16
    subplot(4,4,i);
    plot(mean(minE(:,:,1,1),2),Jarray(1:(size(Jarray,1)-1),2,i),'.b',...
         mean(minE(:,:,2,2),2),Jarray(1:(size(Jarray,1)-1),2,i),'.r',...
         0:0.1:0.4,N(i,15)*(0:0.1:0.4)+N(i,16),'--c',...
         0:0.1:0.4,M(i,15)*(0:0.1:0.4)+M(i,16),'--k',...
         0:0.1:0.4,N(i,24)*(0:0.1:0.4)+N(i,25),':c',...
         0:0.1:0.4,M(i,24)*(0:0.1:0.4)+M(i,25),':k');
end


figure
for i=1:16
    subplot(4,4,i);
    plot(mean(minE(:,:,1,2),2),Jarray(1:(size(Jarray,1)-1),2,i),'.b',...
         mean(minE(:,:,2,1),2),Jarray(1:(size(Jarray,1)-1),2,i),'.r',...
         0:0.1:0.4,N(i,18)*(0:0.1:0.4)+N(i,19),'--c',...
         0:0.1:0.4,M(i,18)*(0:0.1:0.4)+M(i,19),'--k',...
         0:0.1:0.4,N(i,21)*(0:0.1:0.4)+N(i,22),':c',...
         0:0.1:0.4,M(i,21)*(0:0.1:0.4)+M(i,22),':k');
end


figure
for i=1:16
    subplot(4,4,i);
    plot(mean(mnE(:,:,1,1),2),Jarray(1:(size(Jarray,1)-1),2,i),'.b',...
         mean(mnE(:,:,2,2),2),Jarray(1:(size(Jarray,1)-1),2,i),'.r',...
         0.7:0.1:1.1,N(i,27)*(0.7:0.1:1.1)+N(i,28),'--c',...
         0.7:0.1:1.1,M(i,27)*(0.7:0.1:1.1)+M(i,28),'--k',...
         0.7:0.1:1.1,N(i,33)*(0.7:0.1:1.1)+N(i,34),':c',...
         0.7:0.1:1.1,M(i,33)*(0.7:0.1:1.1)+M(i,34),':k');
end

figure
for i=1:16
    subplot(4,4,i);
    plot(mean(mnE(:,:,1,2),2),Jarray(1:(size(Jarray,1)-1),2,i),'.b',...
        0.7:0.1:1.1,N(i,30)*(0.7:0.1:1.1)+N(i,31),'--c',...
        0.7:0.1:1.1,M(i,30)*(0.7:0.1:1.1)+M(i,31),'--k');
end

%}
%{
figure
        plot(mean(minE(:,:,1,2),2),Jarray(1:(size(Jarray,1)-1),2,9),'.b',...
        0:0.1:0.2,M(9,18)*(0:0.1:0.2)+M(9,19),'--k');
        xlabel('m_{12}')
        ylabel('J_P','rot',-360)
        
        figure
        plot(mean(mnE(:,:,1,2),2),Jarray(1:(size(Jarray,1)-1),2,5),'.b',...
        0.7:0.1:1.1,M(5,30)*(0.7:0.1:1.1)+M(5,31),'--k');
        xlabel('d_{12}')
        ylabel('J_P','rot',-360)
        
         figure
        plot(mean(mnE(:,:,1,2),2),Jarray(1:(size(Jarray,1)-1),2,13),'.b',...
        0.7:0.1:1.1,M(13,30)*(0.7:0.1:1.1)+M(13,31),'--k');
        xlabel('d_{12}')
        ylabel('J_P','rot',-360)
%}   
%semilogx(N(:,1),N(:,5:3:34))
%semilogx(N(:,1),N(:,17),'k',N(:,1),N(:,29),'r')
semilogx(N(:,1),N(:,5),'k',N(:,1),N(:,8),'g',N(:,1),N(:,26),'b',N(:,1),N(:,32),'r')
N(:,1)
N(:,5)
N(:,8)
N(:,26)
N(:,32)