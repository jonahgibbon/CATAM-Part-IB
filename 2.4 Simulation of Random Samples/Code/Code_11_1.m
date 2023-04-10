format shortg

%Set parameters
n=100;
mu=0;
SigmaSquared=1;
alpha=0.8;
SampleSize=25;

Counter=1;
Data=zeros(SampleSize,4);

disp([num2str(alpha),'% confidence interval is XBar plus or minus ',...
    num2str(norminv((1+alpha)/2)*sqrt(SigmaSquared/n))])

%Calculates random variables and fills table
while Counter<=SampleSize
    Variables=NormGenerate(n,mu,SigmaSquared);
    
    Data(Counter,1)=sum(Variables)/n;
    Data(Counter,2)=sum(Variables)/n-norminv((1+alpha)/2)*sqrt(SigmaSquared/n);
    Data(Counter,3)=sum(Variables)/n+norminv((1+alpha)/2)*sqrt(SigmaSquared/n);
    Data(Counter,4)=and(mu>=Data(Counter,2),mu<=Data(Counter,3));
    Counter=Counter+1;
end
disp(Data)
% latex(sym(vpa(Data)))

Out=SampleSize-sum(Data(:,4));
disp(['Number of cases outside interval: ',num2str(Out)])

function Variables = NormGenerate(n,mu,SigmaSquared)
    A=rand((n-mod(n,2))/2+mod(n,2),1);
    B=rand((n-mod(n,2))/2+mod(n,2),1);
    Phi=2*pi*A;
    V=-2*log(1-B);
    X=mu+SigmaSquared*sqrt(V).*cos(Phi);
    Y=mu+SigmaSquared*sqrt(V).*sin(Phi);
    Variables=[X;Y(1:n-(n-mod(n,2))/2-mod(n,2),1)];
end