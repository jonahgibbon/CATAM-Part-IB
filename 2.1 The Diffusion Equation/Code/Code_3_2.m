format long g
T=0.1;
N=10;
C=1/2;

%Calculates error for 8 points
W=U2(1,T,0.5*10^(-18));
Counter=0;
Points=zeros(8,2);
while Counter<=7
    Points(Counter+1,1)=1/(N*2^(Counter));
    V=numU(T,N*2^(Counter),C);
    Points(Counter+1,2)=abs(V(end,1)-W);
    Counter=Counter+1;
end
disp(Points)
digits(14)
latex(vpa(sym(Points)))
Points=log(Points);

%Calculates least squares regression line
SXX=dot(Points(:,1),Points(:,1))-(sum(Points(:,1)))^2/8;
SYY=dot(Points(:,2),Points(:,2))-(sum(Points(:,2)))^2/8;
SXY=dot(Points(:,1),Points(:,2))-sum(Points(:,1))*sum(Points(:,2))/8;
YMean=sum(Points(:,2))/8;
XMean=sum(Points(:,1))/8;
PMCC=SXY/(sqrt(SXX*SYY));

X=linspace(Points(1,1),Points(end,1));
Y=(SXY/SXX).*(X-XMean)+YMean;

%Plots scatter plot
figure
scatter(Points(:,1),Points(:,2))
hold on
plot(X,Y)
legend('Error points','Least Squares Regression Line','Location','northwest')
xlabel('log(\delta X)')
ylabel('log(\epsilon)')
print('Image_3_1_Order_Of_Accuracy','-depsc')

disp(SXY/SXX)
    

function answer = numU(T,N,C)
    DeltaT=C/N^2;
    U=zeros(N+1,T/DeltaT+1);
    U(1,:)=1;
    U(1,1)=0.5;
    ColCount=2;
    while ColCount<=T/DeltaT+1
        RowCount=2;
        while RowCount<=N
            U(RowCount,ColCount)=U(RowCount,ColCount-1)+...
                C*(U(RowCount-1,ColCount-1)-2*U(RowCount,ColCount-1)+...
                U(RowCount+1,ColCount-1));
            RowCount=RowCount+1;
        end
        U(RowCount,ColCount)=U(RowCount,ColCount-1)+2*C*(U(RowCount-1,...
            ColCount-1)-U(RowCount,ColCount-1));
        ColCount=ColCount+1;
    end
    answer=U(:,end);
end

function answer = U2(X,T,epsilon)
    answer=0; 
    t=2*lambertw(exp(1)/(pi*epsilon));
    k=1+ceil((1/(pi))*exp(-1/t)*(sqrt(t/(2*T))));
    n=k;
        while n>=1
        g=-4/((2*n-1)*pi)*exp((-((2*n-1)*pi/2)^2)*T);
        h=sin((2*n-1)*pi*X/2);
        answer=answer+g*h;
        n=n-1;
        end
    answer=answer+1;
end