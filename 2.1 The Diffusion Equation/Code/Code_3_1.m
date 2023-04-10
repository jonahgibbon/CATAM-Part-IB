format long g
T=0.05;
N=10;
C=0.5;

Counter=1;
Data=zeros(N+1,4);
Data(:,2)=numU(T,N,C);
while Counter<=N+1
    Data(Counter,1)=(Counter-1)/N;
    Data(Counter,3)=U2(Data(Counter,1),T,5*10^(-18));
    Data(Counter,4)=Data(Counter,2)-Data(Counter,3);
    Counter=Counter+1;
end

disp(Data)
latex(sym(vpa(Data)))

figure
plot(Data(:,1),Data(:,2))
hold on
plot(Data(:,1),Data(:,3))
legend('Numerical Solution','Analytic Solution','Location','northeast')
xlabel('X')
ylabel('U')
print('Image_3_1','-depsc')

%This calculates the numerical approximation to U at time T
function answer = numU(T,N,C)
    DeltaT=C/N^2;
    U=zeros(N+1,floor(T/DeltaT+1));
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

%This calculates the exact solution
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