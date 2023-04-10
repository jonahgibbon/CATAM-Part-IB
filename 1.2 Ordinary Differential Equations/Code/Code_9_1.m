format long
%This searches for the first intervals where the eigenvalues lie
InitialStep=1;
Eigenvalues=zeros(5,1);
NoVals=1;
Counter=1;
InitialSign=sign(g(InitialStep));
while NoVals<6
    x=Counter*InitialStep;
    if sign(g(x))~=InitialSign
        Eigenvalues(NoVals,1)=x;
        NoVals=NoVals+1;
        InitialSign=sign(g(x));
    end
    Counter=Counter+1;
end
disp('First guess of eigenvalues')
disp(Eigenvalues)

%This clarifies that we can use the particular epsilon and h specified in
%the report
NoVals=1;
iterates=zeros(5,1);
gradient=zeros(5,1);
while NoVals<6
    iterates(NoVals)=gError(Eigenvalues(NoVals));
    gradient(NoVals)=g(Eigenvalues(NoVals))-...
        g(Eigenvalues(NoVals)-1);
    NoVals=NoVals+1;
end
disp('This clarifies that -1<g(p)<1 when h=0.1 is used')
disp(iterates)
disp("This calculates when |g'(p)| is smallest")
disp(abs(gradient))
epsilon=abs(g(Eigenvalues(5,1))-g(Eigenvalues(5,1)-1))*5*10^(-6)

%This calculates the accurate eigenvalues using the particular epsilon and
%h
NoVals=1;
while NoVals<6
    LowerP=Eigenvalues(NoVals,1)-1;
    UpperP=Eigenvalues(NoVals,1);
    if g(LowerP)==g(UpperP)
        error('The g value must differ at these bounds')
    end
    if g(LowerP)<g(UpperP)
        gradient=+1;
    else
        gradient=-1;
    end
    while true
        P=(g(UpperP)*LowerP-g(LowerP)*UpperP)/(g(UpperP)-g(LowerP));
        if abs(g(P))<epsilon
            break 
        else
            %This checks what the gradient is at the root
            if g(LowerP)<g(UpperP)
                gradient=+1;
            else
                gradient=-1;
            end
            %This determines which bound should be replaced
            if g(P)*gradient<0
                LowerP=P;
            else
                UpperP=P;
            end
        end
    end
    Eigenvalues(NoVals,1)=P;
    NoVals=NoVals+1;
end
disp('Eigenvalues:')
disp(Eigenvalues)
digits(16)
latex(sym(vpa(Eigenvalues)))

%This plots the normalised eigenfunctions
NoVals=1;
while NoVals<6
    m=RK4Vector(0,1,0,1,0.001,Eigenvalues(NoVals),8);
    A=trapz(m(:,2).^2./((1+m(:,1)).^(8)))*Eigenvalues(NoVals)^2;
    m(:,2)=(1/sqrt(A)).*m(:,2);
    figure
    plot(m(:,1),m(:,2))
    xlabel('x')
    ylabel('y')
    print(strcat('Image_9_',num2str(NoVals)),'-depsc')
    NoVals=NoVals+1;
end

function answer=gError(x)
    answer=RK4Vector(0,1,0,1,0.1,x,8);
    answer=answer(end,2);
end

function answer=g(x)
    answer=RK4Vector(0,1,0,1,0.001,x,8);
    answer=answer(end,2);
end

function Solution = RK4Vector(x_0,x_n,y_0,z_0,h,p,a)
%Creates the table of iterates
X=zeros(ceil(x_n/h)+1,1);
Y=zeros(ceil(x_n/h)+1,2);
Y(1,:)=[y_0,z_0];
X(1)=x_0;
counter=2;
%This fills the table using the Runge-Kutta method
while counter<=ceil((x_n-x_0)/h)+1
    X(counter)=(counter-1)*h;
    K1=h*F(X(counter-1),Y(counter-1,:),p,a);
    K2=h*F(X(counter-1)+h*1/2,Y(counter-1,:)+K1/2,p,a);
    K3=h*F(X(counter-1)+h*1/2,Y(counter-1,:)+K2/2,p,a);
    K4=h*F(X(counter-1)+h,Y(counter-1,:)+K3,p,a);
    Y(counter,:)=Y(counter-1,:)+(K1+2*K2+2*K3+K4)/6;
    counter=counter+1;
end
Solution=[X,Y];
end
%This is f(X,Y)
function out = F(x,Y,p,a)
    out=zeros(1,2);
    %f_1
    out(1)=Y(2);
    %f_2
    out(2)=-p^2*(1+x)^(-a)*Y(1);
end