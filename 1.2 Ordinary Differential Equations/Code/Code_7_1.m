data=RK4Vector(0,pi/2,0,1,pi/100,1,0)
latex(sym(vpa(data)))
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