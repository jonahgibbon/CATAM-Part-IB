format long
p=6;
a=4;
k=0;
data=zeros(13,6);
while k<=12
    h=0.1/(2^k);
    data(k+1,1)=h;
    Matrix=RK4Vector(0,1,0,1,h,p,a);
    data(k+1,2)=Matrix(end,2);
    data(k+1,3)=Matrix(end,3);
    data(k+1,4)=y(1,p);
    data(k+1,5)=Matrix(end,2)-y(1,p);
    data(k+1,6)=data(k+1,5)/h^4;
    k=k+1;
end
data
digits(5)
latex(sym(vpa(data)))

function Solution = RK4Vector(x_0,x_n,y_0,z_0,h,p,a)
X=zeros(ceil(x_n/h)+1,1);
Y=zeros(ceil(x_n/h)+1,2);
Y(1,:)=[y_0,z_0];
X(1)=x_0;
counter=2;

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

function out = F(x,Y,p,a)
    out=zeros(1,2);
    out(1)=Y(2);
    out(2)=-p^2*(1+x)^(-a)*Y(1);
end

function output = y(x,p)
output= (1/p)*(1+x)*sin(p-p/(1+x));
end