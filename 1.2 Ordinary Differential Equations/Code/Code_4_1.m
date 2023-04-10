k=0;
Error=zeros(16,3);
while k<=15
    h=1.6/(2^k);
    data=Euler(0,1.6,0,h);
    Error(k+1,:)=[h,data(end,4),0];
    data=RK4(0,1.6,0,h);
    Error(k+1,3)=data(end,4);
    k=k+1;
end
digits(5)
disp(Error)
latex(sym(vpa(Error)))
logError=log(abs(Error));
meanError=mean(logError);
scatter(logError(:,1),logError(:,2))
hold on
x=xlim;
plot(linspace(x(1),x(2)),linspace(x(1),x(2))+meanError(2)-...
    meanError(1))
% Adds labels and saves the image
legend({'Data','y=x+a'},'location',...
    'northwest')
title('Euler Error as h decreases')
xlabel('log(h)')
ylabel('log(|E_{n}|)')
print('Image_4_1','-depsc')
figure
scatter(logError(:,1),logError(:,3))
hold on
x=xlim;
plot(linspace(x(1),x(2)),4*linspace(x(1),x(2))+meanError(3)-...
    4*meanError(1))
%Adds labels and saves the image
legend({'Data','y=4x+b'},'location',...
    'northwest')
title('Runge-Kutta Error as h decreases')
xlabel('log(h)')
ylabel('log(|E_{n}|)')
print('Image_4_2','-depsc')



function data = Euler(x_0,x_n,y_0,h)
%Creates a table to be filled
    data=zeros(ceil((x_n-x_0)/h)+1,5);
    data(1,:)=[x_0,y_0,y(x_0),y_0-y(x_0),0];
    counter=2;
%Iterates through the rows filling them
    while counter<=ceil((x_n-x_0)/h)+1
        data(counter,:)=[(counter-1)*h+x_0,...
            data(counter-1,2)+h*f((counter-2)*h+x_0,data(counter-1,2)),...
            y((counter-1)*h+x_0),...
            data(counter-1,2)+h*f((counter-2)*h+x_0,data(counter-1,2))-...
            y((counter-1)*h+x_0),...
            (data(counter-1,2)+h*f((counter-2)*h+x_0,data(counter-1,2))-...
            y((counter-1)*h+x_0))/data(counter-1,4)];
        counter=counter+1;
    end
end

function data = RK4(x_0,x_n,y_0,h)
%Follows a similar idea to the Euler function
    data=zeros(ceil((x_n-x_0)/h)+1,5);
    data(1,:)=[x_0,y_0,y(x_0),y_0-y(x_0),0];
    counter=2;
    while counter<=ceil((x_n-x_0)/h)+1
        k1=h*f((counter-2)*h+x_0,data(counter-1,2));
        k2=h*f((counter-2+1/2)*h+x_0,data(counter-1,2)+1/2*k1);
        k3=h*f((counter-2+1/2)*h+x_0,data(counter-1,2)+1/2*k2);
        k4=h*f((counter-2+1)*h+x_0,data(counter-1,2)+k3);
        data(counter,:)=[(counter-1)*h+x_0,...
            data(counter-1,2)+1/6*(k1+2*k2+2*k3+k4),...
            y((counter-1)*h+x_0),...
            data(counter-1,2)+1/6*(k1+2*k2+2*k3+k4)-y((counter-1)*h+x_0),...
            (data(counter-1,2)+1/6*(k1+2*k2+2*k3+k4)-y((counter-1)*h+x_0))/...
            data(counter-1,4)];
        counter=counter+1;
    end   
end

%These two functions were created to stop repeatedly defining them in the
%above
function z = f(x,y)
    z = -4*y+4*exp(-2*x);
end

function z = y(x)
    z= -2*exp(-4*x)+2*exp(-2*x);
end