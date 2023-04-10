%Initial data
p=7121;
a=1;
b=1;
tic
inverses=zeros(p-1,1);
while a<p
    while 1
        %This checks if b is the inverse of a, and if so, stores it then 
        %breaks the loop.
        if mod(a*b,p)==1
            inverses(a)=b;
            b=1;
            break
        end
        b=b+1;
    end
    a=a+1;
end
disp(inverses)
toc