format long
LowerP=0;
UpperP=3;
epsilon=7.3*10^(-7);
%This records all the iterates
iterates=[0 LowerP g(LowerP)];
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
    iterates=[iterates;iterates(end,1)+1 P g(P)];
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

disp(iterates)
digits(8)
latex(sym(vpa(iterates)))

function answer=g(x)
    answer=x^2-3;
end