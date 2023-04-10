%Initial data
m=[4,6,5,2,3,1;5,0,3,0,1,0;1,5,7,1,0,12;5,5,0,3,1,7;2,1,2,4,0,5];
p=2;

a=1;
b=1;
%This creates a p-1x1 zero vector to be filled
inverses=zeros(p-1,1);
while a<p
    while 1
        %This checks if b is the inverse of a, and if so, stores b as the
        %ath entry in inverses, then resets b=1 and breaks out the loop.
        if mod(a*b,p)==1
            inverses(a)=b;
            b=1;
            break
        end
        b=b+1;
    end
    a=a+1;
end

[m,rank]=Echelon(m,p,inverses)

function [m,rank] = Echelon(m,p,inverses)
%First this makes sure all elements of the matrix lie within mod p
m=mod(m,p);
%The dimensions of the matrix are used for the limits of the while loop
dimensions=size(m);
columns=1;
rows=1;
%lowerRows tracks how many rows have already changed to echelon form
rank=0;
while columns<=dimensions(2)
    %For each column, we search downwards until we find a non-zero element.
    %If we are successful, we convert this column and the top most possible 
    %row into echelon form, before converting the 'lower-right' matrix into
    %echelon form. Otherwise we move onto the next column.
    while rows<=dimensions(1)
        if m(rows,columns)~=0
            m=T(m,rows,rank+1);
            m=D(m,rank+1,m(rank+1,columns),inverses,p);
            rows=rank+2;
            while rows<=dimensions(1)
                m=S(m,rows,m(rows,columns),rank+1,p);
                rows=rows+1;
            end
            rank=rank+1;
            break
        end
        rows=rows+1;
    end
    %We continue to convert the 'lower-right' matrix into echelon form
    rows=rank+1;
    columns=columns+1;
end
end

function m = T(m,i,j)
    v=m(i,:);
    m(i,:)=m(j,:);
    m(j,:)=v;
end

%This divides a row i by a by indexing its inverse from 'inverses'
function m = D(m,i,a,inverses,p)
    v=m(i,:);
    v=mod(inverses(a)*v,p);
    m(i,:)=v;
end

%This takes a multiple a of row j and minuses this from row i
function m = S(m,i,a,j,p)
    if i==j
        error("i should not equal j")
    else
        m(i,:)=m(i,:)-a*m(j,:);
    end
    %This makes sure that the matrix still lies in mod(p)
    m=mod(m,p);
end