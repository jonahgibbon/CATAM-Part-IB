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
ker=Ker(m,p,rank)

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

function basis = Ker(m,p,rank)
%This creates a blank array to be filled with vectors for the basis
basis=zeros(size(m,2),size(m,2)-rank);

%This resets the rows,columns and noBasis (which represents the number of
%basis vectors created) counters.
rows=size(m,1);
columns=1;
noBasis=1;
lastColumn=size(m,2)+1;

%This iterates through until it finds a non-zero element
while rows >= 1
    while columns <= size(m,2)
        if m(rows,columns)~=0
            %If there is a 1 in the far right of the matrix, this ensures
            %the last element of the vector is 0
            if columns==size(m,2)
                lastColumn=columns;
            else
                %This iterates between the column where the non-zero
                %element is, and the last column that caused a vector to be
                %produced. It creates new linearly-independant vectors for
                %the basis, before ensuring that the vector lies in the
                %kernal
                thisColumn=columns;
                columns=columns+1;
                while columns<lastColumn
                    basis(columns,noBasis)=1;
                    noBasis=noBasis+1;
                    columns=columns+1;
                end
                basis(thisColumn,:)=-m(rows,:)*basis;
                lastColumn=thisColumn;
            end
            break
        end
        columns=columns+1;
    end
    columns=1;
    rows=rows-1;
end
%This ensures that the vector lies within mod(p)
basis=mod(basis,p);
end