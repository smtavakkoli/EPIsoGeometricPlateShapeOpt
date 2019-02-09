function [dkN]=DerivativeNurbs(P,U,ind,range,k,u)
%%
dkN=zeros(P+1,1);
N=BuildNurbs(P-k,U,ind-P,ind+k,range,u);
for z=0:P
A=zeros(k+1);
if P>=k

indn=ind-P;

A(1,1)=1;
for i=2:k+1
        a=U(ind-i+2);b=U(indn);
A(i,1)=A(i-1,1)/(a-b);
if (isinf(A(i,1))==1 || isnan(A(i,1))==1)
                A(i,1)=0;
end
end
for i=2:k+1
for j=2:i
        if i==j;
                a=U(ind+1);b=U(indn+i-1);
        A(i,j)=-A(i-1,i-1)/(a-b);
        else
           a=U(ind+j-i+1);b=U(indn+j-1);     
        A(i,j)=(A(i-1,j)-A(i-1,j-1))/(a-b);
        end
        if (isinf(A(i,j))==1 || isnan(A(i,j))==1)
                A(i,j)=0;
        end
end
end

dkN(z+1,1)=(factorial(P)/factorial(P-k))*A(k+1,:)*N(1+z:k+z+1); 
ind=ind+1;
end
end