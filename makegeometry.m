function[COCP]=makegeometry(cpx,cpy,nd)
n1=size(cpx,1);
n2=size(cpy,1);
COCP=zeros(n1*n2,3);
if nd==1
        t=1;
for i=1:n2
        for j=1:n1
                COCP(t,1)=(i-1)*n1+j;
                COCP(t,2)=cpx(j);
                COCP(t,3)=cpy(i);
                t=t+1;
        end
end
else
        t=1;
        for i=1:n1
                
                for j=1:n2
                COCP(t,1)=(i-1)*n2+j;
                COCP(t,2)=cpx(i);
                COCP(t,3)=cpy(j);  
                t=t+1;
                end
        end
end
                        
end