function [X,Y] = CP2Phsical(n_p_e,Numbering,NOR1,U1,U2,P1,P2,n1,n2)

index=cumsum(n_p_e);
n_e=index(end);
X=cell(n_e,1);
Y=cell(n_e,1);

for j=1:n_e
    
        patch_index=find(index>=j,1,'first');
        DP1=size(Numbering{j},2)-1;
        DP2=size(Numbering{j},1)-1;
        if patch_index>1
            ind=j-sum(n_p_e(patch_index-1:-1:1));
        else
            ind=j;
        end
        inde=DP1+mod(ind,NOR1(patch_index))*logical(mod(ind,NOR1(patch_index)))...
        +NOR1(patch_index)*logical(~mod(ind,NOR1(patch_index)));
        indn=DP2+ceil(ind/NOR1(patch_index));
        me=U1{patch_index}(inde+1)-U1{patch_index}(inde);
        mn=U2{patch_index}(indn+1)-U2{patch_index}(indn);
        

    for i1=0:n1

        e=U1{patch_index}(inde)+i1*me/n1;
        N=nurbs2(DP1,U1{patch_index},inde-DP1,inde,inde,e);
        for i2=0:n2
                n=U2{patch_index}(indn)+i2*mn/n2;
                M=nurbs2(DP2,U2{patch_index},indn-DP2,indn,indn,n);    
                Px=P1{j};
                Py=P2{j};
                X{j}(i1+1,i2+1)=M'*Px*N;
                Y{j}(i1+1,i2+1)=M'*Py*N;    
        end
    end
end
% X=cell2mat(X);
% Y=cell2mat(Y);
% X(abs(X)<1e-5)=0;
% Y(abs(Y)<1e-5)=0;
end
