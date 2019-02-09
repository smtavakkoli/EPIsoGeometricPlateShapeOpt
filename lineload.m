function [LL]=lineload(ll,LL,DOF,DP1,DP2,U1,U2,Wre,Wrn,re,rn,NOR,Numbering)

for i=1:size(ll,1)


if ll(i,1)==0

    ind1=ll(i,2)-DP1;
    indn=ll(i,3)-DP2;
    
    j=(indn-1)*NOR+ind1;
    switch ind1

            
        case length(U1(DP1+1:end-DP1))
            
            ind1=ind1-1;
            j=(indn-1)*NOR+ind1;
    end
    
    N=nurbs2(DP1,U1,ind1,ind1+DP1,ind1+DP1,U1(ll(i,2)));
    ind=Numbering{j};
    n1=U2(ll(i,3));
    n2=U2(ll(i,3)+1);
    mn=n2-n1;pn=n2+n1;
      
for count=1:length(rn)
     
     n=0.5*(mn*rn(count)+pn);
      M=nurbs2(DP2,U2,indn,indn+DP2,indn+DP2,n);
     
for j=1:size(ind,1)     
     
     
  LL(DOF*(ind(j,:)-1)+1)=LL(DOF*(ind(j,:)-1)+1)...
  +M(j)*N*ll(i,4)*Wrn(count)*0.5*mn;

  LL(DOF*(ind(j,:)-1)+2)=LL(DOF*(ind(j,:)-1)+2)...
  +M(j)*N*ll(i,5)*Wrn(count)*0.5*mn;

   LL(DOF*ind(j,:))=LL(DOF*ind(j,:))...
  +M(j)*N*ll(i,6)*Wrn(count)*0.5*mn;

end
 end
else
    
    ind1=ll(i,2)-DP2;
    inde=ll(i,3)-DP1;
    j=(ind1-1)*NOR+inde;
    
    switch ind1

            
        case length(U2(DP2+1:end-DP2))
            
            ind1=ind1-1;
            j=(ind1-1)*NOR+inde;
    end
    
    M=nurbs2(DP2,U2,ind1,ind1+DP2,ind1+DP2,U2(ll(i,2)));
    ind=Numbering{j};
    e1=U1(ll(i,3));
    e2=U1(ll(i,3)+1);
    me=e2-e1;pe=e2+e1;
    
    for count=1:length(re)
     
     e=0.5*(me*re(count)+pe);
     N=nurbs2(DP1,U1,inde,inde+DP1,inde+DP1,e);
     
for j=1:size(ind,1)     
    
     
     
  LL(DOF*(ind(j,:)-1)+1)=LL(DOF*(ind(j,:)-1)+1)...
  +N(j)*M*ll(i,4)*Wre(count)*0.5*me;

  LL(DOF*(ind(j,:)-1)+2)=LL(DOF*(ind(j,:)-1)+2)...
  +N(j)*M*ll(i,5)*Wre(count)*0.5*me;

   LL(DOF*ind(j,:))=LL(DOF*ind(j,:))...
  +N(j)*M*ll(i,6)*Wre(count)*0.5*me;

end
    end
 
end
end
end
