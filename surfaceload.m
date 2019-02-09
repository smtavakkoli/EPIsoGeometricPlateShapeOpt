function [SL]=surfaceload(sl,SL,Numbering,R,detJ,DOF,re,Wre,rn,Wrn,U1,DP1,U2,DP2,n_dof)
%% Surface Load
if (isempty(sl)~=1)
    SL=zeros(n_dof,1);

nc1=length(U1)-2*DP1-1;

for i=1:size(sl,1)
    
        elementnumber=sl(i,1);
        indn=ceil((elementnumber+DP1)/(nc1+DP1))+DP2;
        inde=mod(elementnumber,nc1)+DP1+~mod(elementnumber,nc1)*nc1;
        me=U1(inde+1)-U1(inde);mn=U2(indn+1)-U2(indn);
        detj=0.25*me*mn;
        Sl=sl(i,2:4);
        index=Numbering{elementnumber};
        index1=index';
        index2=index1(:);
      
        for countx=1:length(re)
            for county=1:length(rn)
                r=DOF*(index2-1)+1;
                RSL=zeros(1,numel(index));
                t=1;
                
                for i1=1:size(index,1)
                    for i2=1:size(index,2)
                RSL(t)=R{countx,county,elementnumber}(i1,i2);
                t=t+1;
                    end
                end
               
                h=(Sl(1)*RSL*Wre(countx)*Wrn(county)*detJ(countx,county,elementnumber)*detj)';
                h1=(Sl(2)*RSL*Wre(countx)*Wrn(county)*detJ(countx,county,elementnumber)*detj)';
                h2=(Sl(3)*RSL*Wre(countx)*Wrn(county)*detJ(countx,county,elementnumber)*detj)';

                SL(r)=SL(r)+h;
                SL(r+1)=SL(r+1)+h1;
                SL(r+2)=SL(r+2)+h2;
                   
            end
        end
end
end