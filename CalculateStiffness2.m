function [K]=CalculateStiffness2(B,D,DOF,DP1,DP2,re,Wre,rn,Wrn,detJ2,detJ1)
k=zeros(DOF*(DP1+1)*(DP2+1));
for countx=1:length(re)
    wre=Wre(countx);
    for county=1:length(rn)
        wrn=Wrn(county);
        kk=(B{countx,county})'*(D{countx,county})*(B{countx,county})*wre*wrn*detJ1(countx,county)*detJ2;
        k=kk+k;
    end
end
        K=k;
end





