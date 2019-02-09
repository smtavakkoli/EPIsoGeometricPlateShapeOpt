function [KK,R,Re,Rn,Xe,Xn,Ye,Yn,Bmat,detJ1,detJ2]=CalculateStiffness(D,DOF,DPx,...
    DPy,re,Wre,rn,Wrn,Ux,Uy,ne,Px,Py,norx)

    R=cell(ne,1);
    Re=R;Rn=R;
    Rx=R;Ry=R;
    Bmat=cell(ne,1);
    Xe=R;Xn=R;
    Ye=R;Yn=R;
    detJ1=cell(ne,1);
    detJ2=zeros(ne,1);
    KK=cell(ne,1);
    
     for i=1:ne
        n_dofe=DOF*numel(Px{i});
        inde=DPx+mod(i,norx)*logical(mod(i,norx))+norx*logical(~mod(i,norx));
        indn=DPy+ceil(i/norx);
        me=Ux(inde+1)-Ux(inde);mn=Uy(indn+1)-Uy(indn);
        pe=Ux(inde+1)+Ux(inde);pn=Uy(indn+1)+Uy(indn);
        detJ2(i)=0.25*me*mn;
        k=zeros(n_dofe);

        for countx=1:length(re)
            e=0.5*(me*re(countx)+pe);
            wre=Wre(countx);
            for county=1:length(rn)
                    n=0.5*(mn*rn(county)+pn);
                    wrn=Wrn(county);
                    [R{i}{countx,county},Re{i}{countx,county},Rn{i}{countx,county},...
                     Ry{i}{countx,county},Rx{i}{countx,county},...
                     Xe{i}(countx,county),Xn{i}(countx,county),...
                     Ye{i}(countx,county),Yn{i}(countx,county),...
                     detJ1{i}(countx,county)]=...
                    DefineSortNurbs(DPx,Ux,inde,e,DPy,Uy,indn,n,Px{i},Py{i});
                    %% Build strain-displacement matrix
                    [Bmat{i}{countx,county}]=buildB(DPx,DPy,DOF,...
                        R{i}{countx,county},Rx{i}{countx,county},Ry{i}{countx,county});
                    %% Build stiffness for i'th element
                     kk=(Bmat{i}{countx,county})'*D{i}*Bmat{i}{countx,county}...
                         *wre*wrn*detJ1{i}(countx,county)*detJ2(i);
                     k=kk+k;
            end
        end
        k=round(k,1);
        if issymmetric(k)
            KK{i}=k;
        end
     end
end





