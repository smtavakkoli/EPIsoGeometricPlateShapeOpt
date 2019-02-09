function [K_Px]= Stiffness_sensitivity_x(wanted,D,B,Re,Rn,Xe,Xn,Ye,Yn,...
    detJ1,detJ2,Numbering,DOF,n_dof,w_guss_x,w_guss_y)
%%

C=[];
for i=1:length(Numbering)
    ind=find(Numbering{i}(:)==wanted);
    if ~isempty(ind)
        C=[C;i];
    end
end
  
C(isempty(C))=[];
 
    %% Calculate K_px{i} and K_py{i}
    Re1=Re(C);Rn1=Rn(C);
    Xe1=Xe(C);Xn1=Xn(C);
    Ye1=Ye(C);Yn1=Yn(C);
    Numbering1=Numbering(C);
    D1=D(C);
    B1=B(C);
    detJ11=detJ1(C);
    detJ2_1=detJ2(C);
    K_px=cell(length(C),1);
%     K_py=cell(length(C),1);
    wgx=w_guss_x(C);
    wgy=w_guss_y(C);
    %%...........................
    for i1=1:length(C)
        
        [r,c]=find(Numbering1{i1}==wanted);
        k_px=zeros(DOF*numel(Numbering1{i1}));
%         k_py=zeros(DOF*numel(Numbering1{i1}));
        DP1=size(Numbering1{i1},1)-1;
        DP2=size(Numbering1{i1},2)-1;
        R=zeros(DP2+1,DP1+1);
   %%................................     
        for countx=1:length(wgx{i1})
            for county=1:length(wgy{i1})

    %%.........................................
                Xe_px=Re1{i1}{countx,county}(r,c);
%                 Xe_py=0;
                Xn_px=Rn1{i1}{countx,county}(r,c);
%                 Xn_py=0;
                Ye_px=0;
%                 Ye_py=Re1{i1}{countx,county}(r,c);
                Yn_px=0;
%                 Yn_py=Rn1{i1}{countx,county}(r,c);
    %%..............................................................................
yn1=Yn1{i1}(countx,county);
ye1=Ye1{i1}(countx,county);
xe1=Xe1{i1}(countx,county);
xn1=Xn1{i1}(countx,county);
                detJ1_px=Xe_px*yn1-Xn_px*ye1;
%                 detJ1_py=Yn_py*xe1-Ye_py*xn1;
   %%.................................................................................
JJ1=detJ11{i1}(countx,county);

    Rx_px=((Yn_px*Re1{i1}{countx,county}-Ye_px*Rn1{i1}{countx,county}*JJ1-...
        detJ1_px*(yn1*Re1{i1}{countx,county}-ye1*Rn1{i1}{countx,county}))/(JJ1^2));
    
%     Rx_py=((Yn_py*Re1{i1}{countx,county}-Ye_py*Rn1{i1}{countx,county})*JJ1-...
        detJ1_py*(yn1*Re1{i1}{countx,county}-ye1*Rn1{i1}{countx,county}))/(JJ1^2);
    
    Ry_px=((Xe_px*Rn1{i1}{countx,county}-Xn_px*Re1{i1}{countx,county})*JJ1-...
        detJ1_px*(-xn1*Re1{i1}{countx,county}+xe1*Rn1{i1}{countx,county}))/(JJ1^2);
    
%     Ry_py=((Xe_py*Rn1{i1}{countx,county}-Xn_py*Re1{i1}{countx,county})*JJ1-...
        detJ1_py*(-xn1*Re1{i1}{countx,county}+xe1*Rn1{i1}{countx,county}))/(JJ1^2);
   

   %%.................................................................................
                       
                B_px=buildB(DP1,DP2,DOF,R,Rx_px,Ry_px);
%                 B_py=buildB(DP1,DP2,DOF,R,Rx_py,Ry_py);
                 
    %%................................................................................. 

    DD=D1{i1}{countx,county};
    BB=B1{i1}{countx,county};
    
    a=(B_px'*DD*BB+BB'*DD*B_px)*JJ1;
    b=BB'*DD*BB*detJ1_px;
    
 k_px=k_px+(a+b)*wgx{i1}(countx)*wgy{i1}(county)*detJ2_1(i1);

%  a=(B_py'*DD*BB+BB'*DD*B_py)*JJ1;
%     b=BB'*DD*BB*detJ1_py;
 
%  k_py=k_py+(a+b)*wgx{i1}(countx)*wgy{i1}(county)*detJ2_1(i1);
  %%.................................................................................              
            end
        end
%%.................................................................................
        K_px{i1}=k_px;
%         K_py{i1}=k_py;
%%.................................................................................
    end
%%.................................................................................
    K_Px=Assembelingstifness(n_dof,Numbering1,K_px,DOF);
%     K_Py=Assembelingstifness(n_dof,Numbering1,K_py,DOF);
%%.................................................................................

end    