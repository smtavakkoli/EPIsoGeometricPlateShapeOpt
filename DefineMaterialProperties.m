function [Dmat,Numbering,Px,Py,KK,R,Re,Rn,Xe,Xn,Ye,Yn,Bmat,detJ1,detJ2,...
    guss_x,w_guss_x,guss_y,w_guss_y,COCP]=DefineMaterialProperties(N,COCP,DOF,...
    E,t,v,ne,DPx,DPy,ncpx,ncpy,re,Wre,rn,Wrn,Ux,Uy,norx)
%% Definition of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmat{i} = Elasticity matrix for i'th element
% Numbering{i} = Numbers of control points for i'th element
% Px{i} = X_coordinates of control points for i'th element
% Py{i} = Y_coordinates of control points for i'th element
% R{i} = Shape function values for i'th element
% Re{i} = Derivative of shape function with respect to e for i'th element
% Rn{i} = Derivative of shape function with respect to e for i'th element
% Xe{i} = Partial derivative x with respect to e for i'th element
% Xn{i} = Partial derivative x with respect to n for i'th element
% Ye{i} = Partial derivative y with respect to e for i'th element
% Yn{i} = Partial derivative y with respect to n for i'th element
% Bmat{i} = Strain-Displacement matrix
% KK{i} = Stiffness matrix for i'th element
% guss_x{i} = Gaussian points at X direction for i'th element
% w_guss_x{i} = Gaussian weights at x direction for i'th element
% guss_y{i} = Gaussian points at y direction for i'th element
% w_guss_y{i} = Gaussian weights at y direction for i'th element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
guss_x=[];
w_guss_x=[];
guss_y=[];
w_guss_y=[];
for i=1:N
    if i==1
        ind_cp=COCP(1:ncpx(i)*ncpy(i),1);
        range=1:ne(i);
        is=0;
        ie=0;
    else
        is=is+ncpx(i-1)*ncpy(i-1);
        ie=is+ncpx(i)*ncpy(i);
        ind_cp=COCP(is+1:ie);
        range=sum(ne(1:i-1))+1:sum(ne(1:i-1))+ne(i);
    end
[Dmat(range,1)]=DefineElasticityMatrix(E(range),t(range),v(range),ne(i));
[Numbering(range,1)]=numbering(DPx(i),DPy(i),ncpx(i),norx(i),ne(i),ind_cp);
[Px(range,1),Py(range,1)]=BuildCoordinate(Numbering(range),COCP);
[KK(range,1),R(range,1),Re(range,1),Rn(range,1),...
 Xe(range,1),Xn(range,1),Ye(range,1),Yn(range,1),...
 Bmat(range,1),detJ1(range,1),detJ2(range,1)]...
=CalculateStiffness(Dmat(range,1),DOF,DPx(i),DPy(i),...
re{i},Wre{i},rn{i},Wrn{i},Ux{i},Uy{i},ne(i),Px(range),Py(range),norx(i));
guss_x=[guss_x;repmat(re(i),ne(i),1)];
w_guss_x=[w_guss_x;repmat(Wre(i),ne(i),1)];
guss_y=[guss_y;repmat(rn(i),ne(i),1)];
w_guss_y=[w_guss_y;repmat(Wrn(i),ne(i),1)];
end
COCP=round(COCP,5);
index=COCP(:,1);
index=unique(index);
COCP=COCP(index,:);
clear re Wre rn Wrn
end