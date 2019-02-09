function [R,Re,Rn,Ry,Rx,Xe,Xn,Ye,Yn,detJ1]=...
    DefineSortNurbs(DPx,Ux,inde,e,DPy,Uy,indn,n,Px,Py)
%% Calculation of nurbs function
% N(i) = Nurbs function value for i'th control point at x direction
% M(j) = Nurbs function value for j'th control point at y direction
% R = M*N'
 
% R(j,i) = shape function value for i'th control point at x direction and
% j'th control point at y direction
 
% X = R*Px'
% Y = R*Py'
 
% dN(i) = Derivative of nurbs function with respect to e ...
% for i'th control point at x direction
 
% dM(j) = Derivative of nurbs function with respect to n ...
% for j'th control point at y direction
 
% Xe = Partial derivative x with respect to e 
% Ye = Partial derivative y with respect to e
% Xn = Partial derivative x with respect to n 
% Yn = Partial derivative y with respect to n
% detJ1 = Determinant of jacobian matrix
 
% Re(j,i) = Derivative of shape function with respect to e for i'th control
% point at x direction and j'th control point at y direction
 
% Rn(j,i) = Derivative of shape function with respect to n for i'th control
% point at x direction and j'th control point at y direction
 
% Rx(j,i) = Derivative of shape function with respect to x for i'th control
% point at x direction and j'th control point at y direction
 
% Ry(j,i) = Derivative of shape function with respect to y for i'th control
% point at x direction and j'th control point at y direction


N=BuildNurbs(DPx,Ux,inde-DPx,inde,inde,e);
M=BuildNurbs(DPy,Uy,indn-DPy,indn,indn,n);
%%
if nargout==1
%Calculate Jacobian
%detJ = determinant of jacobi
R=M*N';
else 
%% Calculate Partial Derivative
R=M*N';
dN=DerivativeNurbs(DPx,Ux,inde,inde,1,e);
dM=DerivativeNurbs(DPy,Uy,indn,indn,1,n);
Re=M*dN';
Rn=dM*N';
Xe=M'*Px*dN;
Xn=dM'*Px*N;
Ye=M'*Py*dN;
Yn=dM'*Py*N;
detJ1=Xe*Yn-Ye*Xn;
Rx=(1/detJ1)*(Yn*Re-Ye*Rn);
Ry=(1/detJ1)*(-Xn*Re+Xe*Rn);
end




 