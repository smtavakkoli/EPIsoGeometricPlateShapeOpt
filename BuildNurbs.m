function [NURBS]=BuildNurbs(P,U,a1,a2,j,u)
%%
% P is polynomial order
% U is open knot vector
% N is nurbs function with P=0
N=eye(max(size(U)));
N(1:P,:)=0;
N(end-P:end,:)=0;
if nargin==5
        N0=N;
 N00=vpa(N0);
 syms u;
elseif (nargin==6 && ischar(u))
        N0=N;
 N00=vpa(N0);
 u=vpa(u);
else
N0=N;
N00=zeros(P+2,1);
end
if (P~=0)
%noi is Number Of Intervals
%N0 and N00 are  copies of N for edit 
% noi=max(size(U))-2*P-1;
 for i=1:P;
   for k=a1:a2;
             
              a=((u-U(k))/(U(k+i)-U(k)));
              if (isinf(a)==1 || isnan(a)==1)
                      a=0;
              end
              b=((-u+U(k+i+1))/(U(k+i+1)-U(k+1)));
              if (isinf(b)==1 || isnan(b)==1)
                      b=0;
              end
             N00(k,j)=a*N0(k,j)+b*N0(k+1,j);
   end
   N0=N00;
   N0(k+1,j)=0;
  end
 NURBS=N00(a1:a2,j);
else
 NURBS=N(a1:a2,j);
 NURBS(isnan(NURBS)==1)=0;
end
