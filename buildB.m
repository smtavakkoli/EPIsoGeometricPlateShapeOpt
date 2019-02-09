function [B]=buildB(DP1,DP2,DOF,R,Rx,Ry)
%% Making matric B
B=zeros(5,(DP1+1)*(DP2+1)*DOF);
 for j=1:DP2+1;
 for i=1:DP1+1;
%Making matric B
     VB=[0,-Rx(j,i),0;0,0,-Ry(j,i);...
         0,-Ry(j,i),-Rx(j,i);Rx(j,i),-R(j,i),0;Ry(j,i),0,-R(j,i)];
      B(:,DOF*(i-1+(DP1+1)*(j-1))+1:DOF*(i+(DP1+1)*(j-1)))=VB;
 end
 end
end