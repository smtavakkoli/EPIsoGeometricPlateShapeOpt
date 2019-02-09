function [D]=DefineElasticityMatrix(E,t,v,n_e)
%% Making matric D

% E mudole of elasticity of element
% t=thickness of element,v=poisson ratio of element

D=cell(n_e,1);

      for i=1:n_e
      d=E(i)*t(i)^3/(12*(1-v(i)^2)); % d=Et^3/(12(1-v^2))
      G=E(i)/(2+2*v(i));
      s=G*t(i)/1.5;      % s=G*t/1.5
      D{i}(1,1)=d;D{i}(2,2)=d;D{i}(1,2)=v(i)*d;D{i}(2,1)=v(i)*d;
      D{i}(3,3)=(1-v(i))*d/2;D{i}(4,4)=s;D{i}(5,5)=s;
      end
      
end