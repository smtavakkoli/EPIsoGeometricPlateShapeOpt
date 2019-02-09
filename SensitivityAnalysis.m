function [dF,dG] =SensitivityAnalysis(SD,SS,DVI,DIS,P,Stiffness,Final_Dmat,...
                                    Bmat,Re,Rn,Xe,Xn,Ye,Yn,detJ1,detJ2,Numbering,...
                                    DOF,ncp,nli,Ndof,fd,sd,w_guss_x,w_guss_y,COCP)
%% Objective Function Sensitivity
% Preallocation Of Variables
    delta_px=cell(1);
    P_px=cell(1);
    U_px=zeros(nli,length(DVI));
    delta_py=cell(1);
    P_py=cell(1);
    U_py=zeros(nli,length(DVI));
    U_Px=zeros(ncp,1);
    U_Py=zeros(ncp,1);
% Choose Sensitivity Direction
% 1 = Calculation Sensitivity At X Direction
% 2 = Calculation Sensitivity At Y Direction
% 3 = Calculation Sensitivity At Boath Directions
switch SD
%................................................................
    case 1
       for counter=1:length(DVI)
            for ii=1:nli
                    if ii==1
                        P_px{ii,counter}=zeros(Ndof,1);
                    end
                    delta_d=DIS(:,ii+1)-DIS(:,ii);
                    k_inv=inv(Stiffness{ii}(sd,sd));
                    kk=Stiffness{ii};
                    [K_Px]= Stiffness_sensitivity_x(DVI(counter),...
                        Final_Dmat{ii},Bmat,Re,Rn,Xe,Xn,Ye,Yn,...
                        detJ1,detJ2,Numbering,DOF,Ndof,w_guss_x,w_guss_y);
                    delta_px{ii,counter}(fd,1)=0;
                    delta_px{ii,counter}(sd,1)=-k_inv*(K_Px(sd,:)*delta_d);
                    a=K_Px*delta_d+kk*delta_px{ii,counter};
                    P_px{ii+1,counter}=a+P_px{ii,counter};
                    a=P_px{ii+1,counter}+P_px{ii,counter};
                    b=(P(:,ii+1)+P(:,ii));
                    c=delta_px{ii,counter};
                    U_px(ii,counter)=0.5*(a'*delta_d+b'*c);
            end
        end
        U_Px(DVI)=(sum(U_px))';
%................................................................
    
    case 2
       for counter=1:length(DVI)
        for ii=1:nli
            if ii==1
                P_py{ii,counter}=zeros(Ndof,1);
            end
            delta_d=DIS(:,ii+1)-DIS(:,ii);
            k_inv=inv(Stiffness{ii}(sd,sd));
            kk=Stiffness{ii};
            [K_Py]= Stiffness_sensitivity_y(DVI(counter),...
                Final_Dmat{ii},Bmat,Re,Rn,Xe,Xn,Ye,Yn,...
                detJ1,detJ2,Numbering,DOF,Ndof,w_guss_x,w_guss_y);
            delta_py{ii,counter}(fd,1)=0;
            delta_py{ii,counter}(sd,1)=-k_inv*(K_Py(sd,:)*delta_d);
            a=K_Py*delta_d+kk*delta_py{ii,counter};
            P_py{ii+1,counter}=a+P_py{ii,counter};
            a=P_py{ii+1,counter}+P_py{ii,counter};
            b=(P(:,ii+1)+P(:,ii));
            c=delta_py{ii,counter};
            U_py(ii,counter)=0.5*(a'*delta_d+b'*c);
        end
      end
      U_Py(DVI)=(sum(U_py))';
        
%..........................................................................
    case 3
      for counter=1:length(DVI)
        for ii=1:nli
            if ii==1
                P_py{ii,counter}=zeros(Ndof,1);
                P_px{ii,counter}=zeros(Ndof,1);
            end
            delta_d=DIS(:,ii+1)-DIS(:,ii);
            k_inv=inv(Stiffness{ii}(sd,sd));
            kk=Stiffness{ii};
            [K_Py,K_Px]= Stiffness_sensitivity(DVI(counter),...
                Final_Dmat{ii},Bmat,Re,Rn,Xe,Xn,Ye,Yn,...
                detJ1,detJ2,Numbering,DOF,Ndof,w_guss_x,w_guss_y);
            delta_px{ii,counter}(fd,1)=0;
            delta_px{ii,counter}(sd,1)=-k_inv*(K_Px(sd,:)*delta_d);
            a=K_Px*delta_d+kk*delta_px{ii,counter};
            P_px{ii+1,counter}=a+P_px{ii,counter};
            a=P_px{ii+1,counter}+P_px{ii,counter};
            b=(P(:,ii+1)+P(:,ii));
            c=delta_px{ii,counter};
            U_px(ii,counter)=0.5*(a'*delta_d+b'*c);
            delta_py{ii,counter}(fd,1)=0;
            delta_py{ii,counter}(sd,1)=-k_inv*(K_Py(sd,:)*delta_d);
            a=K_Py*delta_d+kk*delta_py{ii,counter};
            P_py{ii+1,counter}=a+P_py{ii,counter};
            a=P_py{ii+1,counter}+P_py{ii,counter};
            b=(P(:,ii+1)+P(:,ii));
            c=delta_py{ii,counter};
            U_py(ii,counter)=0.5*(a'*delta_d+b'*c);
        end
      end
      U_Px(DVI)=(sum(U_px))';
      U_Py(DVI)=(sum(U_py))';
%..........................................................................
end
%% Constraint Function Sensitivity
A_Px=zeros(ncp,1);
A_Py=zeros(ncp,1);
switch SD
    case 1
        [A_px]= Area_sensitivity_x(DVI,Re,Rn,Xe,Xn,Ye,Yn,...
        detJ2,Numbering,w_guss_x,w_guss_y);
        A_Px(DVI)=A_px;
    case 2
        [A_py]= Area_sensitivity_y(DVI,Re,Rn,Xe,Xn,Ye,Yn,...
        detJ2,Numbering,w_guss_x,w_guss_y);
        A_Py(DVI)=A_py;
    case 3
        [A_px,A_py]= Area_sensitivity(DVI,Re,Rn,Xe,Xn,Ye,Yn,...
        detJ2,Numbering,w_guss_x,w_guss_y);
        A_Px(DVI)=A_px;
        A_Py(DVI)=A_py;
end
[dF,dG] = Symmetrization(COCP,SS,DVI,U_Px,U_Py,A_Px,A_Py);
end