function [dF,dG] =Symmetrization(COCP,SS,DVI,U_Px,U_Py,A_Px,A_Py)

COCP=round(COCP,5);
switch SS
    
    case 1
        for i=1:length(DVI)
            index_x_1=find(COCP(:,2)==COCP(DVI(i),2));    
            index_x_2= COCP(index_x_1,3)==-COCP(DVI(i),3);
            index_x=setdiff(index_x_1(index_x_2),DVI(i));
            U_Px(index_x)=U_Px(DVI(i));
            U_Py(index_x)=-U_Py(DVI(i));
            A_Px(index_x)=A_Px(DVI(i));
            A_Py(index_x)=-A_Py(DVI(i));
        end
    case 2
        for i=1:length(DVI)
            index_y_1=find(COCP(:,3)==COCP(DVI(i),3));    
            index_y_2= COCP(index_y_1,2)==-COCP(DVI(i),2);
            index_y=setdiff(index_y_1(index_y_2),DVI(i));
            U_Px(index_y)=-U_Px(i);
            U_Py(index_y)=U_Py(DVI(i));
            A_Px(index_y)=-A_Px(DVI(i));
            A_Py(index_y)=A_Py(DVI(i));
        end
    case 3
        for i=1:length(DVI)
            index_x_1=find(COCP(:,2)==COCP(DVI(i),2));    
            index_x_2= COCP(index_x_1,3)==-COCP(DVI(i),3);
            index_x=setdiff(index_x_1(index_x_2),DVI(i));
            index_y_1=find(COCP(:,3)==COCP(DVI(i),3));    
            index_y_2= COCP(index_y_1,2)==-COCP(DVI(i),2);
            index_y=setdiff(index_y_1(index_y_2),DVI(i));
            U_Px(index_x)=U_Px(DVI(i));
            U_Py(index_x)=-U_Py(DVI(i));
            U_Px(index_y)=-U_Px(DVI(i));
            U_Py(index_y)=U_Py(DVI(i));
            U_Px(index_xy)=-U_Px(DVI(i));
            U_Py(index_xy)=-U_Py(DVI(i));
            A_Px(index_x)=A_Px(DVI(i));
            A_Py(index_x)=-A_Py(DVI(i));
            A_Px(index_y)=-A_Px(DVI(i));
            A_Py(index_y)=A_Py(DVI(i));
       end
end
%% Sort Sensitivity Vectors
    dF=[U_Px;U_Py];
    dG=[A_Px;A_Py];
end
        