function [dF,dG] = BuildSensitivity(U_Px,U_Py,A_px,A_py,wanted,COCP,symmetry_state)

    ux=zeros(size(COCP,1),1);
    uy=ux;ax=ux;ay=ux;
    %% Symmetric
    center_x_point=find(COCP(:,2)==0);
    [~,~,ia]=intersect(center_x_point,wanted);
    U_Px(ia)=0;
    A_px(ia)=0;
    center_y_point=find(COCP(:,3)==0);
    [~,~,ib]=intersect(center_y_point,wanted);
    U_Py(ib)=0;
    A_py(ib)=0;

%%
    ux(wanted)=U_Px;
    uy(wanted)=U_Py;
    ax(wanted)=A_px;
    ay(wanted)=A_py;


%%
if symmetry_state==1
    for i=1:length(wanted)
    index_x_1=find(COCP(:,2)==COCP(wanted(i),2));    
    index_x_2= COCP(index_x_1,3)==-COCP(wanted(i),3);
    index_x=setdiff(index_x_1(index_x_2),wanted(i));
    ux(index_x)=U_Px(i);
    uy(index_x)=-U_Py(i);
    ax(index_x)=A_px(i);
    ay(index_x)=-A_py(i);
    end
elseif symmetry_state==2
    for i=1:length(wanted)
    index_y_1=find(COCP(:,3)==COCP(wanted(i),3));    
    index_y_2= COCP(index_y_1,2)==-COCP(wanted(i),2);
    index_y=setdiff(index_y_1(index_y_2),wanted(i));
    ux(index_y)=-U_Px(i);
    uy(index_y)=U_Py(i);
    ax(index_y)=-A_px(i);
    ay(index_y)=A_py(i);
    end
else 
    
    for i=1:length(wanted)
    
    index_x_1=find(COCP(:,2)==COCP(wanted(i),2));    
    index_x_2= COCP(index_x_1,3)==-COCP(wanted(i),3);
    index_x=setdiff(index_x_1(index_x_2),wanted(i));
    index_y_1=find(COCP(:,3)==COCP(wanted(i),3));    
    index_y_2= COCP(index_y_1,2)==-COCP(wanted(i),2);
    index_y=setdiff(index_y_1(index_y_2),wanted(i));
    index_xy1=find(COCP(:,2)==-COCP(wanted(i),2));
    index_xy2=COCP(index_xy1,3)==-COCP(wanted(i),3);
    index_xy=setdiff(index_xy1(index_xy2),wanted(i));
    ux(index_x)=U_Px(i);
    uy(index_x)=-U_Py(i);
    ux(index_y)=-U_Px(i);
    uy(index_y)=U_Py(i);
    ux(index_xy)=-U_Px(i);
    uy(index_xy)=-U_Py(i);
    ax(index_x)=A_px(i);
    ay(index_x)=-A_py(i);
    ax(index_y)=-A_px(i);
    ay(index_y)=A_py(i);
    ax(index_xy)=-A_px(i);
    ay(index_xy)=-A_py(i);
    
    end
    
end

% 
%     ia=find(COCP(:,2)==100);
%     ib=find(COCP(:,2)==-100);
%     uy([ia;ib])=0;
%     ay([ia;ib])=0;
%     ia=find(COCP(:,3)==100);
%     ib=find(COCP(:,3)==-100);
%     ux([ia;ib])=0;
%     ax([ia;ib])=0;
    
%%
    dF=[ux;uy];
    dG=[ax;ay];
%% 
%     th=atan2(COCP(:,3),COCP(:,2));
%     dU=ux.*cos(th)+uy.*sin(th);
%     dA=ax.*cos(th)+ay.*sin(th);
%%

end

