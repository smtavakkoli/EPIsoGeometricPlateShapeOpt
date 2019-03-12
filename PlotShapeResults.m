function []=PlotShapeResults(DVI,FPI,npe,Numbering,norx,Ux,Uy,Px,Py,...
    COCP,Sig,t,nx,ny,fs_ty)
scrsz = get(groot,'ScreenSize');
figure('Position',[50 50 scrsz(3)/1.7 scrsz(4)/1.7]);
Color_plate=[0.5,0.5,0.5];

for i=1:2
    [X,Y] = CP2Phsical(npe,Numbering,norx,Ux,Uy,Px{i},Py{i},nx,ny);
    ax(i)=subplot(2,2,i);
    hold(ax(i),'on');
    xlabel({'X(mm)'},'FontSize',fs_ty);
    ylabel({'Y(mm)'},'FontSize',fs_ty);
    colormap(ax(i),Color_plate);
    ii=1;
    hold on;
    L=X;
    while ii<=sum(npe(:))

     contourf(X{ii},Y{ii},L{ii},'LineStyle','none');
     ii=ii+1;

    end
    LPI=setdiff(COCP{i}(:,1),[DVI;FPI]);
    if ~isempty(LPI)
    scatter(COCP{i}(LPI,2),COCP{i}(LPI,3),'MarkerFaceColor','black',...
        'MarkerEdgeColor','black')
    end
    if ~isempty(FPI)
    scatter(COCP{i}(FPI,2),COCP{i}(FPI,3),'MarkerFaceColor','black',...
        'MarkerEdgeColor','black')
    end
    scatter(COCP{i}(DVI,2),COCP{i}(DVI,3),'MarkerFaceColor','red',...
        'MarkerEdgeColor','red')
    for i_1=1:length(Numbering)

     for z1=1:size(Numbering{i_1},1)
         plot(COCP{i}(Numbering{i_1}(z1,:),2),COCP{i}(Numbering{i_1}(z1,:),3),'k')  
     end

     for z1=1:size(Numbering{i_1},2)
         plot(COCP{i}(Numbering{i_1}(:,z1),2),COCP{i}(Numbering{i_1}(:,z1),3),'k')  
     end
    end
    hold off
end
for i=3:4
    ax(i)=subplot(2,2,i);
    for z=1:sum(npe)
    z_1=size(Sig{i-2}{z,end},1);
    for ij=1:z_1
        z_2=size(Sig{i-2}{z,end},2);
        for j=1:z_2
                
             v_sig=Sig{i-2}{z,end}{ij,j};
             sig(1:3)=v_sig(1:3)*4/(t(ij)^2);
             sig(4:5)=v_sig(4:5)*1.5/t(ij);
             load1{z}(ij,j)=vonmises(sig);
        end
     end
    for ii=1:z_1
        Load1(ii,:)=strech_value(load1{z}(ii,:),ny);
    end
    for ii=1:size(Load1,2)
        Load{z,1}(:,ii)=strech_value(Load1(:,ii),nx);
    end
    end
nx=size(Load{z},1);
ny=size(Load{z},2);
[X,Y] = CP2Phsical(npe,Numbering,norx,Ux,Uy,Px{i-2},Py{i-2},nx-1,ny-1);
hold on;
 ii=1;
 while ii<=size(Load,1)
 contourf(X{ii},Y{ii},Load{ii},'LineStyle','none');
 ii=ii+1;
 end

xlabel({'X(mm)'},'FontSize',fs_ty);
ylabel({'Y(mm)'},'FontSize',fs_ty);
title({'Vonmises Stress(N/mm2)'},'FontSize',fs_ty);
colorbar('peer',ax(i));
 hold off
end
