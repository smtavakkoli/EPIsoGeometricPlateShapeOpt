function [PointLoad,PointDis]=AssemblingLoad_Dis(Ndof,DOF,pointload,pointdis)
PointLoad=zeros(Ndof,1);
PointDis=nan(Ndof,1);
for i=1:size(pointload,1)
    pointload(i,isnan(pointload(i,2:4)))=0;
    index=(pointload(i,1)-1)*DOF+1:pointload(i,1)*DOF;
    PointLoad(index)=pointload(i,2:4)';
end
for i=1:size(pointdis,1)
    pointdis(i,isnan(pointdis(i,2:4)))=0;
    index=(pointdis(i,1)-1)*DOF+1:pointdis(i,1)*DOF;
    PointDis(index)=pointdis(i,2:4)';
end
PointLoad(~isnan(PointDis))=nan;
end