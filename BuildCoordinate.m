function [P1,P2]=BuildCoordinate(Numbering,COCP)
    P1=cell(length(Numbering),1);
    P2=P1;
    for i=1:length(Numbering)
    
        for j=1:size(Numbering{i},1)
            for k=1:size(Numbering{i},2)
           ind=find(COCP(:,1)==Numbering{i}(j,k),1,'first');
           P1{i}(j,k)=COCP(ind,2);
           P2{i}(j,k)=COCP(ind,3);
           
            end
        end
    end