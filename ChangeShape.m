function [cocp]=ChangeShape(N,COCP,d,ncp,SS,ncpx,ncpy)
cocp=COCP;
%% Movement Of Control Points Points
if N==1
    switch SS
        
        case 1
                
            cocp(:,3)=cocp(:,3)+d(ncp+1:end);
            x_vector=unique(cocp(:,2));
            
            for i=1:length(x_vector)
              index=find(cocp(:,2)==x_vector(i));
              index=sort(index);
              max_y=max(cocp(index,3));
              min_y=min(cocp(index,3));
              delta1=(max_y-min_y)/(ncpy-1);
              delta2=delta1*(0:ncpy-1)';
              cocp(index,3)=min_y+delta2;
            end
            
        case 2
            cocp(:,2)=cocp(:,2)+d(1:ncp);
            y_vector=unique(cocp(:,3));
            
            for i=1:length(y_vector)
              index=find(cocp(:,3)==y_vector(i));
              index=sort(index);
              max_x=max(cocp(index,2));
              min_x=min(cocp(index,2));
              delta1=(max_x-min_x)/(ncpx-1);
              delta2=delta1*(0:ncpx-1)';
              cocp(index,3)=min_x+delta2;
            end
            
        case 3
            cocp(:,3)=cocp(:,3)+d(ncp+1:end);
            x_vector=unique(cocp(:,2));
            
            for i=1:length(x_vector)
              index=find(cocp(:,2)==x_vector(i));
              index=sort(index);
              max_y=max(cocp(index,3));
              min_y=min(cocp(index,3));
              delta1=(max_y-min_y)/(ncpy-1);
              delta2=delta1*(0:ncpy-1)';
              cocp(index,3)=min_y+delta2;
            end
            
            cocp(:,2)=cocp(:,2)+d(1:ncp);
            y_vector=unique(cocp(:,3));
            
            for i=1:length(y_vector)
              index=find(cocp(:,3)==y_vector(i));
              index=sort(index);
              max_x=max(cocp(index,2));
              min_x=min(cocp(index,2));
              delta1=(max_x-min_x)/(ncpx-1);
              delta2=delta1*(0:ncpx-1)';
              cocp(index,3)=min_x+delta2;
            end
        case 4
            cocp(:,2)=cocp(:,2)+d(1:ncp);
            cocp(:,3)=cocp(:,3)+d(ncp+1:end);
    end
    
end
end
