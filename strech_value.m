function [c]=strech_value(c,n)
nc=length(c);
cc=[];
while n-nc>=2
    cc(1)=c(1);
    t=2;
  for i=1:nc-1
      
  cc(t)=(c(i)+c(i+1))/2;
  cc(t+1)=c(i+1);
  t=t+2;
  end

  c=cc;
  cc=[];
  nc=length(c);
  
end

end
        