function [FF,dd]=Solveequation(F,d,K)
fd=find(isnan(d)==0);
sd=find(isnan(d)==1);
dd=zeros(length(d),1);
FF=zeros(length(F),1);
dd(fd)=d(fd);
FF(sd)=F(sd);
dd(sd)=K(sd,sd)\(F(sd)-K(sd,fd)*d(fd));
FF(fd)=K(fd,:)*dd;
end
