function [P,DIS,sig_M,epsilon_M,Final_D,Stiffness,fd,sd]=PlasticAnalysis(nli,PointDis,PointLoad,...
    K,KK,detJ1,detJ2,...
    guss_x,w_guss_x,guss_y,w_guss_y,DOF,...
    Numbering,t,D1,B,H,Ndof,noe,effective_stress)
residue1=1000*ones(Ndof,1);
fd= isnan(PointDis)==0;
sd= isnan(PointDis)==1;
DIS=zeros(Ndof,nli+1);
P=DIS;
D=cell(noe,1);
sig_M=cell(noe,nli+1);
sigeffecty1=cell(noe,1);
Stiffness=cell(50,1);
Final_D=cell(50,1);
for i=1:noe
    for j=1:size(B{i},1)
        for z=1:size(B{i},2)
            for ii=1:nli+1
    sig_M{i,ii}{j,z}=zeros(5,1);
            end
    D{i}{j,z}=D1{i};
    sigeffecty1{i}(j,z)=effective_stress(i);
        end
    end
end
dsiggg=sig_M;
epsilon_M=sig_M;
%% start correction
d0=PointDis/nli;
dF=PointLoad/nli;
for ii=1:nli
 residue=residue1;
  count=1;
while (all([any(abs(residue)>10),count<=10]))

fff=zeros(Ndof,1);
if count<=1
[~,d_d]=Solveequation(dF,d0,K);
end

for i=1:noe

number=Numbering{i}';
number=number(:);
df=zeros(DOF*numel(number),1);
dis=zeros(DOF*numel(number),1);
ind=1;
for j=1:numel(number)
   
        i1=(number(j)-1)*DOF+1:number(j)*DOF;
        dis(ind:ind+DOF-1,1)=d_d(i1);
        ind=ind+DOF;
end
for z1=1:size(B{i},1)
for z2=1:size(B{i},2)

depsilon=(B{i}{z1,z2})*dis;
epsilon_M{i,ii+1}{z1,z2}=epsilon_M{i,ii}{z1,z2}+depsilon;
dsig=D{i}{z1,z2}*depsilon;
sig_M{i,ii+1}{z1,z2}=sig_M{i,ii}{z1,z2}+dsig;
% M=sigma*t^2/4 & Q=  sig_M=[sigx;sigy;sigxy]
v_sig=sig_M{i,ii+1}{z1,z2};
sig_Mi1(1:3)=v_sig(1:3)*4/(t(i)^2);
sig_Mi1(4:5)=v_sig(4:5)*1.5/t(i);
sigeffect=vonmises(sig_Mi1);
v_sig=sig_M{i,ii}{z1,z2};
sigi(1:3)=v_sig(1:3)*4/(t(i)^2);
sigi(4:5)=v_sig(4:5)*1.5/t(i);
sigeffect_0=vonmises(sigi);
%%  check for yeilding at level r
error=sigeffect - sigeffecty1{i}(z1,z2);
if (error>0)
     epeffect=0;
     siggs=sig_M{i,ii}{z1,z2};
if sigeffect_0>=sigeffecty1{i}(z1,z2)
% calculate flow vector and D
             delta_e_ep=depsilon;
             v_sig=siggs;
             sig(1:3)=v_sig(1:3)*4/(t(i)^2);
             sig(4:5)=v_sig(4:5)*1.5/t(i);
             a=flowvector(sig);
             D2=D{i}{z1,z2};
             Dep=D2-(D2*(a*a')*D2)/(H(i)+a'*D2*a);
             D{i}{z1,z2}=Dep;
             dsig1=Dep*delta_e_ep;
             siggs=siggs+dsig1;
             v_sig=siggs;
             sig(1:3)=v_sig(1:3)*4/(t(i)^2);
             sig(4:5)=v_sig(4:5)*1.5/t(i);
             sigeffect=vonmises(sig);
             delta_e_p=delta_e_ep-D1{i}\dsig1;
             delta_e_p(1:3)=delta_e_p(1:3)*12/t(i)^3;
             delta_e_p(1:3)=delta_e_p(1:3)*1.5/t(i);
             epeffect=(2/3)*vonmises(delta_e_p);
             sig_M{i,ii+1}{z1,z2}=((sigeffecty1{i}(z1,z2))/sigeffect)*siggs;
             %% check for yeilding of level r-1 answer=no
else
    c=1;
while (error>0.01 && c<=5)

r=(sigeffect-sigeffecty1{i}(z1,z2))/(sigeffect-sigeffect_0);
h=1-r;
             siggs=sig_M{i,ii}{z1,z2}+h*dsig;
             delta_e_ep=r*depsilon;
             sig(1:3)=siggs(1:3)*4/(t(i)^2);
             sig(4:5)=siggs(4:5)*1.5/t(i);
% calculate flow vector and D 
             a=flowvector(sig);
             D2=D{i}{z1,z2};
             Dep=D2-(D2*(a*a')*D2)/(H(i)+a'*D2*a);
             D{i}{z1,z2}=Dep;
% calculate corrected sigma
             dsig1=Dep*delta_e_ep;
             delta_e_p=delta_e_ep-D1{i}\dsig1;
             delta_e_p(1:3)=delta_e_p(1:3)*12/t(i)^3;
             delta_e_p(1:3)=delta_e_p(1:3)*1.5/t(i);
             epeffect=(2/3)*vonmises(delta_e_p);
             sig_M{i,ii+1}{z1,z2}=siggs+dsig1;
             v_sig=sig_M{i,ii+1}{z1,z2};
             sig(1:3)=v_sig(1:3)*4/(t(i)^2);
             sig(4:5)=v_sig(4:5)*1.5/t(i);
             sigeffect=vonmises(sig);
             error=sigeffect - sigeffecty1{i}(z1,z2);
             dsig=sig_M{i,ii+1}{z1,z2}-sig_M{i,ii}{z1,z2};
              c=c+1;
             
end

end
sigeffecty1{i}(z1,z2)=sigeffecty1{i}(z1,z2)+H(i)*epeffect;
end

    dsiggg{i,ii}{z1,z2}=sig_M{i,ii+1}{z1,z2}-sig_M{i,ii}{z1,z2};
    
end
end
%% calculate sigma for each element and adding them to each other

for z2=1:size(B{i},2)
    
    wrn=w_guss_y{i}(z2);
    
    for z1=1:size(B{i},1)
        
        wre=w_guss_x{i}(z1);
        
        df=df+(B{i}{z1,z2})'*dsiggg{i,ii}{z1,z2}*wre*wrn*detJ1{i}(z1,z2)*detJ2(i);  
        
    end
    
end

ind=1;
for j=1:numel(number)
     i1=(number(j)-1)*DOF+1:number(j)*DOF;
     fff(i1)=fff(i1)+df(ind:ind+DOF-1,1);
     ind=ind+DOF;     
end
DP1=size(Numbering{i},2)-1;
DP2=size(Numbering{i},1)-1;
    BB=B{i};DD=D{i};J22=detJ2(i);J11=detJ1{i};
    [k]=CalculateStiffness2(BB,DD,DOF,DP1,DP2,...
        guss_x{i},w_guss_x{i},guss_y{i},w_guss_y{i},J22,J11);
    KK{i}=k;
end

%% Assembeling stifness

[K]=Assembelingstifness(Ndof,Numbering,KK,DOF);

if count<=1
    [ff,delta_d]=Solveequation(dF,d0,K);
    residue=fff-ff;
    count=count+1;
    d_d=delta_d;
end

if (all([max(abs(residue))>10,count>1]))

    residue=fff-ff;

 if max(abs(residue))<=1e5
    drfp=zeros(Ndof,1);
    drsp=zeros(Ndof,1);
    res_fp=drfp;
    res_sp=drfp;
    drfp(sd)=nan;
    drsp(fd)=nan;
    dr=drfp;
    res_fp(sd)=residue(sd);
    res_sp(fd)=residue(fd);
    [~,drfp]=Solveequation(-res_fp,drfp,K);
    [frsp,~]=Solveequation(-res_sp,drsp,K);
    frsp(fd)=0;
    [~,drsp1]=Solveequation(frsp,dr,K);
    d_d=d_d+drfp+drsp1;
    ff(fd)=K(fd,:)*d_d;
    count=count+1;
else
        ff=nan;
        count=11;
end
end
end
if ~any(isnan(ff))
P(:,ii+1)=P(:,ii)+ff;
DIS(:,ii+1)=DIS(:,ii)+d_d;
Stiffness{ii}=K;
Final_D{ii}=D;
else
P(:,ii+1:end)=nan;
DIS(:,ii+1:end)=nan;
break
end
end
end
