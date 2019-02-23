function [Ndof,DOF,DPx,DPy,N,Ux,Uy,ncpx,ncpy,ncp,norx,nory,...
    re,Wre,rn,Wrn,COCP,E,t,v,H,effective_stress,ne,noe,nli,noi,...
    pointload,pointdis,DVI,FPI,SS,SD,nx,ny,font_size]=SortData(inp)
%% Definition of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = Number of patches
% DOF = Degree of freedom for all of nodes
% DPx(i) = Degree of polynomial at x direction for i'th patch
% DPy(i) = Degree of polynomial at y direction for i'th patch
% ncpx(i) = Number of control points at x direction for i'th patch
% ncpy(i) = Number of control points at y direction for i'th patch
% norx(i) = Number of elements at x direction for i'th patch
% nory(i) = Number of elements at y direction for i'th patch
% ngx(i) = Number of Gaussian points at x direction for elements at i'th patch
% ngy(i) = Number of Gaussian points at y direction for elements at i'th patch
% ncp = Total number of control points
% nli = Number of load increment
% noi = Number of optimization iteration
% ne(i) = Number of elements for i'th patch
% noe = Total number of elements
% Ux{i} = Knot vector  at x direction for i'th patch
% Uy{i} = Knot vector  at y direction for i'th patch
% E(j) = Module of elasticity for j'th knot element
% t(j) = Thickness of element for j'th knot element
% V(j) = Poisson ratio of element for j'th knot element
% effective_stress(j) = Effective yield stress for j'th knot element
% H(j) = Strain hardening for j'th knot element
% re = Vector of gaussian points at x direction
% Wre = Gaussian weights at x direction
% rn = Vector of gaussian points at y direction
% Wrn = Gaussian weights at y direction
% COCP = Coordinates of control points
% ny = Number of points for each patch in order to plot shapes at y direction
% nx = Number of points for each patch in order to plot shapes at x direction
% font_size = Font size for plot results
% SS = Symmetry state
% SD = Sensitivity direction (1=X, 2=Y, 3=Both)
% DVI = Design variable index
% FPI = Fixed point index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code
N=inp(1,1);
DPx=zeros(N,1);
DPy=DPx;ncpx=DPx;ncpy=DPx;ngx=DPx;ngy=DPx;
re=cell(N,1);Wre=re;rn=re;Wrn=re;
counter=2;
for i=1:N
    DPx(i)=inp(counter,1);
    DPy(i)=inp(counter+1,1);
    ncpx(i)=inp(counter+2,1);
    ncpy(i)=inp(counter+3,1);
    ngx(i)=inp(counter+4,1);
    ngy(i)=inp(counter+5,1);
   [re{i},Wre{i}]=makegussianpoint(ngx(i));
   [rn{i},Wrn{i}]=makegussianpoint(ngy(i));
    counter=counter+6;
end
ncp=ncpx*ncpy';
norx=ncpx-DPx;
nory=ncpy-DPy;
ne=norx.*nory;
noe=sum(ne);
DOF=inp(counter,1);
nli=inp(counter+1,1);
noi=inp(counter+2,1);
ux=inp(:,2);uy=inp(:,3);
indx=find(isnan(ux),N,'first');indy=find(isnan(uy),N,'first');
Ux=cell(N,1);Uy=cell(N,1);
for i=1:N
    if i==1
        Ux{1}=ux(1:indx(1)-1);
        Uy{1}=uy(1:indy(1)-1);
    else
        Ux{i}=ux(indx(i-1)+1:indx(i)-1);
        Uy{i}=uy(indy(i-1)+1:indy(i)-1);
    end
end
clear ux uy
E=inp(1:noe+1,5);v=inp(1:noe+1,6);t=inp(1:noe+1,7);
effective_stress=inp(1:noe+1,8);H=inp(1:noe+1,9);
cocp=inp(:,10:12);
cocp(isnan(cocp(:,1)),:)=[];
COCP=cocp;
Ndof=max(COCP(:,1))*DOF;
%% Loads And Displacement
bodyforce=inp(:,13:16);
bodyforce(isnan(bodyforce(:,1)),:)=[];
tractionforce=inp(:,17:21);
tractionforce(isnan(tractionforce(:,1)),:)=[];
pointload=inp(:,22:25);
pointload(isnan(pointload(:,1)),:)=[];
pointdis=inp(:,26:29);
pointdis(isnan(pointdis(:,1)),:)=[];
%% Design Variable Index
DVI=inp(:,30);
DVI(isnan(DVI(:,1)),:)=[];
%% Fixed Point Index
FPI=inp(:,30);
FPI(isnan(FPI(:,1)),:)=[];
%% Symmetry State
SS=inp(1,32);
%% Sensitivty Direction
SD=inp(2,32);
%% Number Of Points And FontSize For Plot 
nx=inp(1,33);
ny=inp(2,33);
font_size=inp(3,33);
end
