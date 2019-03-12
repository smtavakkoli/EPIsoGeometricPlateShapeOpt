clear
clc
close all
%% Getting information from excel file
inp1=input('enter input file name: ','s');
inpp=strcat(inp1,'.xlsm');
inp=xlsread(inpp);
inp(end+1,:)=nan;
format short
%% Sorting and preparing data
[Ndof,DOF,DPx,DPy,N,Ux,Uy,ncpx,ncpy,ncp,norx,nory,...
    re,Wre,rn,Wrn,COCP{1},E,t,v,H,effective_stress,ne,noe,...
    nli,noi,pointload,pointdis,...
    DVI,FPI,SS,SD,nx,ny,font_size]=SortData(inp);
%% Define material properties
[Dmat,Numbering,Px{1},Py{1},KK,R,Re,Rn,Xe,Xn,Ye,Yn,Bmat,detJ1,detJ2,...
    guss_x,w_guss_x,guss_y,w_guss_y,COCP{1}]=DefineMaterialProperties...
    (N,COCP{1},DOF,E,t,v,ne,DPx,DPy,ncpx,ncpy,re,Wre,rn,Wrn,Ux,Uy,norx);
%% Assembeling stiffness of elements
[K]=Assembelingstifness(Ndof,Numbering,KK,DOF);
%% Assembling loads and displacements
[PointLoad,PointDis]=AssemblingLoad_Dis(Ndof,DOF,pointload,pointdis);
%% Plastic analysis
% fd = Address of known displacement
% sd = Address of unknown displacement
Force=cell(1);
Dis=cell(1);
Stress=cell(1);
Epsilon=cell(1);
[Force{1},Dis{1},Stress{1},Epsilon{1},Final_Dmat,Stiffness,fd,sd]...
    =PlasticAnalysis(nli,...
    PointDis,PointLoad,K,KK,detJ1,detJ2,...
    guss_x,w_guss_x,guss_y,w_guss_y,DOF,...
    Numbering,t,Dmat,Bmat,H,Ndof,noe,effective_stress);
%% Claculate energy dissipation
U=zeros(noi,1);
U(1)=CalculateEnergy(Force{1,1},Dis{1,1});
%% Calculate Area
Plate_Area=zeros(noi,1);
[Plate_Area(1)]=CalculateArea(detJ1,detJ2,guss_x,w_guss_x,guss_y,w_guss_y);
%% Optimization algoritm
if ~isnan(U(1))

% (X(i+1)=X(i)+d)
% X is vector of design varible by size m x 1
% m is total number of design varibles
% d = alpha*(dF-dG*lambda)
% alpha is step 
% dF is gradient vector of objective function by size m x 1
% dG is gradient matrix of constraint functions by size m x n
% n is total number of constraint functions
% lambda is vector of lagrang multiplier by size n x 1
% ' is transpose sign
    alpha=0.005;
    ii=2;
    G=0;
    while (ii<=noi && alpha<1)
    %% Sensitivity analysis
        [dF,dG] = SensitivityAnalysis(SD,SS,DVI,Dis{ii-1},Force{ii-1},...
            Stiffness,Final_Dmat,Bmat,Re,Rn,Xe,Xn,Ye,Yn,...
            detJ1,detJ2,Numbering,DOF,ncp,nli,...
            Ndof,fd,sd,w_guss_x,w_guss_y,COCP{ii-1});
    %% Calculation of lagrangMultiplier
        [lambda]=LambdaCalculation(dF,dG,G);
        if isnan(lambda)
            break
        else
            %% Update design variables
            while alpha<1
                d=alpha*(dF-dG*lambda);
                %% Different section
                % Because of variety in moving link points
                % this section of code is made just for example 1
                [cocp]=ChangeShape(N,COCP{ii-1},d,ncp,SS,ncpx,ncpy);
                %% Define material properties
                [Dmat,Numbering,px,py,KK,R,Re,Rn,Xe,Xn,Ye,Yn,Bmat,detJ1,detJ2,...
                 guss_x,w_guss_x,guss_y,w_guss_y,cocp]=DefineMaterialProperties(...
                    N,cocp,DOF,E,t,v,ne,...
                    DPx,DPy,ncpx,ncpy,re,Wre,rn,Wrn,Ux,Uy,norx);
                %% Check the constraints
                % Determine jacobi sign
                % If state=1 then jacobi>0
                % Else state=0 then jacobi<0
                state=CheckJacobian(detJ1);
                % calculate Area
                [plate_area]=CalculateArea(detJ1,detJ2,guss_x,w_guss_x,...
                    guss_y,w_guss_y);
                G=plate_area-Plate_Area(1);
                if (state==0 && G<=0)
                    clc
                    alpha=alpha+0.002
                    break
                end
                
                %% Assembeling stiffness of elements
                [K]=Assembelingstifness(Ndof,Numbering,KK,DOF);
                %% Plastic analysis
                [force,dis,sig,epsilon,final_Dmat,stiffness]=PlasticAnalysis(nli,...
                    PointDis,PointLoad,K,KK,detJ1,detJ2,...
                    guss_x,w_guss_x,guss_y,w_guss_y,DOF,...
                    Numbering,t,Dmat,Bmat,H,Ndof,noe,effective_stress);
                %% Claculate energy dissipation
                U_Test=CalculateEnergy(force,dis);
                if (~isnan(U_Test) && U_Test>U(ii-1))
                    U(ii)=U_Test;
                    Plate_Area(ii)=plate_area;
                    COCP{ii}=cocp;
                    Force{ii}=force;
                    Dis{ii}=dis;
                    Stress{ii}=sig;
                    Epsilon{ii}=epsilon;
                    Px{ii}=px;
                    Py{ii}=py;
                    Stiffness=stiffness;
                    Final_Dmat=final_Dmat;
                    %% Report situation
                    clc
                    ii
                    alpha
                    COCP{ii}(DVI,3)
                    %% Optimization iterations
                    ii=ii+1;
                    alpha=0.005;
                    break
                else
                    clc
                    alpha=alpha+0.002
                end
            end
        end

    end
    U(U==0)=U(ii-1);
    output_name=strcat(inp1,'_result');
    save(output_name,'U','Force','Dis','Sig','Epsilon',...
        'Numbering','Px','Py','R','COCP',...
        'Ux','Uy','norx','N','DVI');
    disp('Optimization is Completed');
    %% Plot shapes and contours
    PlotShapeResults(DVI,FPI,ne,Numbering,norx,Ux,Uy,...
        [Px(1),Px(end)],[Py(1),Py(end)],[COCP(1),COCP(end)],...
        [Stress(1),Stress(end)],t,nx,ny,font_size)
end
