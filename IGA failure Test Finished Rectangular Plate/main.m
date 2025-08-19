close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       PRE-PROCESSING                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 100; % Lenght
D = 100; % Width
F = -1; % Force  1:Tension  -1:Compression
Side1      = 'Right'; % Side of applying Force
Side2      = 'Left';  % Side of applying Force

% Composite Laminate Material GR
E1 = 126e3; E2 = 11e3;
nu12 = 0.28; G12 = 6.6e3;
nu21 = nu12*E2/E1;
STRENGHT1 = [1950 1480 48 200 79];
Prop1 = [E1  E2  nu12  G12  nu21];

% Composite Laminate Material GL
E12 = 53.48e3; E22 = 17.7e3;
nu122 = 0.278; G122 = 5.83e3;
nu212 = nu122*E22/E12;
STRENGHT2 = [1140 570 35 114 72];
Prop2 = [E12 E22 nu122 G122 nu212];
%
Prop     = [Prop1 ; Prop2];
STRENGHT = [STRENGHT1 ; STRENGHT2];

Nplies  = 8;
 theta = [0 0 0 0 0 0 0 0]; %TEST 01
% theta = [0 45 0 45 45 0 45 0];% TEST 02
% theta = [0 90 0 90 90 0 90 0];% TEST 03
% theta = [45 -45 45 -45 -45 45 -45 45];% TEST 04
% theta = [0 45 -45 90 90 -45 45 0];% TEST 05
% theta = [0 30 60 90 90 60 30 0];% TEST 06

CLS = [1 2 1 2 2 1 2 1]; % 1 = GR; 2 = Gl.

h   = 1/8; % Ply Thickness
ht  = Nplies * h;

% control points
CtrlPts = zeros(4, 2, 2);
CtrlPts(1 : 3, 1, 1) = [0; 0; 0];
CtrlPts(1 : 3, 2, 1) = [L; 0; 0];
CtrlPts(1 : 3, 1, 2) = [0; D; 0];
CtrlPts(1 : 3, 2, 2) = [L; D; 0];
CtrlPts(4, :, :) = 1; % Wight

KntVect{1} = [0 0 1 1];
KntVect{2} = [0 0 1 1];


%%
% Refinements
% -------------
% degree of basis function
p = 2; q = 2;         % Order elevation --> P-rafinment

% repeated knot value inside knot vector
kx = 1; ky = 1;       % decrease continuity of knot

% number of elements in each direction
nelx = 10;
nely = 10;             % knot insertion  --> H-rafinment

% h+p = k-refinements
Surf = KRefine(KntVect, CtrlPts, [nelx, nely], [p, q], [p-kx, q-ky]);

figure
hold on
axis equal
daspect([1, 1, 1])
PlotGeo(Surf);
PlotKnts(Surf);
PlotCtrlPts(Surf);
PlotCtrlNet(Surf);

% Evaluation Nurbs Parametres
[uKnot, vKnot, noCtrPts,noElems,noElemsU,noElemsV,noDofs,...
    elRangeU, elConnU, elRangeV, elConnV ,element, index,...
    noGPs, controlPts,weights,p, q , noPtsX , noPtsY] = Evaluation_Nurbs_Par(Surf);

for ii = 1 : size(theta,1)
    thetadt = theta(ii,:);
    fprintf('===================================================\n\n')
    fprintf(['Theta Sequense is: [', repmat('%g, ', 1, numel(thetadt)-1), '%g]\n'], thetadt)
    
    [W,Q] = quadrature(  noGPs, 'GAUSS', 2 );
    
    % initialization
    [AMatrix,qbarra] =  ABDMatrix(Nplies, thetadt, h, Prop, CLS);%
    K = sparse(noDofs,noDofs);  % global stiffness matrix
    f = sparse(noDofs,1);        % external force vector
    
    jacob   = zeros(2,2);
    Nxi     = zeros(1,p+1);
    Neta    = zeros(1,q+1);
    dNdxi   = zeros(1,p+1);
    dNdeta  = zeros(1,q+1);
    
    % Loop over elements (knot spans)
    for e=1:noElems
        idu    = index(e,1);
        idv    = index(e,2);
        xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
        connU  = elConnU(idu,:);
        connV  = elConnV(idv,:);
        
        noFnsU = length(connU);
        noFnsV = length(connV);
        
        sctr   = element(e,:);          %  element scatter vector
        sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
        nn     = length(sctr);
        B      = zeros(3,2*nn);
        
        % loop over Gauss points
        for gp = 1:size(W,1)
            pt      = Q(gp,:);
            wt      = W(gp);
            
            % compute coords in parameter space
            Xi      = parent2ParametricSpace(xiE,pt(1));
            Eta     = parent2ParametricSpace(etaE,pt(2));
            J2      = jacobianPaPaMapping(xiE,etaE);
            
            % compute NURBSderivative of basis functions w.r.t parameter coord
            [dRdxi , dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');
            
            % compute the jacobian of physical and parameter domain mapping
            % then the derivative w.r.t spatial physical coordinates
            pts = controlPts(sctr,:);
            
            % Jacobian matrix
            jacob(1,1) = dRdxi  * pts(:,1);
            jacob(1,2) = dRdeta * pts(:,1);
            jacob(2,1) = dRdxi  * pts(:,2);
            jacob(2,2) = dRdeta * pts(:,2);
            
            J1         = det(jacob);
            
            % Jacobian inverse and spatial derivatives
            dRdx       = [dRdxi' dRdeta']/jacob;
            
            % B matrix
            %        _                                      _
            %        |  N_1,x  N_2,x  ...      0      0  ... |
            %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
            %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
            %        -                                      -
            
            B(1,1:nn)       = dRdx(:,1)';
            B(2,nn+1:2*nn)  = dRdx(:,2)';
            B(3,1:nn)       = dRdx(:,2)';
            B(3,nn+1:2*nn)  = dRdx(:,1)';
            
            % compute elementary stiffness matrix and
            % assemble it to the global matrix
            K(sctrB,sctrB) = K(sctrB,sctrB) + B' * AMatrix * B * J1 * J2 * wt; %#ok<*SPRIX>
        end
    end
    %
    % find boundary nodes for bounddary conditions
    [S_Point1,Ed_Point1] =  Boundry_Par(controlPts,noElemsU,...
        noElemsV, Side1,q, p ,L, D);
    [S_Point2,Ed_Point2] =  Boundry_Par(controlPts,noElemsU,...
        noElemsV, Side2,q, p ,L, D);
    
    % Campute Forces
    f = Forces_IGA( f ,F ,  noElemsV , elRangeV , elConnV , S_Point1 , ...
        Ed_Point1 , vKnot , weights, q);
    f = Forces_IGA( f ,-F ,  noElemsV , elRangeV , elConnV , S_Point2 , ...
        Ed_Point2 , vKnot , weights, q);
    
    % Displacment Of Contol Points
    U = K\f;
    
    % Campute Displacment and Stresses
    [stress, disp]  = Camp_STR_DISP(noElems , index, ...
        elRangeU , elRangeV ,element , noCtrPts, controlPts(:,1:2),noDofs, ...
        noPtsX , noPtsY , uKnot , vKnot ,weights , qbarra , U, p , q);
    
    % Visualization
    [node ,elementV] = buildVisualizationMesh ( controlPts(:,1:2), weights,...
        uKnot, vKnot, noPtsX, noPtsY , p, q );
    %
    
    Comp = 3; Ply = 2; 
    Visualisation_STR_DISP(stress, disp, node,Ply, elementV , Comp)
   
    % Recognize Stresses and Displacment of Real Space
    [Str ,Disp] = Nodal_SD(stress,disp ,node,elementV );
    
    %%
    fprintf('==> max_displacement  = %1.10f mm\n\n',max(Disp(:)))
    
    % Maximum Stress Failure
    % ======================
    [rapport_max1,Pli_etude,Mode,failure_C] = Stress_Max_StrISO(  Nplies,Str , STRENGHT, thetadt, CLS );
    
    fprintf('==> Ply Failed = %1.0f  %s \n\n',Pli_etude , Mode)
    
    fprintf (failure_C )
    fprintf('==> Stress_Failure    = %4.3f \n\n',rapport_max1)
    %======================================================================================================
    
    % Maximum Strain Failure
    % ======================
    [rapport_max2,Pli_etude,Mode,failure_C] = Strain_Max_StrISO(  Nplies,Str , STRENGHT, thetadt, CLS , Prop , qbarra);
    fprintf('==> Ply Failed = %1.0f  %s \n\n',Pli_etude , Mode)
    fprintf (failure_C )
    fprintf('==> Strain_Failure    = %4.3f \n\n',rapport_max2)
    %======================================================================================================
    
    % Hoffman_Failure
    % ===============
    [Hoffman_max,Pli_etude,effort_min1,failure_C] = Stress_HoffmanISO(  Nplies,Str , STRENGHT, thetadt, CLS );
    fprintf('==> Ply Failed = %1.0f  %s \n\n',Pli_etude , failure_C)
    fprintf (failure_C )
    fprintf('==> Stress_Failure    = %4.3f \n\n',effort_min1)
    %=======================================================================================================
    
    % Tsai_Hill_Failure
    % ===============
    [Tsai_Hill_max,Pli_etude,effort_min2,failure_C] = Stress_Tsai_HillISO( Nplies,Str , STRENGHT, thetadt, CLS );
    fprintf('==> Ply Failed = %1.0f  %s \n\n',Pli_etude , failure_C)
    fprintf (failure_C )
    fprintf('==> Stress_Failure    = %4.3f \n\n',effort_min2)
    %=======================================================================================================
    
    % Tsai_Wu_Failure
    % ===============
    [Tsai_Wu_max,Pli_etude,effort_min3,failure_C] = Stress_Tsai_WuISO(  Nplies,Str , STRENGHT, thetadt, CLS );
    fprintf('==> Ply Failed = %1.0f  %s \n\n',Pli_etude , failure_C)
    fprintf (failure_C )
    fprintf('==> Stress_Failure    = %4.3f \n\n',effort_min3)
    %=======================================================================================================
    Fald =[rapport_max1;rapport_max2;effort_min1;effort_min2;effort_min3];
    for i = 1: size(Str,3)
        Pr1(i) = Str(Pli_etude,1,i); %#ok<*SAGROW>
        Pr2(i) = Str(Pli_etude,2,i);
        Pr3(i) = Str(Pli_etude,3,i);
    end
    
    Ds = 1:size(Str,3);
    
    figure
    hold on
    plot(Ds,Pr1,'-o');plot(Ds,Pr2,'-o');plot(Ds,Pr3,'-*');
    title('stress', 'FontSize', 20);
    legend('Sigma X', 'Sigma Y','Sigma XY');
    hold off
end
H = gca;
H.LineWidth = 1.5; 

