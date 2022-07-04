% .............................................
% Finite Element Analysis Kirchhoff plate 5 DOF :
% ================================================

% clear memory :
% =============
clear; close all; warning('off');clc

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
elemType = 'Q4';
Nplies   = 8;

% theta = [ 0    0    0    0    0    0    0    0]; % TEST 01
% theta = [ 0   45    0   45   45    0   45    0]; % TEST 02
% theta = [45  -45   45  -45  -45   45  -45   45]; % TEST 03
% theta = [ 0   90    0   90   90    0   90    0]; % TEST 04
theta = [ 0   45  -45   90   90  -45   45    0]; % TEST 05
% theta = [ 0   30   60   90   90   60   30    0]; % TEST 06

CLS = [1 2 1 2 2 1 2 1];

h   = 1/8; % Ply Thickness
ht  = Nplies * h;

% load
% ------------------------------------------
%                  V
%           _______2_________
%          |                 |
%  1  N <- |                 | -> N   2
%          |_________________|
%                   1
%                   V
% 1 for applying load on side 1
% 2 for applying load on side 2
% 3 for applying load on side 1 and 2
% ------------------------------------

N = [1,  3]; % Along U loding
dof_per_node = 2; % number of kinematic parameters per node

% mesh generation
L = 100; W = 100;
Nx = 20;
Ny = Nx;
numberElements  = Nx*Ny;
[nodeCoordinates, elementNodes] = MeshRectanglularPlate(L,W,Nx,Ny);
numberNodes = size(nodeCoordinates,1);

% Plot Mesh
PlotMesh(nodeCoordinates,elementNodes,L,W)

% GDof: Global Number of Degrees of Freedom
GDof = dof_per_node * numberNodes;

for ii = 1 : size(theta,1)
    thetadt = theta(ii,:);
    fprintf('===================================================\n\n')
    fprintf(['Theta Sequense is: [', repmat('%g, ', 1, numel(thetadt)-1), '%g]\n'], thetadt)
    
    % ABD Matrix
    [AMatrix,qbarra] =  ABDMatrix(Nplies, thetadt, h, Prop, CLS);%
    
    % Stiffness Matrix
    [stiffness] = formStiffnessMatrixK(GDof,numberElements, ...
        elementNodes,numberNodes,nodeCoordinates, ...
        AMatrix,'complete',elemType);
    
    % Membrane Forces
    force  = zeros(GDof,1);
    force  = formForceVectorKUV(L,nodeCoordinates,force,N,Ny);
    displacements  = stiffness\force;
    fprintf('==> max_displacement  = %1.20f mm\n\n',max(displacements))
    
    % Maximum Stress Failure
    % ======================
    [stress_layer,rapport_max,Pli_etude,Mode,failure_C] = Stress_Max_Str( 8 , Nplies , numberElements,...
        elementNodes,numberNodes,nodeCoordinates,qbarra,...
        displacements,STRENGHT,thetadt,CLS,'complete');
    fprintf('==> Ply Failed = %1.0f  %s \n\n',Pli_etude , Mode)
    disp (failure_C)
    fprintf('==> Stress_Failure    = %4.3f \n\n',rapport_max)
    %========================================================
    
    % Maximum Strain Failure
    % ======================
    [rapport_max,Pli_etude,Mode,failure_C] = Strain_Max_Str( 8 , Nplies , numberElements,...
        elementNodes,numberNodes,nodeCoordinates,...
        displacements,STRENGHT,thetadt,CLS,'complete', Prop);
    
    disp (failure_C)
    fprintf('==> Stress_Failure    = %4.3f \n\n',rapport_max)
    %========================================================
    
    % Hoffman_Failure
    % ===============
    [~,Hoffman_max,~,effort_min,failure_C] = Stress_Hoffman( 8 , Nplies , numberElements,...
        elementNodes,numberNodes,nodeCoordinates,qbarra,...
        displacements,STRENGHT,thetadt,CLS,'complete');
    
    disp (failure_C)
    fprintf('==> Stress_Failure    = %4.3f \n\n',effort_min)
   %========================================================
    
    % Tsai_Hill_Failure
    % ===============
    [~,Tsai_Hill_max,~,effort_min,failure_C] = Stress_Tsai_Hill( 8 , Nplies , numberElements,...
        elementNodes,numberNodes,nodeCoordinates,qbarra,...
        displacements,STRENGHT,thetadt,CLS,'complete');
    disp (failure_C)
    fprintf('==> Stress_Failure    = %4.3f \n\n',effort_min)
    %=======================================================
    
    % Tsai_Wu_Failure
    % ===============
    [~,Tsai_Wu_max,~,effort_min,failure_C] = Stress_Tsai_Wu( 8 , Nplies , numberElements,...
        elementNodes,numberNodes,nodeCoordinates,qbarra,...
        displacements,STRENGHT,thetadt,CLS,'complete');
    disp (failure_C)
    fprintf('==> Stress_Failure    = %4.3f \n\n',effort_min)
    %=======================================================
    
    for i = 1: size(stress_layer,2)
        Pr1(i) = stress_layer(1,i,Pli_etude); %#ok<*SAGROW>
        Pr2(i) = stress_layer(2,i,Pli_etude);
        Pr3(i) = stress_layer(3,i,Pli_etude);
    end
    
    Ds = 1:size(stress_layer,2);
    figure(ii+1);
    hold on
    plot(Ds,Pr1,'-o');plot(Ds,Pr2,'-o');plot(Ds,Pr3,'-*');
    title('stress', 'FontSize', 20);
    legend('Sigma X', 'Sigma Y','Sigma XY');
end
