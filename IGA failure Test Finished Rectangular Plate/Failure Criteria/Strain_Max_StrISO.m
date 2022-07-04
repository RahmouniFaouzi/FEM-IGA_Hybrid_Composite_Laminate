function [rapport_max,Pli_etude,Mode,failure_C] = Strain_Max_StrISO( Nlyer,STR, ...
    strength,theta,CLS ,Prop , C )

Strength_ratio =0;
Pli_etude =0;
Mode =0;
MK1 =zeros(3,Nlyer);
rapport_max =0;

% Strain Strength
Div = [Prop(1,1),Prop(1,1),Prop(1,2),Prop(1,2),Prop(1,4);
       Prop(2,1),Prop(2,1),Prop(2,2),Prop(2,2),Prop(2,4)];
strength = strength./Div;

% Strain is the same on all layers
Str_G = STR(1,:,1)';
STRAIN = C(:,:,1)\Str_G;


% Failures Criteria
[Strength_ratio,Pli_etude,Mode,~,rapport_max,failure_C] = Maximum_strain_failureISO...
    (Nlyer,CLS,STRAIN,strength, theta, Strength_ratio ,Pli_etude,Mode,MK1,rapport_max);


end

% function [rapport_max,Pli_etude,Mode,failure_C] = Strain_Max_Str( ndof,Nlyer,...
%     numberElements,elementNodes,numberNodes,nodeCoordinates,...
%     U,strength,theta,CLS, quadTypeB , Prop)
% 
% 
% 
% Strength_ratio =0;
% Pli_etude =0;
% Mode =0;
% MK1 =zeros(3,Nlyer);
% rapport_max =0;
% 
% Div = [Prop(1,1),Prop(1,1),Prop(1,2),Prop(1,2),Prop(1,4);
%        Prop(2,1),Prop(2,1),Prop(2,2),Prop(2,2),Prop(2,4)];
% strength = strength./Div;
% 
% 
% % Gauss quadrature for bending part
% [gaussWeights,gaussLocations] = gaussQuadrature(quadTypeB);
% 
% % cycle for element
% for e = 1:numberElements
%     
%     indice = elementNodes(e,:);
%     indiceB = [indice indice+numberNodes];
%     
%     % cycle for Gauss point
%     for q = 1:size(gaussWeights,1)
%         pt = gaussLocations(q,:);
%         xi = pt(1);
%         eta = pt(2);
%         
%         % Membrane
%         B_m = call_Bm( xi,eta,'Q4', nodeCoordinates(indice,:), ndof);
%         
%         % Strainess
%         strain = B_m*U(indiceB);
%         
%         % Failures Criteria
%         [Strength_ratio,Pli_etude,Mode,~,rapport_max,failure_C] = Maximum_strain_failure...
%         (Nlyer,CLS,strain,strength, theta, Strength_ratio ,Pli_etude,Mode,MK1,rapport_max);
%     
%     end % end Gauss point loop
% end % end element loop
% end