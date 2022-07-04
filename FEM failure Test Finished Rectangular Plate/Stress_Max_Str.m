function [stress_layer,rapport_max,Pli_etude,Mode,failure_C] = Stress_Max_Str( ndof,Nlyer,...
    numberElements,elementNodes,numberNodes,nodeCoordinates,...
    qbarra,U,strength,theta,CLS, quadTypeB)


Str_G = zeros(3,Nlyer);
Strength_ratio =0;
Pli_etude =0;
Mode =0;
MK1 =zeros(3,Nlyer);
rapport_max =0;


count = 1;

% normal stresses in each layer
stress_layer = zeros(3,Nlyer,numberNodes);

% Gauss quadrature for bending part
[gaussWeights,gaussLocations] = gaussQuadrature(quadTypeB);

% cycle for element
for e = 1:numberElements
    
    indice = elementNodes(e,:);
    indiceB = [indice indice+numberNodes];
    
    % cycle for Gauss point
    for q = 1:size(gaussWeights,1)
        pt = gaussLocations(q,:);
        xi = pt(1);
        eta = pt(2);
        
        % Membrane
        B_m = call_Bm( xi,eta,'Q4', nodeCoordinates(indice,:), ndof);
        disp(B_m)
        % Strainess
        strain = B_m*U(indiceB);
        %disp(strain)
        
        % Stress
        for J = 1 : Nlyer
            stress_layer(:,J,count) = qbarra(1:3,1:3,J)*strain;
        end
        Str_G = stress_layer(:,:,count);
        count = count +1;
        
        % Failures Criteria
        [Strength_ratio,Pli_etude,Mode,~,rapport_max,failure_C] = Maximum_stress_failure...
        (Nlyer,CLS,Str_G,strength, theta, Strength_ratio ,Pli_etude,Mode,MK1,rapport_max);
    
    end % end Gauss point loop
end % end element loop
end
    
  
