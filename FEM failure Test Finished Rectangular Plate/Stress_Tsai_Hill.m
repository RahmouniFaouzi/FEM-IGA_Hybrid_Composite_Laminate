function [stress_layer,Tsai_Hill_max,Pli_etude,effort_min,failure_C] = Stress_Tsai_Hill( ndof,Nlyer,...
    numberElements,elementNodes,numberNodes,nodeCoordinates,...
    qbarra,U,strength,theta,CLS, quadTypeB)

Tsai_Hill_max=0;
effort_min=0;
Pli_etude = 0;
MK3=zeros(1,Nlyer);
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
        
        % Strainess
        strain = B_m*U(indiceB);
        
        % Stress
        for J = 1 : Nlyer
            stress_layer(:,J,count) = qbarra(1:3,1:3,J)*strain;
        end
        Str_G = stress_layer(:,:,count);
        count = count +1;
        
        % Failures Criteria
        [Tsai_Hill_max,Pli_etude,MK3,effort_min, failure_C]=Tsai_Hill_failure(Nlyer,Str_G,CLS,strength, theta,Tsai_Hill_max,Pli_etude,MK3,effort_min);
      
    end % end Gauss point loop
end % end element loop
end
    

