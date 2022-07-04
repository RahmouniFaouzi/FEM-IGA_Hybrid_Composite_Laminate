function [K] =  formStiffnessMatrixK(GDof,numberElements, ...
                elementNodes,numberNodes,nodeCoordinates, ...
                AMatrix,quadType,elemType)
            
% computation of stiffness matrix for Kirchhoff plate element
% K: stiffness matrix
K = zeros(GDof);

% Gauss quadrature 
[gaussWeights,gaussLocations] = gaussQuadrature(quadType);
% cycle for element
for e = 1:numberElements
    % indice: nodal connectivities for each element
    % elementDof: element degrees of freedom
    indice = elementNodes(e,:);
    elementDof = [indice indice+numberNodes];
    ndof = length(elementDof);
    
    % cycle for Gauss point
    for q = 1:size(gaussWeights,1)
        GaussPoint = gaussLocations(q,:);
        xi = GaussPoint(1);
        eta = GaussPoint(2);
        % part related to the mapping
        % shape functions and derivatives
        [~,natDerQ4] = shapeFunctionsQ(xi,eta,elemType);
        
        % Jacobian matrix, inverse of Jacobian,
        
        [Jacob,~,~ ] = JacobianK(nodeCoordinates(indice,:),natDerQ4);
        
        
        % Membrane
        B_m = call_Bm( xi,eta,elemType, nodeCoordinates(indice,:), ndof);
        
        % Stiffness Matrix
        %-----------------
        % ... membrane-membrane
        K(elementDof,elementDof) = K(elementDof,elementDof) + ...
            B_m'*AMatrix*B_m*gaussWeights(q)*det(Jacob);
        
        
    end % Gauss point
end % element
end

