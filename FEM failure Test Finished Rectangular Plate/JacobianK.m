function [JacobianMatrix,invJacobian,XYDerivatives] = ...
    JacobianK(nodeCoordinates,naturalDerivatives)
% JacobianMatrix: Jacobian matrix
% ===============================

JacobianMatrix = nodeCoordinates'*naturalDerivatives(:,1:2);
invJacobian = inv(JacobianMatrix);
XYDerivatives = nodeCoordinates'*naturalDerivatives;

end % end function Jacobian

