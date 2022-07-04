function [S_Point,Ed_Point] =  Boundry_Par(controlPts,noElemsU,...
    noElemsV, Side , q, p ,L, D)

RightNodes = find(controlPts(:,1)==L)';
LeftNodes  = find(controlPts(:,1)==0)';
UpNodes    = find(controlPts(:,2)==D)';
DownNodes  = find(controlPts(:,2)==0)';

rightPoints     = controlPts(RightNodes,:);
leftPoints      = controlPts(LeftNodes,:);
upPoints        = controlPts(UpNodes,:);
downPoints      = controlPts(DownNodes,:);

rightEdgeMesh = zeros(noElemsV,q+1);
leftEdgeMesh  = zeros(noElemsV,q+1);
upEdgeMesh    = zeros(noElemsU,p+1);
downEdgeMesh  = zeros(noElemsU,p+1);
for i=1:noElemsV
    rightEdgeMesh(i,:)= RightNodes(i:i+q);
end
for i= 1 : noElemsV
    leftEdgeMesh(i,:) = LeftNodes(i:i+q);
end
for i= 1 : noElemsU
    upEdgeMesh(i,:)   = UpNodes(i:i+q);
end
for i= 1 : noElemsU
    downEdgeMesh(i,:) = DownNodes(i:i+q);
end

if strcmp(Side , 'Right')
    S_Point = rightPoints;
    Ed_Point= rightEdgeMesh;
    
elseif strcmp(Side , 'Left')
    S_Point = leftPoints;
    Ed_Point= leftEdgeMesh;
    
elseif strcmp(Side , 'Up')
    S_Point = upPoints;
    Ed_Point= upEdgeMesh;
    
elseif strcmp(Side , 'Down')
    S_Point = downPoints;
    Ed_Point= downEdgeMesh;
else
    disp ('error -->> INVALID NPUT')
end
end


