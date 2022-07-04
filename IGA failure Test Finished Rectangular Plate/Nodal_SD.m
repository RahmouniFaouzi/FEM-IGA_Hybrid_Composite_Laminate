function [Str ,Disp] = Nodal_SD(stress,disp ,node,elementV )

dispX = zeros(size(node,1),1);
dispY = zeros(size(node,1),1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid = connect(in);
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
    end
end

Disp = [dispX , dispY];


%%
InD = size (stress, 4);
for M = 1: InD
    sigmaXX = zeros(size(node,1),1);
    sigmaYY = zeros(size(node,1),1);
    sigmaXY = zeros(size(node,1),1);
    
    for e=1:size(elementV,1)
        connect = elementV(e,:);
        for in=1:4
            nid = connect(in);
            sigmaXX(nid) = stress(e,in,1, M);
            sigmaYY(nid) = stress(e,in,2, M);
            sigmaXY(nid) = stress(e,in,3, M);
        end
    end
    Str(:,:,M)  = [sigmaXX, sigmaYY,sigmaXY];
end

end