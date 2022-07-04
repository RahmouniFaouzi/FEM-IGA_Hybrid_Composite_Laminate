function [prescribedDof,activeDof] = BoundaryCond(Nodes, t, NBNd)


U  = unique(sort(Nodes ));
V  = unique(sort(Nodes+NBNd));

if     (t == 1)
    prescribedDof = U ;
    
elseif (t == 2)
    prescribedDof = [U V];
else 
    disp('erorr -->> Invalid Input')
end

activeDof = setdiff([1:NBNd*2]' , prescribedDof);

end