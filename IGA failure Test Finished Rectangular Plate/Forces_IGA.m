function f = Forces_IGA(f ,F , noElemsV , elRangeV , elConnV , S_Point , ...
    Ed_Point , vKnot , weights, q)

[W1,Q1] = quadrature(q+1, 'GAUSS', 1 ); 

for e=1:noElemsV
    
   xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
   conn  = elConnV(e,:);   
   pts   = S_Point(conn,1:2);
   sctrx = Ed_Point(e,:);
   
   % loop over Gauss points 
    for gp=1:size(W1,1)                        
      xi      = Q1(gp,:);                          
      wt      = W1(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      [N, dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights');
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      
      jacob1   = dNdxi * pts;
      J1       = norm (jacob1);
      
      f(sctrx) = f(sctrx) + N' * F * J1 * J2 * wt;
      
    end
end