function [Strength_ratio,Pli_etude,Mode,MK1,rapport_max,failure_C]=Maximum_stress_failure(Npli,MATP,Sig,K , theta,...
                                                            Strength_ratio ,Pli_etude,Mode,MK1,rapport_max)

% For each direction of solicitation MK1 gathers the strength ratio for each ply of the laminate 
failure_C  = 'Maximum_stress_failure';
Sig_mat = zeros(3,Npli);

for Colonne=1:Npli
    
    angle = theta(Colonne);
    
    % the angle of the lamina
    s = sind(angle); c = cosd(angle);
    
    % Inverse Transformation Matrix
    T=[c^2, s^2 , 2*s*c;
       s^2, c^2 ,-2*s*c;
      -s*c, s*c ,(c^2)-(s^2)];
    
    % Calculating the local strains
    Sig_mat(:,Colonne) = T*Sig(:,Colonne); 
    
    % Traction/compression along X
    if Sig_mat(1,Colonne)> 0
        MK1(1,Colonne)= abs(Sig_mat(1,Colonne)/K(MATP(Colonne),1));
        rapport=abs(K(MATP(Colonne),1)/Sig_mat(1,Colonne));
        mode_int='traction along X';
    else
        MK1(1,Colonne)= abs(Sig_mat(1,Colonne)/K(MATP(Colonne),2));
        mode_int='compression along X';
        rapport=abs(K(MATP(Colonne),2)/Sig_mat(1,Colonne));
    end
    
    if MK1(1,Colonne)>Strength_ratio
        Strength_ratio=MK1(1,Colonne);
        Pli_etude=Colonne;
        Mode=mode_int;
        rapport_max=rapport;
    end
    
    % Traction/compression along Y
    
    if Sig_mat(2,Colonne)> 0
        MK1(2,Colonne)= abs(Sig_mat(2,Colonne)/K(MATP(Colonne),3));
        mode_int='traction along Y';
        rapport=abs(K(MATP(Colonne),3)/Sig_mat(2,Colonne));
    else
        MK1(2,Colonne)= abs(Sig_mat(2,Colonne)/K(MATP(Colonne),4));
        mode_int='compression along Y';
        rapport=abs(K(MATP(Colonne),4)/Sig_mat(2,Colonne));
    end
    
    if MK1(2,Colonne)>Strength_ratio
        Strength_ratio=MK1(2,Colonne);
        Pli_etude=Colonne;
        Mode=mode_int;
        rapport_max=rapport;
    end
    
    % Shear
    
    if Sig_mat(3,Colonne)> 0
        MK1(3,Colonne)= abs(Sig_mat(3,Colonne)/K(MATP(Colonne),5));
        rapport=abs(K(MATP(Colonne),5)/Sig_mat(3,Colonne));
    else
        MK1(3,Colonne)= abs(Sig_mat(3,Colonne)/K(MATP(Colonne),5));
        rapport=abs(K(MATP(Colonne),5)/Sig_mat(3,Colonne));
    end
    
    if MK1(3,Colonne)>Strength_ratio
        Strength_ratio=MK1(3,Colonne);
        Pli_etude=Colonne;
        Mode='shear';
        rapport_max=rapport;
    end
end

end