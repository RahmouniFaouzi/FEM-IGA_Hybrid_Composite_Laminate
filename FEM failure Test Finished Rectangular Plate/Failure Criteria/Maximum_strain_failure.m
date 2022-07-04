function [Strength_ratio,Pli_etude,Mode,MK1,rapport_max,failure_C] = Maximum_strain_failure(Npli,MATP,Strain_L,K , theta,...
                                                            Strength_ratio ,Pli_etude,Mode,MK1,rapport_max)

% For each direction of solicitation MK1 gathers the strength ratio for each ply of the laminate 
failure_C  = 'Maximum_strain_failure';
Strain = zeros(3,Npli);

for Colonne=1:Npli
    
    angle = theta(Colonne);
    
    % the angle of the lamina
    s = sind(angle); c = cosd(angle);
    
    % Inverse Transformation Matrix
    T=[c^2, s^2 , 2*s*c;
       s^2, c^2 ,-2*s*c;
      -s*c, s*c ,(c^2)-(s^2)];
    
    % Reuter Matrix
    R=[1 0 0; 0 1 0; 0 0 2];
    
    % Calculating the local strains
    Strain = R*T*inv(R) * Strain_L;
    
    % Traction/compression along X
    if Strain(1)> 0
        MK1(1)= abs(Strain(1)/K(MATP(Colonne),1));
        rapport=abs(K(MATP(Colonne),1)/Strain(1));
        mode_int='traction along X';
    else
        MK1(1)= abs(Strain(1)/K(MATP(Colonne),2));
        mode_int='compression along X';
        rapport=abs(K(MATP(Colonne),2)/Strain(1));
    end
    
    if MK1(1)>Strength_ratio
        Strength_ratio=MK1(1);
        Pli_etude=Colonne;
        Mode=mode_int;
        rapport_max=rapport;
    end
    
    % Traction/compression along Y
    
    if Strain(2)> 0
        MK1(2)= abs(Strain(2)/K(MATP(Colonne),3));
        mode_int='traction along Y';
        rapport=abs(K(MATP(Colonne),3)/Strain(2));
    else
        MK1(2)= abs(Strain(2)/K(MATP(Colonne),4));
        mode_int='compression along Y';
        rapport=abs(K(MATP(Colonne),4)/Strain(2));
    end
    
    if MK1(2)>Strength_ratio
        Strength_ratio=MK1(2);
        Pli_etude=Colonne;
        Mode=mode_int;
        rapport_max=rapport;
    end
    
    % Shear
    
    if Strain(3)> 0
        MK1(3)= abs(Strain(3)/K(MATP(Colonne),5));
        rapport=abs(K(MATP(Colonne),5)/Strain(3));
    else
        MK1(3)= abs(Strain(3)/K(MATP(Colonne),5));
        rapport=abs(K(MATP(Colonne),5)/Strain(3));
    end
    
    if MK1(3)>Strength_ratio
        Strength_ratio=MK1(3);
        Pli_etude=Colonne;
        Mode='shear';
        rapport_max=rapport;
    end
end

end