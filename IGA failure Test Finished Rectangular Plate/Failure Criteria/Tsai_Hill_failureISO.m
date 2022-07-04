function [Tsai_Hill_max,Pli_etude,MK3,effort_min,failure_C ]=Tsai_Hill_failure(Npli,Sig,MATP,K,theta,Tsai_Hill_max,Pli_etude,MK3,effort_min)

failure_C  = 'Tsai_Hill_failure \n';
Sig_mat = zeros(3,Npli);

for Colonne=1:Npli
    
    angle = theta(Colonne);
    
    % the angle of the lamina
    s=sind(angle); c=cosd(angle);
    
    % Inverse Transformation Matrix
    T=[c^2, s^2 , 2*s*c;
       s^2, c^2 ,-2*s*c;
      -s*c, s*c ,(c^2)-(s^2)];
    
    % Calculating the local strains
    Sig_mat(:,Colonne) = T*Sig(:,Colonne);
    
                             % S1^2/X^2 - S1*S2/X^2 + S2^2/Y^2 + S12^2/S^2
    if Sig_mat(1,Colonne)>0     
        if Sig_mat(2,Colonne)>0  % S1^2/Xt^2 - S1*S2/Xt^2 + S2^2/Yt^2 + S12^2/S^2 % 1 1 3 5
            MK3(1,Colonne)=(Sig_mat(1,Colonne)^2)/(K(MATP(Colonne),1)^2)-Sig_mat(1,Colonne)*Sig_mat(2,Colonne)/K(MATP(Colonne),1)^2+Sig_mat(2,Colonne)^2/K(MATP(Colonne),3)^2+Sig_mat(3,Colonne)^2/K(MATP(Colonne),5)^2;
            effort=1/sqrt((Sig_mat(1,Colonne)^2)/(K(MATP(Colonne),1)^2)-Sig_mat(1,Colonne)*Sig_mat(2,Colonne)/K(MATP(Colonne),1)^2+Sig_mat(2,Colonne)^2/K(MATP(Colonne),3)^2+Sig_mat(3,Colonne)^2/K(MATP(Colonne),5)^2);
        else                 % S1^2/Xt^2 - S1*S2/Xt^2 + S2^2/Yc^2 + S12^2/S^2 % 1 1 4 5
            MK3(1,Colonne)=(Sig_mat(1,Colonne)^2)/(K(MATP(Colonne),1)^2)-Sig_mat(1,Colonne)*Sig_mat(2,Colonne)/K(MATP(Colonne),2)^2+Sig_mat(2,Colonne)^2/K(MATP(Colonne),4)^2+Sig_mat(3,Colonne)^2/K(MATP(Colonne),5)^2;
            effort=1/sqrt((Sig_mat(1,Colonne)^2)/(K(MATP(Colonne),1)^2)-Sig_mat(1,Colonne)*Sig_mat(2,Colonne)/K(MATP(Colonne),1)^2+Sig_mat(2,Colonne)^2/K(MATP(Colonne),4)^2+Sig_mat(3,Colonne)^2/K(MATP(Colonne),5)^2);
        end
    else
        if Sig_mat(2,Colonne)>0   % S1^2/Xc^2 - S1*S2/Xc^2 + S2^2/Yt^2 + S12^2/S^2 % 2 2 3 5
            MK3(1,Colonne)=(Sig_mat(1,Colonne)^2)/(K(MATP(Colonne),2)^2)-Sig_mat(1,Colonne)*Sig_mat(2,Colonne)/K(MATP(Colonne),2)^2+Sig_mat(2,Colonne)^2/K(MATP(Colonne),3)^2+Sig_mat(3,Colonne)^2/K(MATP(Colonne),5)^2;
            effort=1/sqrt((Sig_mat(1,Colonne)^2)/(K(MATP(Colonne),2)^2)-Sig_mat(1,Colonne)*Sig_mat(2,Colonne)/K(MATP(Colonne),2)^2+Sig_mat(2,Colonne)^2/K(MATP(Colonne),3)^2+Sig_mat(3,Colonne)^2/K(MATP(Colonne),5)^2);
        else                  % S1^2/Xc^2 - S1*S2/Xc^2 + S2^2/Yc^2 + S12^2/S^2 % 2 2 4 5
            MK3(1,Colonne)=(Sig_mat(1,Colonne)^2)/(K(MATP(Colonne),2)^2)-Sig_mat(1,Colonne)*Sig_mat(2,Colonne)/K(MATP(Colonne),1)^2+Sig_mat(2,Colonne)^2/K(MATP(Colonne),4)^2+Sig_mat(3,Colonne)^2/K(MATP(Colonne),5)^2;
            effort=1/sqrt((Sig_mat(1,Colonne)^2)/(K(MATP(Colonne),2)^2)-Sig_mat(1,Colonne)*Sig_mat(2,Colonne)/K(MATP(Colonne),2)^2+Sig_mat(2,Colonne)^2/K(MATP(Colonne),4)^2+Sig_mat(3,Colonne)^2/K(MATP(Colonne),5)^2);
        end
    end
    if MK3(1,Colonne)>Tsai_Hill_max
        Tsai_Hill_max=MK3(1,Colonne);
        Pli_etude=Colonne;
        effort_min=effort;
        
    end
end

end