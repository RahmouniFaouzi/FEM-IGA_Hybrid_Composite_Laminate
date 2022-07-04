function [Tsai_Wu_max,Pli_etude,MK5,effort_min,failure_C]=Tsai_Wu_failure(Npli,Sig,MATP,K, theta,Tsai_Wu_max,Pli_etude,effort_min)

MK5=zeros(1,Npli);
failure_C  = 'Tsai_Wu_failure';
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
    
    F1=((1/K(MATP(Colonne),1))-(1/K(MATP(Colonne),2)))*Sig_mat(1,Colonne);
    F2=((1/K(MATP(Colonne),3))-(1/K(MATP(Colonne),4)))*Sig_mat(2,Colonne);
    F4=(1/(K(MATP(Colonne),1)*K(MATP(Colonne),2)))*Sig_mat(1,Colonne)^2;
    F5=(1/(K(MATP(Colonne),3)*K(MATP(Colonne),4)))*Sig_mat(2,Colonne)^2;
    F6=(1/(K(MATP(Colonne),5)^2))*Sig_mat(3,Colonne)^2;
    F7=-(1/(K(MATP(Colonne),1)*K(MATP(Colonne),2)*K(MATP(Colonne),3)*K(MATP(Colonne),4))^0.5)*Sig_mat(1,Colonne)*Sig_mat(2,Colonne);
    MK5(1,Colonne)=F1+F2+F4+F5+F6+F7;
    
    a=F4+F5+F6+F7;
    b=F1+F2;
    c=-1;
    
    delta=b^2-4*a*c;
    if delta>0
        r1=(-b-sqrt(delta))/(2*a);
        r2=(-b+sqrt(delta))/(2*a);
        r=max(r1,r2);
    elseif delta==0
        r=-b/(2*a);
    elseif delta<0
        disp('No reel solution...')
    end
    
    if effort_min>r
        Tsai_Wu_max=MK5(1,Colonne);
        Pli_etude=Colonne;
        effort_min=r;
    end
    
end