function [rapport_max,Pli_etude,Mode,failure_C] = Stress_Max_StrISO( Nlyer,STR, ...
    strength,theta,CLS )

Str_G = zeros(3,Nlyer);
Strength_ratio =0;
Pli_etude =0;
Mode =0;
MK1 =zeros(3,Nlyer);
rapport_max =0;

% cycle for element
for J = 1 : size (STR,1)
    count = 1;
    for e = 1 : Nlyer
        Str_G(:,count) = STR(J,:,e)';
        count = count + 1;
    end
    
    % Failures Criteria
    [Strength_ratio,Pli_etude,Mode,~,rapport_max,failure_C] = Maximum_stress_failureISO...
        (Nlyer,CLS,Str_G,strength, theta, Strength_ratio ,Pli_etude,Mode,MK1,rapport_max);
    
end % end Gauss point loop
end % end element loop