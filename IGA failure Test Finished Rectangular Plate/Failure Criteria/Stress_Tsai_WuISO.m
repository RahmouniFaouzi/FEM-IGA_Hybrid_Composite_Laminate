function [Tsai_Wu_max,Pli_etude,effort_min,failure_C] = Stress_Tsai_WuISO( Nlyer,STR,...
    strength,theta,CLS)

Tsai_Wu_max=0;
effort_min=99999;
Pli_etude = 0;

for J = 1 : size (STR,1)
    count = 1;
    for e = 1 : Nlyer
        Str_G(:,count) = STR(J,:,e)';
        count = count + 1;
    end
    % Failures Criteria
    [Tsai_Wu_max,  Pli_etude,~,effort_min, failure_C] = Tsai_Wu_failureISO (Nlyer,Str_G,CLS,strength, theta,Tsai_Wu_max,Pli_etude,effort_min);
     
end % end Gauss point loop
end % end element loop