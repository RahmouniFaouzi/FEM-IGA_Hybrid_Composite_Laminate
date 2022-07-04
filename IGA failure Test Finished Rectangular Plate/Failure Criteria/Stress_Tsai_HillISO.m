function [Tsai_Hill_max,Pli_etude,effort_min,failure_C] = Stress_Tsai_HillISO( Nlyer,STR, ...
    strength,theta,CLS )

Tsai_Hill_max=0;
effort_min=0;
Pli_etude = 0;
MK3=zeros(1,Nlyer);

% cycle for element
for J = 1 : size (STR,1)
    count = 1;
    for e = 1 : Nlyer
        Str_G(:,count) = STR(J,:,e)';
        count = count + 1;
    end
    
    % Failures Criteria
    [Tsai_Hill_max,Pli_etude,MK3,effort_min, failure_C]=Tsai_Hill_failureISO(Nlyer,Str_G,CLS,strength, theta,Tsai_Hill_max,Pli_etude,MK3,effort_min);
	
      end % end Gauss point loop
end % end element loop
