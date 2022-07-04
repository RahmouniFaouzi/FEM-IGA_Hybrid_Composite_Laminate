function [Hoffman_max,Pli_etude,effort_min,failure_C] = Stress_HoffmanISO( Nlyer,STR,...
    strength,theta,CLS)
	
Hoffman_max=0;
effort_min=9999999;
Pli_etude = 0;
for J = 1 : size (STR,1)
    count = 1;
    for e = 1 : Nlyer
        Str_G(:,count) = STR(J,:,e)';
        count = count + 1;
    end
    
    % Failures Criteria
        [Hoffman_max,Pli_etude,~,effort_min, failure_C] = Hoffman_failureISO(Nlyer,CLS,Str_G,strength, theta,Hoffman_max,Pli_etude,effort_min);
        
end % end Gauss point loop
end % end element loop