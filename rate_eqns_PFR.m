%code for Attainable region for Plug Fow Reactor


% Van de Vusse (reversible)

function [f] = SLPcode(t,x)

k1 = 1;
k2 = 5;
k3 = 10;
k4 = 100;

f = zeros(4,1);

f(1) = -k1*x(1) + k2*x(2) - k4*x(1)^2;

f(2) = k1*x(1) - k2*x(2) - k3*x(2);

f(3) = k3*x(2);

f(4) = k4*x(1)^2;

end

% Van de Vusse (irreversible)

% function [f] = SLPcode(t,x)
% 
% k1 = 1;
% k2 = 2;
% k3 = 10;
% 
% f = zeros(4,1);
% 
% f(1) = -k1*x(1) - 2*k3*x(1)^2;
% 
% f(2) =  k1*x(1) - k2*x(2) ;
% 
% f(3) = k2*x(2);
% 
% f(4) = k3*x(1)^2;
% 
% end


%Denbigh

% function [f] = SLPcode(t,x)
% 
% k1 = 1;
% k2 = 0.1;
% k3 = 0.1;
% k4 = 0.1;
% 
% f = zeros(5,1);
% 
% f(1) = -k1*x(1)^2 - k3*x(1);
% 
% f(2) = k1*x(1)^2 - k2*x(2) - k4*x(2)^2;
% 
% f(3) = k2*x(2);
% 
% f(4) = k3*x(1);
% 
% f(5) = k4*x(2)^2;
% 
% end


%Trambouze

% function [f] = SLPcode(t,x)
% 
% k1 = 0.025;
% k2 = 0.2;
% k3 = 0.4;
% 
% f = zeros(4,1);
% 
% f(1) = -k1 - k2*x(1)^1 - k3*x(1)^2;
% 
% f(2) = k1;
% 
% f(3) = k2*x(1)^1;
% 
% f(4) = k3*x(1)^2;
% 
% end

%Pinene
% % function [f] = SLPcode(t,x)
% % 
% % k1 = 0.33384;
% % k2 = 0.26687;
% % k3 = 0.14940;
% % k4 = 0.18957;
% % k5 = 0.009598;
% % k6 = 0.29425;
% % k7 = 0.011932;
% % 
% % f = zeros(5,1);
% % 
% % f(1) = -k1*x(1) - k2*x(1) - 2*k5*x(1)^2;
% % 
% % f(2) = -k6*x(2) + k3*x(4);
% % 
% % f(3) = -k7*x(3) + k5*x(1)^2 + k4*x(4)^2;
% % 
% % f(4) = -k3*x(4) - 2*k4*x(4)^2 + k2*x(1) + 2*k7*x(3) + k6*x(2);
% % 
% % f(5) = k1*x(1);
% % 
% % end

%  4 independent reactions

% function [f] = SLPcode(t,x)
% 
% k1 = 0.025;
% k2 = 0.2;
% k3 = 0.4;
% k4 = 0.75;
% 
% f = zeros(5,1);
% 
% f(1) = -k1 - k2*x(1) - k3*x(1)^2 - k4*x(1);
% 
% f(2) = k1;
% 
% f(3) = k2*x(1);
% 
% f(4) = k3*x(1)^2;
% 
% f(5) = k4*x(1);
% end

%A->B->C
% function [f] = rate_eqns_PFR(t,x)
% 
% k1 = 1;
% k2 = 1;
% 
% f = zeros(3,1);
% 
% f(1) = -k1*x(1);
% 
% f(2) = k1*x(1) - k2*x(2);
% 
% f(3) = k2*x(2);
% 
% end

% A -> B
% function [f] = SLPcode(t,x)
% 
% k1 = 5;
% 
% f = zeros(2,1);
% 
% f(1) = -k1*x(1)^0 ;
% 
% f(2) = k1*x(1)^0 ;
% 
% end

%Denbigh (section 4.4)

% function [f] = SLPcode(t,x)
% 
% k1 = 1;
% k2 = 0.6;
% k3 = 0.6;
% k4 = 0.1;
% 
% f = zeros(5,1);
% 
% f(1) = -k1*x(1)^2 - k3*x(1);
% 
% f(2) = 0.5*k1*x(1)^2 - k2*x(2) - k4*x(2)^2;
% 
% f(3) = k2*x(2);
% 
% f(4) = k3*x(1);
% 
% f(5) = k4*x(2)^2;
% 
% end
% 
