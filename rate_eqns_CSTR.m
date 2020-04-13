%code for Attainable region for CSTR


%Van de Vusse 

function [f] = SLPcode2(x,tau,feed)


k1 = 1;
k2 = 5;
k3 = 10;
k4 = 100;

f = zeros(4,1);

f(1) = feed(1) - x(1) + (-k1*x(1) + k2*x(2) - k4*x(1)^2)*tau;

f(2) = feed(2) - x(2) + (k1*x(1) - k2*x(2) - k3*x(2))*tau;

f(3) = feed(3) - x(3) + (k3*x(2))*tau;

f(4) = feed(4) - x(4) + (k4*x(1)^2)*tau;

end


%Van de Vusse (irreversible)
% function [f] = SLPcode2(x,tau,feed)
% 
% 
% k1 = 1;
% k2 = 2;
% k3 = 10;
% 
% f = zeros(4,1);
% 
% f(1) = feed(1) - x(1) + (-k1*x(1) - 2*k3*x(1)^2)*tau;
% 
% f(2) = feed(2) - x(2) + (k1*x(1) - k2*x(2))*tau;
% 
% f(3) = feed(3) - x(3) + (k2*x(2))*tau;
% 
% f(4) = feed(4) - x(4) + (k3*x(1)^2)*tau;
% 
% end



%Denbigh

% function [f] = SLPcode2(x,tau,feed)
% 
% 
% k1 = 1;
% k2 = 0.1;
% k3 = 0.1;
% k4 = 0.1;
% 
% f = zeros(5,1);
% 
% f(1) = feed(1) - x(1) + (-k1*x(1)^2 - k3*x(1))*tau;
% 
% f(2) = feed(2) - x(2) + (k1*x(1)^2 - k2*x(2) - k4*x(2)^2)*tau;
% 
% f(3) = feed(3) - x(3) + (k2*x(2))*tau;
% 
% f(4) = feed(4) - x(4) + (k3*x(1))*tau;
% 
% f(5) = feed(5) - x(5) + (k4*x(2)^2)*tau;
% 
% end


%Trambouze
% 
% function [f] = SLPcode2(x,tau,feed)
% 
% k1 = 0.025;
% k2 = 0.2;
% k3 = 0.4;
% 
% f = zeros(4,1);
% 
% f(1) = feed(1) - x(1) + (-k1*x(1)^0 - k2*x(1)^1 - k3*x(1)^2)*tau;
% 
% f(2) = feed(2) - x(2) + (k1*x(1)^0)*tau;
% 
% f(3) = feed(3) - x(3) + (k2*x(1)^1)*tau;
% 
% f(4) = feed(4) - x(4) + (k3*x(1)^2)*tau;
% 
% end

%zero order

% function[f] = SLPcode2(x,tau,feed)
% 
% k1 = 16;
% k2 = 8;
% 
% f = zeros(3,1);
% 
% f(1) = feed(1) - x(1) + (-k1)*tau;
% 
% f(2) = feed(2) - x(2) + (-k2+k1)*tau;
% 
% f(3) = feed(3) - x(3) + (k2)*tau;
% 
% end

% Pinene
% function [f] = SLPcode2(x,tau,feed)
% 
% k1 = 0.33384;
% k2 = 0.26687;
% k3 = 0.14940;
% k4 = 0.18957;
% k5 = 0.009598;
% k6 = 0.29425;
% k7 = 0.011932;
% 
% f = zeros(5,1);
% 
% f(1) = feed(1) - x(1) + (-k1*x(1) - k2*x(1) - 2*k5*x(1)^2)*tau;
% 
% f(2) = feed(2) - x(2) + (-k6*x(2) + k3*x(4))*tau;
% 
% f(3) = feed(3) - x(3) + (-k7*x(3) + k5*x(1)^2 + k4*x(4)^2)*tau;
% 
% f(4) = feed(4) - x(4) + (-k3*x(4) - 2*k4*x(4)^2 + k2*x(1) + 2*k7*x(3) + k6*x(2))*tau;
% 
% f(5) = feed(5) - x(5) + (k1*x(1))*tau;
% end

%4 independent reactions

% function [f] = SLPcode2(x,tau,feed)
% 
% k1 = 0.025;
% k2 = 0.2;
% k3 = 0.4;
% k4 = 0.75;
% 
% f = zeros(5,1);
% 
% f(1) = feed(1) - x(1) + (-k1*x(1)^0 - k2*x(1) - k3*x(1)^2 - k4*x(1))*tau;
% 
% f(2) = feed(2) - x(2) + (k1*x(1)^0)*tau;
% 
% f(3) = feed(3) - x(3) + (k2*x(1))*tau;
% 
% f(4) = feed(4) - x(4) + (k3*x(1)^2)*tau;
% 
% f(5) = feed(5) - x(5) + (k4*x(1))*tau;
% end


%A->B->C
% function [f] = SLPcode2(x,tau,feed)
% 
% 
% k1 = 1;
% k2 = 1;
% 
% f = zeros(3,1);
% 
% f(1) = feed(1) - x(1) + (-k1*x(1)^2)*tau;
% 
% f(2) = feed(2) - x(2) + (k1*x(1)^2 - k2*x(2))*tau;
% 
% f(3) = feed(3) - x(3) + (k2*x(2))*tau;
% 
% 
% end


%Denbigh (section 4.4)

% function [f] = SLPcode2(x,tau,feed)
% 
% 
% k1 = 1;
% k2 = 0.6;
% k3 = 0.6;
% k4 = 0.1;
% 
% f = zeros(5,1);
% 
% f(1) = feed(1) - x(1) + (-k1*x(1)^2 - k3*x(1))*tau;
% 
% f(2) = feed(2) - x(2) + (0.5*k1*x(1)^2 - k2*x(2) - k4*x(2)^2)*tau;
% 
% f(3) = feed(3) - x(3) + (k2*x(2))*tau;
% 
% f(4) = feed(4) - x(4) + (k3*x(1))*tau;
% 
% f(5) = feed(5) - x(5) + (k4*x(2)^2)*tau;
% 
% end
