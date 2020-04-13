%ORN7

clear all
clc
close all

tic;

% point = [0.2 0.3 0.4 0.5];
% initial feed condition
% the term 'ic' has been used for PFR feed points; 'feed' for CSTR feed points
ic = [1 0 0 0];
p = [1 2];
elem = 1;

l = length(ic);
% ------------------------------------------------------------------------
% STAGE I STARTS


%solving for PFR at intitial feed condition
options = odeset('NonNegative',1:l);
[TT,XX] = ode15s(@(t,x) rate_eqns_PFR(t,x),[0,20000],ic,options);

X = [];
T = [];
if elem == 1
    for i = 1:size(XX,1)
        %     if ((XX(i,1) + XX(i,2) + XX(i,3) + XX(i,4) + XX(i,5)) >= 0.99) && ((XX(i,1) + XX(i,2) + XX(i,3) + XX(i,4) + XX(i,5)) <= 1.01) && (XX(i,1) >= 0) && (XX(i,2) >= 0) && (XX(i,3) >= 0) && (XX(i,4) >= 0) && (XX(i,5) >= 0)
        if sum(XX(i,:)) >= 0.99*sum(ic(1,:)) && sum(XX(i,:)) <= 1.01*sum(ic(1,:)) && (XX(i,1) >= 0) && (XX(i,2) >= 0) && (XX(i,3) >= 0)  && (XX(i,4) >= 0) 
            X = [X;XX(i,:)];
            T = [T;TT(i,:)];
            %         X(i,2) = XX(i,2);
            %         X(i,3) = XX(i,3);
            %         X(i,4) = XX(i,4);
            %         X(i,5) = XX(i,5);
        end
    end
else
    X = XX;
    T = TT;
end

%transfering X into Y just to maintain consistency with future nomenclature
% Y is the concentration vector
Y(:,1:l) = X(:,1:l);
Zn = Y(:,1:l);

%plot Y
figure(1)
hold on
plot(Y(:,p(1)),Y(:,p(2)),'b')

% if condition to check whether Y is convex
% convhull gives an error if the curve is already convex
% convhullfunc is a function which gives the convexified points and also gives the feed conditions for future AR construction stages (Here, feed2)
% Yh is the nomenclature for convexified Y

% if sum(Y(1:size(Y,1),1)) ~= 0
%     [Yh,feed2] = convhullfunc(Y);
% end

Z = [];
feed2= [];
if sum(abs(diff(Y(:,p(1))))) ~= 0       %a check on errors
    %     h = convhull(Y(:,p(1)),Y(:,p(2)));
    %     hunique = unique(h);
    %         if sum(diff(hunique)) ~= length(diff(hunique))  %a check on errors
    [Yh,feed2] = convhullfunc2(Y,p(1),p(2),l);
    
    
    if isempty(feed2) == 1
        return
    end
end


Z = [Z ;Yh];

figure(1)
plot(Yh(:,p(1)),Yh(:,p(2)),'g')

% end

Oc1 = [X; X(1,:)];

if isempty(feed2) == 1
    return
end

%plot convexified Y
% figure(1)
% plot(Yh(:,1),Yh(:,2),'g')

[hz1,A1] = convhull(Z(:,p(1)),Z(:,p(2)));

for i = 1:length(hz1)
    Zh1(i,1:l) = Z(hz1(i),1:l);
end
figure(1)
plot(Zh1(:,p(1)),Zh1(:,p(2)),'k')

figure(2)
hold on
plot(Zh1(:,p(1)),Zh1(:,p(2)),'b')

length(Z)
% Area = A1;

%Y is normal curve
%Yh is convexified curve

X1 = X;
Y1 = Yh;

% STAGE I ENDS
% ------------------------------------------------------------------------


% STAGE II STARTS

% range of residence time
% tau = 0:0.02:20;
% tau = [tau 20000];

tau = 0:0.0001:1;
tau = [tau 1:0.02:20 20000];

%trambouze
% tau = 0:0.001:10;
% tau = [tau 10:0.02:20:100 20000];

feed2 = unique(feed2,'rows');

ic3 = [];
ic3record = [];

for j = 1:size(feed2,1)
    for i = 1:length(tau)
        options = optimset('Display','off');
        YY2(i,1:l) = fsolve(@(x) rate_eqns_CSTR(x,tau(i),feed2(j,:)),zeros(1,l),options);
    end
    %     Y2(end+1,1:4) = [0 0 0.9 0.1];
    
    YYY2 = [];
    TTT2 = [];
    for i = 1:size(YY2,1)
        %         if ((YY2(i,1) + YY2(i,2) + YY2(i,3) + YY2(i,4)  + YY2(i,5)) >= 0.99) && ((YY2(i,1) + YY2(i,2) + YY2(i,3) + YY2(i,4)+ YY2(i,5)) <= 1.01) && (YY2(i,1) >= 0) && (YY2(i,2) >= 0) && (YY2(i,3) >= 0) && (YY2(i,4) >= 0) && (YY2(i,5) >= 0)
        if (YY2(i,1) >= 0) && (YY2(i,2) >= 0) && (YY2(i,3) >= 0) && (YY2(i,4) >= 0) 
            YYY2 = [YYY2;YY2(i,:)];
            %             Y2(i,2) = YY2(i,2);
            %             Y2(i,3) = YY2(i,3);
            %             Y2(i,4) = YY2(i,4);
            %                         Y2(i,5) = YY2(i,5);
        end
    end
    
    Y2 = [];
    T2 = [];
    if elem == 1
        for i = 1:size(YYY2,1)
            %         if ((YY2(i,1) + YY2(i,2) + YY2(i,3) + YY2(i,4)  + YY2(i,5)) >= 0.99) && ((YY2(i,1) + YY2(i,2) + YY2(i,3) + YY2(i,4)+ YY2(i,5)) <= 1.01) && (YY2(i,1) >= 0) && (YY2(i,2) >= 0) && (YY2(i,3) >= 0) && (YY2(i,4) >= 0) && (YY2(i,5) >= 0)
            if sum(YYY2(i,:)) >= 0.99*sum(ic(1,:)) && sum(YYY2(i,:)) <= 1.01*sum(ic(1,:))
                Y2 = [Y2;YYY2(i,:)];
                T2 = [T2;tau(i)];
                %             Y2(i,2) = YY2(i,2);
                %             Y2(i,3) = YY2(i,3);
                %             Y2(i,4) = YY2(i,4);
                %                         Y2(i,5) = YY2(i,5);
            end
        end
    else
        Y2 = YYY2;
        T2 = TTT2;
    end
    
    Zn = [Zn;Y2];
    On2{j} = Y2;
    Onn2{j} = YY2;
    
    %plotting the CSTR ARc's
    figure(1)
    plot(Y2(:,p(1)),Y2(:,p(2)),'b')
    
    if sum(abs(diff(Y2(:,p(1))))) ~= 0 && sum(abs(diff(Y2(:,p(2))))) ~= 0
        h2 = convhull(Y2(:,p(1)),Y2(:,p(2)));
        hunique2 = unique(h2);
        %                 if sum(diff(hunique2)) ~= length(diff(hunique2))
        [Y2h,ic3] = convhullfunc2(Y2,p(1),p(2),l);
        
        if isempty(ic3) == 0
            ic3record = ic3;
        end
        
        %the code below cancels any repetition in storage of intial
        %conditions
        ss2 = [];
        for s2 = 1:size(ic3,1)
            for t2 = 1:size(ic,1)
                if ic3(s2,1:l) == ic(t2,1:l)
                    ss2 = [ss2 s2];
                end
            end
        end
        ic3(ss2,:) = [];
        
        
        Z = [Z ;Y2h];
        Oh2{j} = Y2h;
        Oc2{j} = [On2{j} ; On2{j}(1,:)];
        
        figure(1)
        plot(Y2h(:,p(1)),Y2h(:,p(2)),'g')
        
        %                 end
        
        if isempty(ic3) == 0
            P{j} = ic3;
        end
        
        if isempty(ic3record) == 0
            Precord{j} = ic3record;
        end
    end
    
    % P stores the ic3 for future reference
    
    
    % plotting convexified CSTR curves
    %     figure(1)
    %     plot(Y2h(:,1),Y2h(:,2),'g')
    
    %     Area2(j)= area;
    
    clear ic3;
    clear ic3record
end

clear ic3
ic3 = P{1,1};
for i = 2:length(P)
    ic3 = [ic3 ;P{1,i}];
end

if isempty(ic3) == 1
    [hz2,A2] = convhull(Z(:,p(1)),Z(:,p(2)));
    for i = 1:length(hz2)
        Zh4(i,1:l) = Z(hz2(i),1:l);
    end
    return
end

[hz2,A2] = convhull(Z(:,p(1)),Z(:,p(2)));

for i = 1:length(hz2)
    Zh2(i,1:l) = Z(hz2(i),1:l);
end
figure(1)
plot(Zh2(:,p(1)),Zh2(:,p(2)),'k')

figure(2)
plot(Zh2(:,p(1)),Zh2(:,p(2)),'g')

length(Z)


X2 = On2;
Y2 = Oh2;

% STAGE II ENDS
% ------------------------------------------------------------------------


% STAGE III STARTS

%regenerating ic3 from P
clear ic3
ic3 = P{1,1};
for i = 2:length(P)
    ic3 = [ic3 ;P{1,i}];
end

ic3record = Precord{1,1};
for i = 2:length(Precord)
    ic3record = [ic3record ;P{1,i}];
end

% further extending the AR by constructing PFRs at points of non-convexification on the various CSTRs
feed4 = [];
feed4record = [];

ic3 = unique(ic3,'rows');

for j = 1:size(ic3,1)
    
    options = odeset('NonNegative',1:l);
    [TT3{j},XX3] = ode15s(@(t,x) rate_eqns_PFR(t,x),[0,20000],ic3(j,:),options);
    
    X3 = [];
    T3{j} = [];
    if elem == 1
        for i = 1:size(XX3,1)
            %         if ((XX3(i,1) + XX3(i,2) + XX3(i,3) + XX3(i,4)  + XX3(i,5)) >= 0.99) && ((XX3(i,1) + XX3(i,2) + XX3(i,3) + XX3(i,4)) +  XX3(i,5)) <= 1.01 && (XX3(i,1) >= 0) && (XX3(i,2) >= 0) && (XX3(i,3) >= 0) && (XX3(i,4) >= 0) && (XX3(i,5) >= 0)
            if sum(XX3(i,:)) >= 0.99*sum(ic(1,:)) && sum(XX3(i,:)) <= 1.01*sum(ic(1,:)) && (XX3(i,1) >= 0) && (XX3(i,2) >= 0) && (XX3(i,3) >= 0) && (XX3(i,4) >= 0) 
                X3 = [X3;XX3(i,:)];
                T3{j} = [T3{j};TT3{j}(i)];
                %             X3(i,2) = XX3(i,2);
                %             X3(i,3) = XX3(i,3);
                %             X3(i,4) = XX3(i,4);
                %                         X3(i,5) = XX3(i,5);
            end
        end
    else
        X3 = XX3;
        T3{j} = TT3{j};
    end
    
    if isempty(X3) == 0
        
        Y3(:,1:l) = X3(:,1:l);
        Zn = [Zn;Y3];
        
        On3{j} = X3;
        Onn3{j} = XX3;
        Oc3{j} = [On3{j} ; On3{j}(1,:)];
        
        figure(1)
        plot(Y3(:,p(1)),Y3(:,p(2)),'b')
        
        %convexification and generation of feed points for next stage construction
        if sum(abs(diff(Y3(:,p(1))))) ~= 0 && sum(abs(diff(Y3(:,p(2))))) ~= 0
            Y3unique = unique(Y3(:,1));
            if size(Y3unique,p(1)) ~= 2
                h3 = convhull(Y3(:,p(1)),Y3(:,p(2)));
                hunique3 = unique(h3);
                %         if sum(diff(hunique3)) ~= length(diff(hunique3))
                [Y3h,feed4] = convhullfunc2(Y3,p(1),p(2),l);
                feed4record = feed4;
                ss3 = [];
                for s3 = 1:size(feed4,1)
                    for t3 = 1:size(feed2,1)
                        if feed4(s3,1:l) == feed2(t2,1:l)
                            ss3 = [ss3 s3];
                        end
                    end
                end
                feed4(ss3,:) = [];
                
                Z = [Z ;Y3h];
                
                figure(1)
                plot(Y3h(:,p(1)),Y3h(:,p(2)),'g')
                
                % Q stores feed4 for future reference
                Q{j} = feed4;
                Qrecord{j} = feed4record;
                Oh3{j} = Y3h;
                %     figure(1)
                %     plot(Y3h(:,1),Y3h(:,2),'g')
                %         end
                if isempty(feed4) == 0
                    P4{j} = feed4;
                end
                
                if isempty(feed4record) == 0
                    P4record{j} = feed4record;
                end
            end
            
            
            
            %clearing T3,X3,Y3 as they are of different sizes for different iterations since we have used [0,tspan] and not [0:tspan]
            %     clear T3
            
        end
    end
    clear X3
    clear Y3
    clear feed4
    clear feed4record
end

% regenerating feed4 from Q
feed4 = Q{1,1};
for i = 2:length(Q)
    feed4 = [feed4 ;Q{1,i}];
end


if isempty(feed4) == 1
    [hz3,A3] = convhull(Z(:,p(1)),Z(:,p(2)));
    for i = 1:length(hz3)
        Zh4(i,1:l) = Z(hz3(i),1:l);
    end
    return
end


[hz3,A3] = convhull(Z(:,p(1)),Z(:,p(2)));

for i = 1:length(hz3)
    Zh3(i,1:l) = Z(hz3(i),1:l);
end

figure(1)
plot(Zh3(:,p(1)),Zh3(:,p(2)),'k')

figure(2)
plot(Zh3(:,p(1)),Zh3(:,p(2)),'k')

length(Z)

X3 = On3;
Y3 = Oh3;

% STAGE III ENDS
% ------------------------------------------------------------------------

%
% % % % STAGE IV STARTS
% %
% % regenerating feed4 from Q
% feed4 = Q{1,1};
% for i = 2:length(Q)
%     feed4 = [feed4 ;Q{1,i}];
% end
% 
% feed4record = Qrecord{1,1};
% for i = 2:length(Qrecord)
%     feed4record = [feed4record ;Q{1,i}];
% end
% %
% % last stage; constructing CSTRs from points of non-convexification from the various PFRs
% 
% ic5 = [];
% ic5record = [];
% R = [];
% 
% for j = 1:size(feed4,1)
% 
%     for i = 1:length(tau)
%         options = optimset('Display','off');
%         YY4(i,1:l) = fsolve(@(x) rate_eqns_CSTR(x,tau(i),feed4(j,:)),zeros(1,l),options);
%     end
% 
%     YYY4 = [];
%     TTT4 = [];
%     for i = 1:size(YY4,1)
%         %         if ((YY2(i,1) + YY2(i,2) + YY2(i,3) + YY2(i,4)  + YY2(i,5)) >= 0.99) && ((YY2(i,1) + YY2(i,2) + YY2(i,3) + YY2(i,4)+ YY2(i,5)) <= 1.01) && (YY2(i,1) >= 0) && (YY2(i,2) >= 0) && (YY2(i,3) >= 0) && (YY2(i,4) >= 0) && (YY2(i,5) >= 0)
%         if (YY4(i,1) >= 0) && (YY4(i,2) >= 0) && (YY4(i,3) >= 0) && (YY4(i,4) >= 0) 
%             YYY4 = [YYY4;YY4(i,:)];
%             %             Y2(i,2) = YY2(i,2);
%             %             Y2(i,3) = YY2(i,3);
%             %             Y2(i,4) = YY2(i,4);
%             %                         Y2(i,5) = YY2(i,5);
%         end
%     end
% 
%     Y4 = [];
%     T4 = [];
%       if elem == 1
%         for i = 1:size(YYY4,1)
%             %         if ((YY2(i,1) + YY2(i,2) + YY2(i,3) + YY2(i,4)  + YY2(i,5)) >= 0.99) && ((YY2(i,1) + YY2(i,2) + YY2(i,3) + YY2(i,4)+ YY2(i,5)) <= 1.01) && (YY2(i,1) >= 0) && (YY2(i,2) >= 0) && (YY2(i,3) >= 0) && (YY2(i,4) >= 0) && (YY2(i,5) >= 0)
%             if sum(YYY4(i,:)) >= 0.99*sum(ic(1,:)) && sum(YYY4(i,:)) <= 1.01*sum(ic(1,:))
%                 Y4 = [Y4;YYY4(i,:)];
%                 T4 = [T4;tau(i)];
%                 %             Y2(i,2) = YY2(i,2);
%                 %             Y2(i,3) = YY2(i,3);
%                 %             Y2(i,4) = YY2(i,4);
%                 %                         Y2(i,5) = YY2(i,5);
%             end
%         end
%     else
%         Y4 = YYY4;
%         T4 = TTT4;
%     end
% 
% 
%     Zn = [Zn;Y4];
%     On4{j} = Y4;
%     Onn4{j} = YY4;
% 
%     figure(1)
%     plot(Y4(:,p(1)),Y4(:,p(2)),'b')
% 
%     %     ic5 = [];
% 
%     %     if sum(Y4(1:size(Y,1),1)) ~= 0
%     %     if diff(Y4(:,1)) ~= zeros(size(Y4,1)-1,1)
%     %         h4 = convhull(Y4(:,1),Y4(:,2));
%     %         if diff(h4) ~= ones(size(diff(h4),1),1)
%     if (sum(abs(diff(Y4(:,p(1))))) ~= 0 && sum(abs(diff(Y4(:,p(2))))) ~= 0) && size(Y4,1) > 2
%         h4 = convhull(Y4(:,p(1)),Y4(:,p(2)));
%         hunique4 = unique(h4);
% %         if sum(diff(hunique4)) ~= length(diff(hunique4))
%             [Y4h,ic5] = convhullfunc2(Y4,p(1),p(2),l);
% 
%             if isempty(ic5) == 0
%                 ic5record = ic5;
%             end
% 
%             ss4 = [];
%             for s4 = 1:size(ic5,1)
%                 for t4 = 1:size(ic3,1)
%                     if ic5(s4,1:l) == ic3(t4,1:l)
%                         ss4 = [ss4 s4];
%                     end
%                 end
%             end
%             ic5(ss4,:) = [];
% 
%             Z = [Z ;Y4h];
%             Oh4{j} = Y4h;
%             Oc4{j} = [On4{j} ; On4{j}(1,:)];
% 
%             figure(1)
%             plot(Y4h(:,p(1)),Y4h(:,p(2)),'g')
% 
%             if isempty(ic5) == 0
%             R{j} = ic5;
%         end
%         
%         if isempty(ic5record) == 0
%             Rrecord{j} = ic5record;
%         end
% 
%         end
% %         end
% 
% clear ic5;
% clear ic5record
%     end
% 
% 
% ic5 = [];
% if isempty(R) ~= 1
% ic5 = R{1,1};
% for i = 2:length(R)
%     ic5 = [ic5 ;R{1,i}];
% end
% end
% 
% if isempty(ic5) == 1
%     [hz4,A4] = convhull(Z(:,p(1)),Z(:,p(2)));
%     for i = 1:length(hz2)
%         Zh4(i,1:l) = Z(hz4(i),1:l);
%     end
%     return
% end
% 
% [hz4,A4] = convhull(Z(:,p(1)),Z(:,p(2)));
% 
% for i = 1:length(hz4)
%     Zh4(i,1:l) = Z(hz4(i),1:l);
% end
% 
% figure(1)
% plot(Zh4(:,p(1)),Zh4(:,p(2)),'k')
% 
% figure(2)
% plot(Zh4(:,p(1)),Zh4(:,p(2)),'r')
% 
% length(Z)
% 
% % STAGE IV ENDS
% 
% %------------------------------------------------------------------------
% 
% % STAGE V STARTS
% 
% %regenerating ic5 from P
% clear ic5
% ic5 = R{1,1};
% for i = 2:length(R)
%     ic5 = [ic5 ;R{1,i}];
% end
% 
% ic5record = Rrecord{1,1};
% for i = 2:length(Rrecord)
%     ic5record = [ic5record ;R{1,i}];
% end
% 
% % further extending the AR by constructing PFRs at points of non-convexification on the various CSTRs
% feed6 = [];
% feed6record = [];
% for j = 1:size(ic5,1)
%     
%     options = odeset('NonNegative',1:l);
%     [TT5{j},XX5] = ode15s(@(t,x) rate_eqns_PFR(t,x),[0,200],ic5(j,:),options);
%     
%     X5 = [];
%     T5{j} = [];
%     if elem == 1
%         for i = 1:size(XX5,1)
%             %         if ((XX3(i,1) + XX3(i,2) + XX3(i,3) + XX3(i,4)  + XX3(i,5)) >= 0.99) && ((XX3(i,1) + XX3(i,2) + XX3(i,3) + XX3(i,4)) +  XX3(i,5)) <= 1.01 && (XX3(i,1) >= 0) && (XX3(i,2) >= 0) && (XX3(i,3) >= 0) && (XX3(i,4) >= 0) && (XX3(i,5) >= 0)
%             if sum(XX5(i,:)) >= 0.99*sum(ic(1,:)) && sum(XX5(i,:)) <= 1.01*sum(ic(1,:)) && (XX5(i,1) >= 0) && (XX5(i,2) >= 0) && (XX5(i,3) >= 0)  
%                 X5 = [X5;XX5(i,:)];
%                 T5{j} = [T5{j};TT5{j}(i)];
%                 
%             end
%         end
%     else
%         X5 = XX5;
%         T5{j} = TT5{j};
%     end
%     
%     if isempty(X5) == 0
%         
%         Y5(:,1:l) = X5(:,1:l);
%         Zn = [Zn;Y5];
%         
%         On5{j} = X5;
%         Onn5{j} = XX5;
%         Oc5{j} = [On5{j} ; On5{j}(1,:)];
%         
%         figure(1)
%         plot(Y5(:,p(1)),Y5(:,p(2)),'b')
%         
%         %convexification and generation of feed points for next stage construction
%         if sum(abs(diff(Y5(:,p(1))))) ~= 0 && sum(abs(diff(Y5(:,p(2))))) ~= 0
%             Y5unique = unique(Y5(:,1));
%             if size(Y5unique,p(1)) ~= 2
%                 h5 = convhull(Y5(:,p(1)),Y5(:,p(2)));
%                 hunique3 = unique(h5);
%                 %         if sum(diff(hunique3)) ~= length(diff(hunique3))
%                 [Y5h,feed6] = convhullfunc2(Y5,p(1),p(2),l);
%                 feed6record = feed6;
%                 ss5 = [];
%                 for s5 = 1:size(feed6,1)
%                     for t5 = 1:size(feed4,1)
%                         if feed6(s5,1:l) == feed4(t5,1:l)
%                             ss5 = [ss5 s5];
%                         end
%                     end
%                 end
%                 feed6(ss5,:) = [];
%                 
%                 Z = [Z ;Y3h];
%                 
%                 figure(1)
%                 plot(Y5h(:,p(1)),Y5h(:,p(2)),'g')
%                 
%                 % Q stores feed4 for future reference
%                 S{j} = feed6;
%                 Srecord{j} = feed6record;
%                 Oh5{j} = Y5h;
%                 %     figure(15
%                 %     plot(Y3h(:,1),Y3h(:,2),'g')
%                 %         end
%                 if isempty(feed6) == 0
%                     R4{j} = feed6;
%                 end
%                 
%                 if isempty(feed6record) == 0
%                     R4record{j} = feed6record;
%                 end
%             end
%             
%             
%             
%             %clearing T3,X3,Y3 as they are of different sizes for different iterations since we have used [0,tspan] and not [0:tspan]
%             %     clear T3
%             
%         end
%     end
%     clear X5
%     clear Y5
%     clear feed6
%     clear feed6record
% end
% 
% % regenerating feed6 from S
% feed6 = S{1,1};
% for i = 2:length(S)
%     feed6 = [feed6 ;S{1,i}];
% end
% 
% 
% if isempty(feed6) == 1
%     [hz5,A5] = convhull(Z(:,p(1)),Z(:,p(2)));
%     for i = 1:length(hz5)
%         Zh6(i,1:l) = Z(hz5(i),1:l);
%     end
%     return
% end
% 
% 
% [hz5,A5] = convhull(Z(:,p(1)),Z(:,p(2)));
% 
% for i = 1:length(hz5)
%     Zh5(i,1:l) = Z(hz5(i),1:l);
% end
% 
% figure(1)
% plot(Zh5(:,p(1)),Zh5(:,p(2)),'k')
% 
% figure(2)
% plot(Zh5(:,p(1)),Zh5(:,p(2)),'k')
% 
% length(Z)
% 
% % STAGE V ENDS


toc;