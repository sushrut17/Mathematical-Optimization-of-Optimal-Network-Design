%code for finding convex hull

function [Y_c,feed] = convhullfunc2(Y,p1,p2,l)


[h,Area] = convhull(Y(:,p1),Y(:,p2));
temp = -1;
% l = 5;
feed = [];

for i = 1:length(h)
    Y_c(i,1:l) = Y(h(i),1:l);
end

feed = [];
for i = 1:length(h)-1
    if abs(h(i) - h(i+1)) ~= 1
%         for j =1:4
%             Y_c(h(i):h(i+1),j) = [linspace(Y(h(i),j),Y(h(i+1),j),h(i+1)-h(i)+1)];
%         end
        temp = temp + 2;
        feed = [feed;Y(h(i),1:l);Y(h(i+1),1:l)];
%         feed = [feed;Y(h(i):h(i+1),1:l)];
%         
       
    end
end


end