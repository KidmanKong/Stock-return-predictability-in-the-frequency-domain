function dmspe=DMSPE(predict,actual,theta)

% INPUTS:
%      predict   - a t by k matrix
%      actual    - a t by 1 matrix
% 
% OUTPUTS:
%      dmspe     - a t by 1 matrix

[t,k] = size(predict);
dmspe = nan(t,1);
a = nan(t,1);
for i = 1:t
    a(i) = theta^(i-1);
end
b = ((actual-predict).^2).*a;
c = nan(t,k);
c(1,:) = b(1,:);
for i = 2:t
    c(i,:) = c(i-1,:)+b(i,:);
end
d = c.^(-1);
w = nan(t,k);
for i = 1:t
    s = sum(d(i,:));
    for j = 1:k
        w(i,j) = d(i,j)/s;
    end
end
f = ones(1,k)/k;
w = [f;w(1:end-1,:)];
e = predict.*w;
for i = 1:t
    dmspe(i) = sum(e(i,:));
end        

