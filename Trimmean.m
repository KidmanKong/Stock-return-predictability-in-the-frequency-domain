function y = Trimmean(x)
s = sort(x);
t = s(2:end-1);
y = sum(t)/(length(x)-2);
end

