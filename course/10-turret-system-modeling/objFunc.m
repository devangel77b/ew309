function metric = objFunc(x,dat)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

s = tf('s');
H = x(1)/((s+x(2))*(s+x(3)));
out = lsim(H,dat(:,1),dat(:,2));
metric = sqrt(mean((dat(:,3)-out).^2));

end

