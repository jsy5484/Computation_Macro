function [G0, G1] = linearize(func, Xss, deriv)
% This function linearize func with respect to two input vector around the
% Steady-State level, and return Jacobian matrices.

Fss = func(Xss, Xss);

n_states = length(Xss);

% Find G0 matrix by taking derivatives with respect to Xt
G0 = zeros(n_states,n_states);
for i = 1:n_states
	Xu = Xss;
	h = Xss(i) * deriv;
	if (h < deriv)
		h = deriv;
	end
	Xu(i) = Xu(i) + h;
	Fu = func(Xu, Xss);
	G0(:,i) = (Fu - Fss)./h;
end

% Find G1 matrix by taking derivatives with respect to Xt_1
G1 = zeros(n_states,n_states);
for i = 1:n_states
	Xu = Xss;
	h = Xss(i) * deriv;
	if (h < deriv)
		h = deriv;
	end
	Xu(i) = Xu(i) + h;
	Fu = func(Xss, Xu);
	G1(:,i) = (Fu - Fss)./h;
end
G1 = -G1;


end