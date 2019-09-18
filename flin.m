function [x, isopt] = flin(c,  a)

	asq = a.^2;
	csq = c.^2;
	
	l = sqrt(csq' * asq);
	
	x = -c .* asq / l;
	
	isopt = norm(c + l * x ./ asq)
	
endfunction