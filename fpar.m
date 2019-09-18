function [x_opt, f_x_opt] = fpar(c, a)
	
	x_opt = (-c(1:end-1) .* (a.^2)) ./ (2 * c(end));
	x_opt(end + 1) = sum(c(1:end-1).^2 .* a.^2) / (4 * c(end)^2);

	f_x_opt = c' * x_opt;

endfunction