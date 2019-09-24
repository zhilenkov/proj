function [x, l, f_opt] = gproj(z, a)

asq = a.^2;

f = @(l) sum(( z(1:end-1).^2 ./ ( ( (1 + ( l ./ asq) ) ).^2 .* asq))) - l/2 - z(end);

l_mn = 0;
l_mx = 1e8;

[l, fv_max, info, outp] = fzero(f, [l_mn, l_mx])

x = z(1:end-1) .* (1 + l * asq.^(-1)).^(-1);
x(end + 1) = z(end) + l/2;

if rows(x) < 2
	x = x';
endif

f_opt = (z - x)'*(z - x);

endfunction