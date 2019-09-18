function [l_max, x_p, f_x, info, out_m, time] = fproj(z, a)

asq = a.^2;

f = @(l) sum(( z.^2 ./ ( ( (1 + ( l ./ asq) ) ).^2 .* asq))) - 1;

a_min = min(a);
a_max = max(a);

l_mn = 0;

y = z*min(a)/norm(z);
l_mx = (z - y)' * (z - y)

[l_max, fv_max, info_max, outp_max] = fzero(f, [l_mn, l_mx]);

x_p = z ./ (1 + l_max ./ asq);
f_x = (z - x_p)' * (z - x_p);
info = info_max;
out_m = outp_max;

endfunction