
lst = 3;

mean_conv = zeros(1, length(lst));

for n = lst

format long;

a = randi([1, 2], n - 1, 1)
z = randi([100, 150], n, 1)

[x_opt, l, f_opt] = gproj(z, a)

asq = a.^2;

f = @(l) sum(( z(1:end-1).^2 ./ ( ( (1 + ( l ./ asq) ) ).^2 .* asq))) - l/2 - z(end);

x0 = zeros(n, 1);
x0(end) = z(end) + 1;
c = x0 - z

x = fpar(c, a);
S = [x0 x];

diff_opt = [];

m = 0;

for i = 1:100
	x_p = ptp(S - z, 1000, 1e-7, 0, [], []);
	c = x_p;
	x = fpar(c, a);
	diff_opt = [diff_opt norm(x - x_opt)];
	S = [S x];
endfor

conver = [];

for i = 1:50
	div_norm = norm(S(:, i+1) - x_opt)/norm(S(:, i) - x_opt);
	% printf('norm(x_%d - x_opt) = %e, norm(x_%d - x_opt)/norm(x_%d - x_opt) = %e \n', i, norm(S(:, i) - x_opt), i + 1, i, div_norm)
	printf('norm(x_%d - x_opt) = %e, norm(x_%d - x_opt)/norm(x_%d - x_opt) = %e \n', i, norm(S(:, i) - x_opt), i + 1, i, div_norm)
	conver = [conver; i, norm(S(:, i) - x_opt)];
endfor

endfor

plot(conver(:, 1), log(conver(:, 2)), 'or')