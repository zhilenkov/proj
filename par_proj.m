
lst = 500;

mean_conv = zeros(1, length(lst));

for n = lst

format long;

a = randi([1, 2], n - 1, 1);
z = randi([6, 7], n, 1);

% a = 1
% z = [6; 5]

[x_opt, l, f_opt] = gproj(z, a);

asq = a.^2;

f = @(l) sum(( z(1:end-1).^2 ./ ( ( (1 + ( l ./ asq) ) ).^2 .* asq))) - l/2 - z(end);

x0 = zeros(n, 1);
x0(end) = z(end) + 1;
c = x0 - z;

x = fpar(c, a);
S = [x0 x];

diff_opt = [];

m = 0;
dst = [];

for i = 1:50
	x_p = ptp(S - z, 1000, 1e-7, 0, [], []);
	c = x_p;
	dst = [dst, norm(c)];
	x = fpar(c, a);
	diff_opt = [diff_opt norm(x - x_opt)];
	S = [S x];
endfor

conver = [];

for i = 1:50
	div_norm = norm(S(:, i+1) - x_opt)/norm(S(:, i) - x_opt);
	% printf('norm(x_%d - x_opt) = %e, norm(x_%d - x_opt)/norm(x_%d - x_opt) = %e \n', i, norm(S(:, i) - x_opt), i + 1, i, div_norm)
	printf('norm(x_%d - x_opt) = %e\n', i, norm(S(:, i) - x_opt));
	conver = [conver; i, norm(S(:, i) - x_opt)];
endfor

endfor

%plot(conver(:, 1), log(conver(:, 2)), 'or')

% t = linspace(-7, 7, 100);
% tt = t .^ 2 / a.^2;

% plot(t, tt, '-r');
% hold on;
% plot(S(1, 1:20), S(2, 1:20), '-ob');
% plot(z(1), z(2), '*g')
% plot(x_opt(1), x_opt(2), '*g')
% hold off;

%printf('x(2) = x(1)^2 / a^2, a = %f\n', a);
%printf('projected point z = (%f, %f)\n', z(1), z(2));
%printf('optimal point x* = (%f, %f)\n', x_opt(1), x_opt(2));
% printf('parameters a = (')
% for i = 1:(length(a) - 1)
	% printf('%15.15f, ', a(i))
% endfor
% printf('%15.15f)\n', a(end))
% printf('projected point z = (')
% for i = 1:(length(z) - 1)
	% printf('%15.15f, ', z(i))
% endfor
% printf('%15.15f)\n', z(end))
% printf('optimal point x* = (')
% for i = 1:(length(x_opt) - 1)
	% printf('%15.15f, ', x_opt(i))
% endfor
% printf('%15.15f)\n', x_opt(end))
% x_last = S(:, end)'
% printf('dual l = %15.15f\n', l);
% printf('||z - x_opt|| = %15.15f\n', norm(z - x_opt));

% for i = 1:16
	% printf('x(%d) = (%6.15f, %6.15f)\n', i, S(1, i), S(2, i))
% endfor

function gg = grg(x, a)
	gg = 2*x(1:end-1) ./ a.^2;
	gg(end + 1) = -1;
	if rows(gg) < 2
		gg = gg';
	endif
endfunction

function gf = grf(x, z)
	gf = -2*(z - x);
endfunction

printf("norm(l*grad(g(x_opt)) + grad(f(x_opt))) = %e\n", norm(l*grg(x_opt, a) + grf(x_opt, z)))
for i = 1:length(dst)
	printf('iter = %d, dist(z, conv(S)) = %15.16f\n', i, dst(i));
endfor

rel = grg(x_opt, a) ./ c;
min(rel)
max(rel)