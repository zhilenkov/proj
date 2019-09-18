%rand('state', 3)

%lst = [25, 50, 100, 150, 200, 250, 300];
lst = 500;
av_harm_means = zeros(1, length(lst));

t = 30;

for k = 1:t

harmonic_means = [];

for n = lst
%n = 50;

a = randi([1, 3], n, 1);
z = randi([4, 6], n, 1)

% z = zeros(n, 1);
% [mx imx] = max(a);
% z(imx) = mx + 1e-2;
% z(imx + 1) = 1e-2;

format long

[l_opt, x_opt] = fproj(z, a);

asq = a.^2;

f = @(l) sum(((z.^2) ./( ( (1 + ( l ./ asq) ).^2 ) .* asq))) - 1;

c = -z;

[x isopt] = flin(c, a);
S = [zeros(length(x), 1) x];

diff_opt = [];

m = 0;

for i = 1:50
	x_p = ptp(S - z, 1000, 1e-8, 0, [], []);
	c = x_p;
	x = flin(c, a);
	diff_opt = [diff_opt norm(x - x_opt)];
	S = [S x];
endfor

j = 0;
prod = 1;

for i = 1:30
	div_norm = norm(S(:, i+1) - x_opt, 'inf')/norm(S(:, i) - x_opt, 'inf');
	printf('norm(x_%d - x_opt) = %e, norm(x_%d - x_opt)/norm(x_%d - x_opt) = %e \n', i, norm(S(:, i) - x_opt, 'inf'), i + 1, i, div_norm)
	if div_norm < 1
		prod *= div_norm;
		j++;
	endif
endfor

printf('mean_harmonic = %e\n', prod^(1/j));

harmonic_means = [harmonic_means, prod^(1/j)];

% pkg load geometry;
% figure(1)
% t = [0:0.1:2*pi];
% x = 2 * sin(t);
% y = cos(t);
% plot(x, y, '-b');
% hold on;
% plot(S(1, 1:4), S(2, 1:4), 'or');
% plot(4, 4, '*g');
% plot(x_opt(1), x_opt(2), '*g');
% axis ([-3, 5, -3, 5]);
% hold off;

% figure(2)
% scatter(1:length(diff_opt), log(diff_opt))

endfor

av_harm_means += harmonic_means

endfor

av_harm_means *= 1/t;

i = 1
for n = lst
	printf('n = %d, av_harm_means = %e\n', n, av_harm_means(i));
	i++;
endfor
scatter(lst, av_harm_means)