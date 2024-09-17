trans_steps = 10;
num_p = 20;
initial_x = 0.1;
initial_y = 0.1;
a_k = 1.2;
b_k = 0.4;

x_old = initial_x;
y_old = initial_y;

for j = 1:trans_steps
    x_new = a_k - x_old.^2 + b_k * y_old;
    y_new = x_old;
    x_old = x_new;
    y_old = y_new;
end

x_p = zeros(num_p, 1);
y_p = zeros(num_p, 1);
x_p(1) = x_new;
y_p(1) = y_new;

for j = 1:num_p - 1
    x_p(j+1)=a_k - x_p(j)^2 + b_k * y_p(j);
    y_p(j+1)=x_p(j);
end

dist_mtrx = sparse(num_p, num_p);

for j = 1:num_p
    for i = j + 1:num_p
        dist = (x_p(i)-x_p(j))^2+(y_p(i)-y_p(j))^2;
        dist_mtrx(i,j)=dist;
    end
end

dist_mtrx = sqrt(dist_mtrx);

min_dist = min(min(dist_mtrx + (1000*dist_mtrx == 0)));
max_dist = max(max(dist_mtrx));
max_dist = 2^ceil(log(max_dist) / log(2));
num_div = floor(log(max_dist/min_dist)/log(2));
num_ranges = num_div + 1;
range_vector = max_dist * 2.^(-((1:num_ranges)' - 1));
num_pairs = num_p * (num_p - 1) / 2;

correl_values = zeros(num_ranges, 1);

for j = 1:num_ranges
    range = range_vector(j);
    num_pairs_within_range = sum(sum(dist_mtrx < range & dist_mtrx > 0));
    correl_values(j) = num_pairs_within_range / num_pairs;
end

correl_integral_val = sum(correl_values) / num_ranges;
disp(['Correlation integral value = ', num2str(correl_integral_val)]);
figure(2)
plot(range_vector, correl_values, 'o-');
xlabel('r');
ylabel('C(r)');
grid on


discard = 3;
n1 = discard + 1;
n2 = num_ranges - discard;
inside_range = n1:n2;
log_range = log(range_vector) / log(2);
log_correl_values = log(correl_values) / log(2);
log_inside_range = log(inside_range);
log_correl_inside = log_correl_values(inside_range);

koef = polyfit(log_inside_range, log_correl_inside, 1);
fractal_dimens = koef(1);
fitted_values = fractal_dimens * log_range + koef(2);

figure(1);
plot(log_range,log_correl_values, 'o-');
hold on
plot(log_range, fitted_values, 'r-');
axis tight
plot([log_range(n1), log_range(n1), [-30 30], 'k--']);
plot([log_range(n2), log_range(n2), [-30 30], 'k--']);
xlabel('log_2(r)');
ylabel('log_2(C(r))');
title(['D_c = ', num2str(fractal_dimens)]);
gridon

