N = 300;
tvec = linspace(-pi + 2*pi/N, pi, N);
rvec = 3 + cos(4.*tvec + pi);
rprimvec = -4*sin(4.*tvec + pi);
rbisvec = -16*cos(4.*tvec + pi);
y1 = rvec .* cos(tvec);
y2 = rvec .* sin(tvec);
nu1 = rvec .* cos(tvec) + rprimvec .* sin(tvec);
nu2 = rvec .* sin(tvec) - rprimvec .* cos(tvec);
nu1 = nu1 ./ sqrt(rvec.^2 + rprimvec.^2);
nu2 = nu2 ./ sqrt(rvec.^2 + rprimvec.^2);
vecdsdt = sqrt(rprimvec.^2 + rvec.^2);

index = 1;
steps = 20001;
iterations = linspace(0,2,steps);
cond_vals = zeros(1,steps); % Store condition numbers

for k = iterations
    A_k = zeros(N, N);
    % Compute A_k
    for i = 1:N
        for j = 1:N
                nu_i = [nu1(i), nu2(i)];
                r_j = [y1(j), y2(j)];
                r_i = [y1(i), y2(i)];
                difference = r_i - r_j;
                hankel_func = besselh(1, 1, k*norm(difference));
                taljare = (1i*k/4)*hankel_func;
                namnare = norm(difference);
                A_k(i, j) = dot(nu_i, difference) * taljare / namnare;
        end
    end

    % Diagonal elements of A_k
    for i = 1:N
        taljare = rprimvec(i)^2 - 0.5*rbisvec(i)*rvec(i) + 0.5*rvec(i)^2;
        namnare = 2*pi*(rprimvec(i)^2 + rvec(i)^2)^(3/2);
        A_k(i, i) = taljare / namnare;
    end

    % Compute inverse_matrix
    Kmat = (-eye(N)/2 + 2*pi/N * A_k * diag(vecdsdt));

    % Compute and store the condition number
    cond_vals(index) = log10(cond(Kmat));
    index = index + 1;
end
%%
% Plot the condition number as a function of k
% cond_values = log10(cond_values);

figure;
plot(times_to_loop, cond_values, 'LineWidth', 0.6);
xlabel('k-value');
ylabel('Condition Number');
title('logrithmic plot Condition Number for each k');
grid on;
