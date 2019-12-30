function J = costFunction_rvc(logk0, data)
% costFunction_mcmc calculates the error between model and experiment

tspanmax = 6000;
k = 10.^logk0;
c = 1;
J_arr = zeros(length(data.Edim), 1);

for i = 1:length(data.Edim)
    Vary_dim = data.Edim(i);
    y = zeros(Vary_dim, 1);
    for j = 1:Vary_dim
        [yt, tt, ~] = simulation(data.C0(c, :), tspanmax, k);
        [y(j), ~, ~] = calcUbrate(yt, tt);
        c = c + 1;
    end
    y0 = data.ydata(:, i);
    y0 = y0(y0 ~= 0);
    J_arr(i) = 1/Vary_dim*sum((y - y0).^2)*1E19;
end

J = sum(J_arr);

end
