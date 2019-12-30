function [initial_rate, Ubrate, output_t] = calcUbrate(y, tpoint)
% clcUbrate takes output of simulation and calculate
% percentile rate/ inital rate

Ubconc = y(:, 1);

% percentile rate
dt = tpoint(2:end) - tpoint(1:end-1);
dy = Ubconc(2:end) - Ubconc(1:end-1);
Ubrate = -dy./dt;
output_t = tpoint(1:end-1);

% t is when Ub concentration drops to midway (50%) between 
% beginning and the first minimum concentration
% if t is shorter than 10 min, t is set to 10 min. This is in line with the 
% initial reaction rate calculation from FRET measurement data
final_conc = 0.5*(Ubconc(1)+Ubconc(end));
for j = 1:length(Ubconc)
    if Ubconc(j) < final_conc
        break;
    end
end

for t = 1:length(tpoint)
    if t > 10*60
        break;
    end
end

j = max(t, j);

initial_rate = (Ubconc(1) - Ubconc(j))/tpoint(j);

end

