function [val] = f_obj(param)

global t;
global G1;

loc_t = t;

G_test = tf([0, param(1)], [param(2), 1], 'InputDelay', param(3));
%G_test = G_test * pade(param(3), 3)

% Calculate the step response of G1 and G_test
step_response_G1 = step(G1, t);
step_response_G_test = step(G_test, t);

% Calculate the RMS error point-wise
error = step_response_G1 - step_response_G_test;
rms_error = sqrt(mean(error.^2)) + max(abs(error .* loc_t'));

val = rms_error;
end