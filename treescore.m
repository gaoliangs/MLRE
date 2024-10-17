function [F] = treescore()
clear all

objective = @(x) myObjective(x);
x0 = ones(1,62);
lb = [zeros(1,31),ones(1,31)*-10];
ub = [ones(1,31),ones(1,31)*10];

options = optimoptions('fmincon', 'Display', 'off','MaxFunctionEvaluations',50000,'MaxIterations',10000);
[x_optimal,fval,exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);

disp('Optimal solution t:');
disp(-log(x_optimal(1:31)));
disp('Optimal solution c:');
disp(exp(x_optimal(32:62)));
disp('Optimal objective function value:');
disp(exp(vpa(-fval)));
F= -fval;

end

function result = myObjective(x)
	variables
	lambda_formula
	lambda_formula_trivial
	suml=-exp(c1)*log(t1) - exp(c10)*log(t10) - exp(c11)*log(t11) - exp(c12)*log(t12) - exp(c13)*log(t13) - exp(c14)*log(t14) - exp(c15)*log(t15) - exp(c16)*log(t16) - exp(c17)*log(t17) - exp(c18)*log(t18) - exp(c19)*log(t19) - exp(c2)*log(t2) - exp(c20)*log(t20) - exp(c21)*log(t21) - exp(c22)*log(t22) - exp(c23)*log(t23) - exp(c24)*log(t24) - exp(c25)*log(t25) - exp(c26)*log(t26) - exp(c27)*log(t27) - exp(c28)*log(t28) - exp(c29)*log(t29) - exp(c3)*log(t3) - exp(c30)*log(t30) - exp(c31)*log(t31) - exp(c4)*log(t4) - exp(c5)*log(t5) - exp(c6)*log(t6) - exp(c7)*log(t7) - exp(c8)*log(t8) - exp(c9)*log(t9) + 1-sum(la);
	F = sum(log(l))-length(l)*log(suml);
	result = -F;
end
