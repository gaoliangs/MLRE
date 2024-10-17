function [F] = treescore()
clear all

objective = @(x) myObjective(x);
x0 = ones(1,8);
lb = [zeros(1,4),ones(1,4)*-10];
ub = [ones(1,4),ones(1,4)*10];

options = optimoptions('fmincon', 'Display', 'off','MaxFunctionEvaluations',50000,'MaxIterations',10000);
[x_optimal,fval,exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);

disp('Optimal solution t:');
disp(-log(x_optimal(1:4)));
disp('Optimal solution c:');
disp(exp(x_optimal(5:8)));
disp('Optimal objective function value:');
disp(exp(vpa(-fval)));
F= -fval;

end

function result = myObjective(x)
	variables
	lambda_formula
	lambda_formula_trivial
	suml=-exp(c1)*log(t1) - exp(c2)*log(t2) - exp(c3)*log(t3) - exp(c4)*log(t4) + 1-sum(la);
	F = sum(log(l))-length(l)*log(suml);
	result = -F;
end
