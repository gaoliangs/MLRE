from sympy import symbols, Matrix, log, exp, Sum, Mul
from functools import reduce
import re
import sys
import os

current_directory = os.getcwd()

topo = sys.argv[1]
taxaorder = topo.replace("(", "").replace(")", "").split(',')


data = []
file_name = sys.argv[2]
input_file_path = os.path.join(current_directory, file_name)

with open(input_file_path, 'r') as mark:
   line = mark.readlines()
for i in range(len(line)):
   line1 = line[i][:-1]
   marker = line1.split(",")
   data.append(marker)



taxaname= data[0]
taxadata = data[1:]


#grouptaxa,singletaxa的taxa名
singletaxa = []
grouptaxa = []
for col_index in range(len(taxaname)):
   # 判断该列是否存在 'B'
   if any(row[col_index] == '01' for row in data):
      grouptaxa.append(taxaname[col_index])
   else:
      singletaxa.append(taxaname[col_index])


smloc = [taxaorder.index(i)+1 for i in singletaxa]
gm = [taxaname.index(i)+1 for i in grouptaxa]
gmloc = [taxaorder.index(i)+1 for i in grouptaxa]
#print(smloc)


markerdata = []
for mi in taxadata:
	orginial = topo
	for i,name in enumerate(taxaname):
		pattern = re.compile(r'\b' + re.escape(name) + r'\b')
		orginial = pattern.sub(mi[i], orginial, count=1)
	markerdata.append(orginial)
#print(markerdata)



def generate_trivial_marker(taxanum,gm):
    lists = []
    # 1. 全为0
    lists.append(['0'] * taxanum)
    # 2. 全为1
    lists.append(['1'] * taxanum)
    # 3. 只有一个1，其余全为0
    for i in range(taxanum):
        current_list = ['0'] * taxanum
        current_list[i] = '1'
        lists.append(current_list)
    # 4. 只有一个01，其余全为0，这个只有len(gm)个
    for i in range(len(gm)):
        current_list = ['0'] * taxanum
        current_list[gm[i]-1] = '01'
        lists.append(current_list)
    return lists

trivial_marker = generate_trivial_marker(len(taxaname), gm)

trivialtopo = []
for ti in trivial_marker:
   tritopo = topo
   for j,name in enumerate(taxaname):
      pattern = re.compile(r'\b' + re.escape(name) + r'\b')
      tritopo = pattern.sub(ti[j], tritopo, count=1)
   trivialtopo.append(tritopo)
#print(trivialtopo)



#number of edges
n_edge =topo.count('(')-1
#number of taxa
n_leaf = len(taxaname)

# 定义符号变量 x1 到 x10
t = symbols('t:{}'.format(n_edge+1))
c = symbols('c:{}'.format(n_edge+1))

# 定义符号变量 xl1 到 xl5
tl = symbols('tl:{}'.format(n_leaf+1))


# 叶节点的矩阵
matrix_leaf_group = Matrix([
    [0, 0, 0, 0],
    [0, 1, 0, 0],  
    [0, 0, 1, 0],
    [0, (1 - tl[0]) / 2, (1 - tl[0]) / 2, tl[0]]
])


# 内部节点的矩阵
matrix_internal = Matrix([
    [1, (-log(t[0])-(1-t[0]))*c[0]/2, (-log(t[0])-(1-t[0]))*c[0]/2,c[0]*(1-t[0])],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, (1-t[0])/2, (1-t[0])/2, t[0]]
])



matrix_leaf = Matrix([
    [0, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 1/2, 1/2, 1]
])

'''
matrix_leaf_group = Matrix([
    [0, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, (1 - exp(-tl[0]))/ 2, (1 - exp(-tl[0])) / 2, exp(-tl[0])]
])

matrix_internal = Matrix([
    [1, (t[0]-(1-exp(-t[0])))*c[0]/2, (t[0]-(1-exp(-t[0])))*c[0]/2,c[0]*(1-exp(-t[0]))],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, (1-exp(-t[0]))/2, (1-exp(-t[0]))/2, exp(-t[0])]
])
'''

def findsplit(s,edge_index,taxa_index):
	result = []
	stack = []
	edge_index = edge_index[0]
	taxa_index = taxa_index[0]
	current_part = ""

	for char in s[1:-1]:
		if char == ',' and not stack:
			string = current_part.strip()
			taxanum = len(string.replace('(','').replace(')','').split(','))
			edgenum = string.count('(')
			result.append([string,list(range(edge_index+1,edge_index+edgenum+1)),list(range(taxa_index,taxa_index+taxanum))])
			edge_index += edgenum
			taxa_index += taxanum
			current_part = ""
		else:
			current_part += char
			if char == '(':
				stack.append(char)
			elif char == ')':
				if stack:
					stack.pop()

	string = current_part.strip()
	taxanum = len(string.replace('(','').replace(')','').split(','))
	edgenum = string.count('(')
	result.append([string,list(range(edge_index+1,edge_index+edgenum+1)),list(range(taxa_index,taxa_index+taxanum))])

	return result




def combine(vector):
	v_not_invented = [i[0] for i in vector]
	v_lost= [i[1] for i in vector]
	v_fixed= [i[2] for i in vector]
	v_polymorphic = [i[3] for i in vector]

	v_l = 1 if all(element == 1 for element in v_lost) else 0
	v_f = 1 if all(element == 1 for element in v_fixed) else 0
	v_p = reduce(lambda x, y: x * y, v_polymorphic)

	if all(my_element == 0 for my_element in v_lost):
		v_n=0
	elif all(my_element == 1 for my_element in v_lost):
		v_n = reduce(lambda x, y: x + y, v_not_invented)
	else:
		introduce_index = [index for index, value in enumerate(v_lost) if value == 0]
		v_n = sum(v_not_invented[i] for i in introduce_index)

	return Matrix([v_n,v_l,v_f,v_p])


#print(combine([v1,v2]))


def recursive_split(s,edge_num,leaf_index):
	if '(' in s:
		parts = findsplit(s,edge_num,leaf_index)
		result = []
		for part in parts:
			if '(' not in part[0]:
				v = recursive_split(part[0],part[1],part[2])
				result.append(v)
				#print(result)
			else:
				#print('p',part)
				v = recursive_split(part[0],part[1],part[2])
				#print('s',v)
				e = part[1][0]
				current_matrix = matrix_internal.subs({t[0]: t[e]}).subs({c[0]: c[e]})
				#print('pv',current_matrix*v)
				result.append(current_matrix*v)

		vectors = combine(result)
		return vectors
	else:
		if leaf_index[0] in smloc:
			current_matrix = matrix_leaf.subs({tl[0]: tl[leaf_index[0]]})
		else:
			current_matrix = matrix_leaf_group.subs({tl[0]: tl[leaf_index[0]]})
		if s == '0':
			return current_matrix*Matrix([1,1,0,0])
		if s == '1':
			return current_matrix*Matrix([0,0,1,0])
		if s == '01':
			return current_matrix*Matrix([0,0,0,1])
		if s == '?':
			return Matrix([1,1,1,1])

#print(recursive_split('(((1,1),0),(0,0))',[0,1,2,3],[1,2,3,4,5]))

sumlambda=1
for i in range(1,n_leaf-1):
	sumlambda =sumlambda+c[i]*-log(t[i])
#print(sumlambda)


file_name_lambda_formula = 'lambda_formula.m'
output_file_path = os.path.join(current_directory, file_name_lambda_formula)
file = open(output_file_path, "w+")
for i,m in enumerate(markerdata):
	en = list(range(n_edge+1))
	ln = list(range(1,n_leaf+1))
	f = recursive_split(m,en,ln)
	fx = str(f[0]+f[3])
	file.write(f'l({i+1}) = {fx};')
	file.write('\r\n')

file_name_lambda_formula_trivial = 'lambda_formula_trivial.m'
output_file_path = os.path.join(current_directory, file_name_lambda_formula_trivial)
file = open(output_file_path, "w+")
for i,m in enumerate(trivialtopo):
	en = list(range(n_edge+1))
	ln = list(range(1,n_leaf+1))
	f = recursive_split(m,en,ln)
	fx = str(f[0]+f[3])
	file.write(f'la({i+1}) = {fx};')
	file.write('\r\n')



user_choice = sys.argv[3]
if user_choice == 'yes':
	file_name_variables = 'variables.m'
	output_file_path = os.path.join(current_directory, file_name_variables)
	file = open(output_file_path, "w+")
	for i in range(1,len(t)):
		file.write(f't{i}=x({i});')
		file.write('\r\n')
	for j,g in enumerate(gmloc):
		file.write(f'tl{g}=x({j+len(t)});')
		file.write('\r\n')
	for k in range(1,len(c)):
		file.write(f'c{k}=x({k+len(t)+len(gmloc)-1});')
		file.write('\r\n')
	file.close()


	file_name_lambda_formula = 'treescore.m'
	output_file_path = os.path.join(current_directory, file_name_lambda_formula)
	file = open(output_file_path, "w+")
	file.write('function [F] = treescore()\nclear all\n\nobjective = @(x) myObjective(x);\n')
	file.write(f'x0 = ones(1,{n_edge*2+len(gm)});\nlb = [zeros(1,{n_edge+len(gm)}),ones(1,{n_edge})*0.001];\nub = [ones(1,{n_edge+len(gm)}),ones(1,{n_edge})*1000];\n\n')
	file.write("options = optimoptions('fmincon', 'Display', 'off','MaxFunctionEvaluations',50000,'MaxIterations',10000);\n")
	file.write("[x_optimal,fval,exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);\n\n")
	file.write("F= -fval\n\nend\n\n")
	file.write("function result = myObjective(x)\n\tvariables\n\tlambda_formula\n\tlambda_formula_trivial\n\t")
	file.write(f"suml={sumlambda}-sum(la);\n\t")
	file.write('F = sum(log(l))-length(l)*log(suml);\n\tresult = -F;\nend\n')
	file.close()



else:
	file_name_variables = 'variables.m'
	output_file_path = os.path.join(current_directory, file_name_variables)
	file = open(output_file_path, "w+")
	for i in range(1,len(t)):
		file.write(f't{i}=x({i});')
		file.write('\r\n')
	for j,g in enumerate(gmloc):
		file.write(f'tl{g}=x({j+len(t)});')
		file.write('\r\n')
	for k in range(1,len(c)):
		file.write(f'c{k}=1;')
		file.write('\r\n')
	file.close()


	file_name_lambda_formula = 'treescore.m'
	output_file_path = os.path.join(current_directory, file_name_lambda_formula)
	file = open(output_file_path, "w+")
	file.write('function [F] = treescore()\nclear all\n\nobjective = @(x) myObjective(x);\n')
	file.write(f'x0 = ones(1,{n_edge+len(gm)});\nlb = zeros(1,{n_edge+len(gm)});\nub = ones(1,{n_edge+len(gm)});\n\n')
	file.write("options = optimoptions('fmincon', 'Display', 'off','MaxFunctionEvaluations',50000,'MaxIterations',10000);\n")
	file.write("[x_optimal,fval,exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);\n\n")
	file.write("F= -fval\n\nend\n\n")
	file.write("function result = myObjective(x)\n\tvariables\n\tlambda_formula\n\tlambda_formula_trivial\n\t")
	file.write(f"suml={sumlambda}-sum(la);\n\t")
	file.write('F = sum(log(l))-length(l)*log(suml);\n\tresult = -F;\nend\n')
	file.close()


