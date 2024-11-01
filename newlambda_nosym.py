from sympy import symbols, Matrix, log, exp, Sum, Mul
from functools import reduce
import re
import sys
import os

current_directory = os.getcwd()

topo = sys.argv[1]
#topo = '(((((((((FULGL,PYGAD_APTFO),(((EGRGA,PELCR),NIPNI),PHACA)),GAVST),((CAPCA,CHAPE_CALAN),((COLST,(((CATAU,HALLE_HALAL),TYTAL),(LEPDI,(APAVI,(BUCRH,PICPU_MERNU))))),(((TAEGU_MANVI_GEOFO_CORBR_ACACH,NESNO_MELUN),FALPE),CARCR)))),((OPHHO,(EURHE,PHALE)),CHAVO)),(((CHLUN,TAUER),CUCCA),BALRE)),PODCR_PHORU),COLLI),(PTEGU,MESUN))'
#topo = '((COLST,(newdataHH_CATAU,(newdataPM_BUCRH_APAVI_LEPDI,TYTAL))),newdataCAGTM_newdataMN_FALPE_CARCR)'
taxaorder = topo.replace("(", "").replace(")", "").split(',')


data = []
file_name = sys.argv[2]
#file_name = 'landbird.csv'
#file_name = 'newseq1.csv'
input_file_path = os.path.join(current_directory, file_name)


with open(input_file_path, 'r') as mark:
   line = mark.readlines()
data.append(line[0][:-1].split(","))
for i in range(1,len(line)):
   line1 = line[i][:-1].replace('B','01')
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
variables_t = [f"t{i}" for i in range(1, n_edge + 1)]
variables_c = [f"c{i}" for i in range(1, n_edge + 1)]
variables_tl = [f"tl{i}" for i in range(1, n_leaf + 1)]
#print(variables_t)



# 叶节点的矩阵
matrix_leaf_group = [
    ["0", "0", "0", "0"],
    ["0", "1", "0", "0"],  
    ["0", "0", "1", "0"],
    ["0", "(1 - tl0) / 2", "(1 - tl0) / 2", "tl0"]
]


# 内部节点的矩阵
matrix_internal = [
    ["1", "(-log(t0)-(1-t0))*exp(c0)/2", "(-log(t0)-(1-t0))*exp(c0)/2","exp(c0)*(1-t0)"],
    ["0", "1", "0", "0"],
    ["0", "0", "1", "0"],
    ["0", "(1-t0)/2", "(1-t0)/2", "t0"]
]


matrix_leaf = [
    ["0", "0", "0", "0"],
    ["0", "1", "0", "0"],
    ["0", "0", "1", "0"],
    ["0", "1/2", "1/2", "1"]
]



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
	#print(vector)
	v_not_invented = [i[0] for i in vector]
	v_lost= [i[1] for i in vector]
	v_fixed= [i[2] for i in vector]
	v_polymorphic = [i[3] for i in vector]
	#print('ser;',v_not_invented,v_lost,v_fixed,v_polymorphic)

	v_l = ["1"] if all(element[0] == "1" for element in v_lost) else ["0"]
	v_f = ["1"] if all(element[0] == "1" for element in v_fixed) else ["0"]
	v_p = reduce(lambda x, y: ["0"] if x[0] == "0" or y[0] == "0" else [f"({x[0]} * {y[0]})"],v_polymorphic)
	#print('sdf',v_p)

	if all(my_element[0] == "0" for my_element in v_lost):
		v_n= ["0"]
	elif all(my_element[0] == "1" for my_element in v_lost):
		v_n = reduce(lambda x, y: [f"({x[0]} + {y[0]})"], v_not_invented)
	else:
		introduce_index = [index for index, value in enumerate(v_lost) if value == ["0"]]
		v_n_i = [v_not_invented[i] for i in introduce_index]

		non_zero_elements = [element[0] for element in v_n_i if element[0] != "0"]
		v_n = [" + ".join(non_zero_elements)] if non_zero_elements else ["0"]

	#print('sdfsdf',v_n)
	return [v_n,v_l,v_f,v_p]


#print(combine([v1,v2]))

def generate_matrix_multiplication_expression(matrix_A, matrix_B):
    rows_A = len(matrix_A)
    cols_A = len(matrix_A[0])
    rows_B = len(matrix_B)
    cols_B = len(matrix_B[0])
    #print('iiii',matrix_A,matrix_B)

    expressions = []

    # 判断矩阵 B 是否为 4x1 或 4x4
    if cols_B == 1: 
        #print('sdaaaa') # 如果 B 是 4x1 矩阵，结果为 4x1
        for i in range(rows_A):
            terms = [
                matrix_B[j][0] if matrix_A[i][j] == "1" else
                matrix_A[i][j] if matrix_B[j][0] == "1" else
                f"{matrix_A[i][j]} * {matrix_B[j][0]}"
                for j in range(cols_A)
                if matrix_A[i][j] != "0" and matrix_B[j][0] != "0"
            ]
            if not terms:  # 如果没有有效项，则结果为 "0"
                expressions.append(["0"])
            elif all(matrix_A[i][j] == "1" and matrix_B[j][0] == "1" for j in range(cols_A) if matrix_A[i][j] != "0"):
                expressions.append(["1"])  # 所有有效项都是 "1"
            else:
                expre = "+".join(terms)
                #rint('sdf',['('+expre])
                expressions.append(['('+expre+')'])
        #print('werwytrueyu',expressions)

    else:  # 如果 B 是 4x4 矩阵，结果为 4x4
        for i in range(rows_A):
            row_expr = []
            for j in range(cols_B):
                terms = [
                    f"{matrix_A[i][k]} * {matrix_B[k][j]}"
                    for k in range(cols_A)
                    if matrix_A[i][k] != "0" and matrix_B[k][j] != "0"
                ]
                if not terms:  # 如果没有有效项，则结果为 "0"
                    row_expr.append("0")
                elif all(matrix_A[i][k] == "1" and matrix_B[k][j] == "1" for k in range(cols_A) if matrix_A[i][k] != "0"):
                    row_expr.append("1")  # 所有有效项都是 "1"
                else:
                    row_expr.append(" + ".join(terms))

            expressions.append(row_expr)
    #print('ssssss nn',expressions)

    return expressions





def recursive_split(s,edge_num,leaf_index):
	if '(' in s:
		#print('ekk')
		parts = findsplit(s,edge_num,leaf_index)
		result = []
		for part in parts:
			if '(' not in part[0]:
				v = recursive_split(part[0],part[1],part[2])
				#print('1111',v)
				result.append(v)
				#print('werwerw',result)
			else:
				#print('p',part)
				v = recursive_split(part[0],part[1],part[2])
				#print('s',v)
				e = part[1][0]
				current_matrix = [[elem.replace('t0', f"t{e}").replace("c0", f"c{e}") for elem in row]  for row in matrix_internal]
				#print('pv',current_matrix,v)
				result.append(generate_matrix_multiplication_expression(current_matrix,v))

		#print('rrr',result)

		vectors = combine(result)
		#print('qweqweqwewqwe',vectors)
		return vectors
	else:
		#print('ss')
		if leaf_index[0] in smloc:
			current_matrix = [[elem.replace("tl0", f"tl{leaf_index[0]}") for elem in row]for row in matrix_leaf]
		else:
			current_matrix = [[elem.replace("tl0", f"tl{leaf_index[0]}") for elem in row]for row in matrix_leaf_group]
		if s == '0':
			#print('1')
			#print('sdf',generate_matrix_multiplication_expression(current_matrix,[["1"],["1"],["0"],["0"]]))
			return generate_matrix_multiplication_expression(current_matrix,[["1"],["1"],["0"],["0"]])
		if s == '1':
			#print('2')
			#print('sdf',generate_matrix_multiplication_expression(current_matrix,[["0"],["0"],["1"],["0"]]))
			return generate_matrix_multiplication_expression(current_matrix,[["0"],["0"],["1"],["0"]])
		if s == '01':
			#print('3')
			return generate_matrix_multiplication_expression(current_matrix,[["0"],["0"],["0"],["1"]])
		if s == '?':
			return [["1"],["1"],["1"],["1"]]

#print(recursive_split('(((0,0),0),(0,0))',[0,1,2,3],[1,2,3,4,5]))
#print(recursive_split('(1,1)', [2], [1, 2]))




sumlambda="1"
for i in range(1,n_leaf-1):
	sumlambda = sumlambda+f"-exp(c{i})*log(t{i})"
#print(sumlambda)




file_name_lambda_formula = 'lambda_formula.m'
output_file_path = os.path.join(current_directory, file_name_lambda_formula)
file = open(output_file_path, "w+")
for i,m in enumerate(markerdata):
	en = list(range(n_edge+1))
	ln = list(range(1,n_leaf+1))
	f = recursive_split(m,en,ln)
	fx = '+'.join([f[0][0],f[3][0]])
	file.write(f'l({i+1}) = {fx};')
	file.write('\r\n')



file_name_lambda_formula_trivial = 'lambda_formula_trivial.m'
output_file_path = os.path.join(current_directory, file_name_lambda_formula_trivial)
file = open(output_file_path, "w+")
for i,m in enumerate(trivialtopo):
	en = list(range(n_edge+1))
	ln = list(range(1,n_leaf+1))
	f = recursive_split(m,en,ln)
	fx = '+'.join([f[0][0],f[3][0]])
	file.write(f'la({i+1}) = {fx};')
	file.write('\r\n')


user_choice = sys.argv[3]
#user_choice = 'yes'
if user_choice == 'yes':
	file_name_variables = 'variables.m'
	output_file_path = os.path.join(current_directory, file_name_variables)
	file = open(output_file_path, "w+")
	for i in range(len(variables_t)):
		file.write(f't{i+1}=x({i+1});')
		file.write('\r\n')
	for j,g in enumerate(gmloc):
		file.write(f'tl{g}=x({j+len(variables_t)+1});')
		file.write('\r\n')
	for k in range(1,len(variables_c)+1):
		file.write(f'c{k}=x({k+len(variables_t)+1+len(gmloc)-1});')
		file.write('\r\n')
	file.close()


	file_name_lambda_formula = 'treescore.m'
	output_file_path = os.path.join(current_directory, file_name_lambda_formula)
	file = open(output_file_path, "w+")
	file.write('function [F] = treescore()\nclear all\n\nobjective = @(x) myObjective(x);\n')
	file.write(f'x0 = ones(1,{n_edge*2+len(gm)});\nlb = [zeros(1,{n_edge+len(gm)}),ones(1,{n_edge})*-10];\nub = [ones(1,{n_edge+len(gm)}),ones(1,{n_edge})*10];\n\n')
	file.write("options = optimoptions('fmincon', 'Display', 'off','MaxFunctionEvaluations',50000,'MaxIterations',10000);\n")
	file.write("[x_optimal,fval,exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);\n\n")
	file.write(f"disp('Optimal solution t:');\ndisp(-log(x_optimal(1:{n_edge+len(gm)})));\ndisp('Optimal solution c:');\ndisp(exp(x_optimal({1+n_edge+len(gm)}:{n_edge+n_edge+len(gm)})));\ndisp('Optimal objective function value:');\ndisp(exp(vpa(-fval)));\n")
	file.write("F= -fval;\n\nend\n\n")
	file.write("function result = myObjective(x)\n\tvariables\n\tlambda_formula\n\tlambda_formula_trivial\n\t")
	file.write(f"suml={sumlambda}-sum(la);\n\t")
	file.write('F = sum(log(l))-length(l)*log(suml);\n\tresult = -F;\nend\n')
	file.close()



else:
	file_name_variables = 'variables.m'
	output_file_path = os.path.join(current_directory, file_name_variables)
	file = open(output_file_path, "w+")
	for i in range(len(variables_t)):
		file.write(f't{i+1}=x({i+1});')
		file.write('\r\n')
	for j,g in enumerate(gmloc):
		file.write(f'tl{g}=x({j+len(variables_t)+1});')
		file.write('\r\n')
	for k in range(1,len(variables_c)+1):
		file.write(f'c{k}=0;')
		file.write('\r\n')
	file.close()


	file_name_lambda_formula = 'treescore.m'
	output_file_path = os.path.join(current_directory, file_name_lambda_formula)
	file = open(output_file_path, "w+")
	file.write('function [F] = treescore()\nclear all\n\nobjective = @(x) myObjective(x);\n')
	file.write(f'x0 = ones(1,{n_edge+len(gm)});\nlb = zeros(1,{n_edge+len(gm)});\nub = ones(1,{n_edge+len(gm)});\n\n')
	file.write("options = optimoptions('fmincon', 'Display', 'off','MaxFunctionEvaluations',50000,'MaxIterations',10000);\n")
	file.write("[x_optimal,fval,exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);\n\n")
	file.write(f"disp('Optimal solution:');\ndisp(-log(x_optimal(1:{n_edge+len(gm)})));\ndisp('Optimal objective function value:');\ndisp(exp(vpa(-fval)));\n")
	file.write("F= -fval;\n\nend\n\n")
	file.write("function result = myObjective(x)\n\tvariables\n\tlambda_formula\n\tlambda_formula_trivial\n\t")
	file.write(f"suml={sumlambda}-sum(la);\n\t")
	file.write('F = sum(log(l))-length(l)*log(suml);\n\tresult = -F;\nend\n')
	file.close()



