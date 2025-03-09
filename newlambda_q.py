from functools import reduce
import itertools
import re
import sys




topo = sys.argv[1]
taxaorder = topo.replace("(", "").replace(")", "").split(',')

allmarker = list(itertools.product(["0", "1","01",'?'], repeat=len(taxaorder)))
realmarker = [list(cha) for cha in allmarker if (1 < cha.count("1")+cha.count("01")) & (cha.count("1")+cha.count('?')<len(cha) )]
def sort_key(s):
    return ('?' in s, s)
realmarker = sorted(realmarker, key=sort_key)



file_name = sys.argv[2]
with open(file_name, 'r') as csv_file:
    csv_content = csv_file.readlines()
taxaname = csv_content[0].strip().split(',')
taxadata = [line.strip().replace('B','01').split(',') for line in csv_content[1:]]



#grouptaxa,singletaxa的taxa名
singletaxa = []
grouptaxa = []
for col_index in range(len(taxaname)):
   # 判断该列是否存在 'B'
   if any(row[col_index] == '01' for row in taxadata):
      grouptaxa.append(taxaname[col_index])
   else:
      singletaxa.append(taxaname[col_index])


smloc = [taxaorder.index(i)+1 for i in singletaxa]
gm = [taxaname.index(i)+1 for i in grouptaxa]
gmloc = [taxaorder.index(i)+1 for i in grouptaxa]
#print(smloc)

realm = [
    sublist for sublist in realmarker
    if all(sublist[i] != '01' for i in [taxaname.index(i) for i in singletaxa])
]


markerdata = []
for mi in realm:
	orginial = topo
	for i,name in enumerate(taxaname):
		pattern = re.compile(r'\b' + re.escape(name) + r'\b')
		orginial = pattern.sub(mi[i], orginial, count=1)
	markerdata.append(orginial)


kvector = [taxadata.count(item) for item in realm]


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
variables_q = [f"q{i}" for i in range(1, n_leaf + 1)]

# 叶节点的矩阵
matrix_leaf_group = [
    ["0", "0", "0", "0"],
    ["0", "1", "0", "0"],  
    ["0", "0", "1", "0"],
    ["0", "(1-tl0)*0.5", "(1-tl0)*0.5", "tl0"]
]


# 内部节点的矩阵
matrix_internal = [
    ["1", "(-log(t0)-(1-t0))*exp(c0)*0.5", "(-log(t0)-(1-t0))*exp(c0)*0.5","exp(c0)*(1-t0)"],
    ["0", "1", "0", "0"],
    ["0", "0", "1", "0"],
    ["0", "(1-t0)*0.5", "(1-t0)*0.5", "t0"]
]


matrix_leaf = [
    ["0", "0", "0", "0"],
    ["0", "1", "0", "0"],
    ["0", "0", "1", "0"],
    ["0", "0.5", "0.5", "1"]
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
	v_not_invented = [i[0] for i in vector]
	v_lost= [i[1] for i in vector]
	v_fixed= [i[2] for i in vector]
	v_polymorphic = [i[3] for i in vector]
	

	v_l = ["1"] if all(element[0] == "1" for element in v_lost) else ["0"]
	v_f = ["1"] if all(element[0] == "1" for element in v_fixed) else ["0"]
	v_p = reduce(lambda x, y: ["0"] if x[0] == "0" or y[0] == "0" else [f"{x[0]}*{y[0]}"],v_polymorphic)

	if all(my_element[0] == "0" for my_element in v_lost):
		v_n= ["0"]
	elif all(my_element[0] == "1" for my_element in v_lost):
		#v_n = reduce(lambda x, y: [f"({x[0]}+{y[0]})"], v_not_invented)
		v_n = reduce(lambda x, y: ["0"] if x[0] == "0" and y[0] == "0" else [f"({x[0]}+{y[0]})"], v_not_invented)
	else:
		introduce_index = [index for index, value in enumerate(v_lost) if value == ["0"]]
		v_n_i = [v_not_invented[i] for i in introduce_index]

		non_zero_elements = [element[0] for element in v_n_i if element[0] != "0"]
		v_n = [" + ".join(non_zero_elements)] if non_zero_elements else ["0"]

	return [v_n,v_l,v_f,v_p]



def matrix_vector(matrix_A, vector_B):
    rows_A = len(matrix_A)
    cols_A = len(matrix_A[0])
    rows_B = len(vector_B)

    expressions = []

    for i in range(rows_A):
        terms = [
            vector_B[j][0] if matrix_A[i][j] == "1" else
            matrix_A[i][j] if vector_B[j][0] == "1" else
            f"{matrix_A[i][j]}*{vector_B[j][0]}"
            for j in range(cols_A)
            if matrix_A[i][j] != "0" and vector_B[j][0] != "0"
        ]
        if not terms:  # 如果没有有效项，则结果为 "0"
            expressions.append(["0"])
        #elif all(matrix_A[i][j] == "1" and matrix_B[j][0] == "1" for j in range(cols_A) if matrix_A[i][j] != "0"):
        #    expressions.append(["1"])  # 所有有效项都是 "1"
        elif len(terms) == 1:
            expressions.append(terms)
        else:
            expre = "+".join(terms)
            expressions.append([f"({expre})"])


    return expressions





def recursive_split(s,edge_num,leaf_index):
	if '(' in s:
		parts = findsplit(s,edge_num,leaf_index)
		result = []
		for part in parts:
			if '(' not in part[0]:
				v = recursive_split(part[0],part[1],part[2])
				result.append(v)
			else:
				#print('p',part)
				v = recursive_split(part[0],part[1],part[2])
				e = part[1][0]
				current_matrix = [[elem.replace('t0', f"t{e}").replace("c0", f"c{e}") for elem in row]  for row in matrix_internal]
				#print('pv',current_matrix,v)
				result.append(matrix_vector(current_matrix,v))

		vectors = combine(result)
		return vectors
	else:
		if leaf_index[0] in smloc:
			current_matrix = [[elem.replace("tl0", f"tl{leaf_index[0]}") for elem in row]for row in matrix_leaf]
		else:
			current_matrix = [[elem.replace("tl0", f"tl{leaf_index[0]}") for elem in row]for row in matrix_leaf_group]
		if s == '0':
			return matrix_vector(current_matrix,[["1"],["1"],["0"],["0"]])
		if s == '1':
			return matrix_vector(current_matrix,[["0"],["0"],["1"],["0"]])
		if s == '01':
			return matrix_vector(current_matrix,[["0"],["0"],["0"],["1"]])
		if s == '?':
			return [["1"],["1"],["1"],["1"]]

#print(recursive_split('(((0,0),0),(0,0))',[0,1,2,3],[1,2,3,4,5]))


def question_mark(m,f):
	for n,char in enumerate(m):
		if char == '?':
			f = [[item[0] + '*' +variables_q[n]] for item in f]
		else:
			f = [[item[0] + '*(1-' +variables_q[n]+')'] for item in f]
	return f


sumlambda="1"
for i in range(1,n_edge+1):
	sumlambda = sumlambda+f"-exp(x({i+n_edge+len(gmloc)}))*log(x({i}))"
#print(sumlambda)


file_name_lambda_formula = 'lambda_formula.cpp'
file = open(file_name_lambda_formula, "w+")

file.write("#include <mex.h>\n#include <vector>\n\n")
file.write("void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[])\n{\n\n")
file.write("\tdouble "+", ".join(variables_t)+';\n')
if len(gmloc)>0:
	file.write("\tdouble "+','.join([f"tl{i}" for i in gmloc])+';\n')
file.write("\tdouble "+", ".join(variables_c)+';\n')
file.write("\tdouble "+", ".join(variables_q)+';\n')
file.write(f"\tconst mxArray* ma = prhs[0];\n\tif(mxGetN(ma) * mxGetM(ma) != {2*n_edge+len(gm)+n_leaf})\n")
file.write('\t\tmexErrMsgTxt("error!\\n");\n\tdouble* x = mxGetDoubles(ma);\n\t')
t_list = ";".join([f"t{i+1}=x[{i}]" for i in range(n_edge)]) + ";\n\t"
file.write(t_list)
if len(gmloc)>0:
	tl_list = ";".join([f"tl{j}=x[{i+n_edge}]" for i,j in enumerate(gmloc)]) + ";\n\t"
	file.write(tl_list)
q_list = ";".join([f"q{i+1}=x[{i+n_edge+len(gmloc)}]" for i in range(n_leaf)]) + ";\n\t"
file.write(q_list)	
c_list = ";".join([f"c{i+1}=x[{i+n_edge+len(gmloc)+n_leaf}]" for i in range(n_edge)]) + ";\n\t"
file.write(c_list)
file.write(f"\n\n\tstd::vector<double> l({len(realm)});\n")



for i,m in enumerate(markerdata):
	en = list(range(n_edge+1))
	ln = list(range(1,n_leaf+1))
	f = recursive_split(m,en,ln)
	ma = m.replace('(','').replace(')','').split(',')
	f = question_mark(ma,f)
	fx = '+'.join([f[0][0],f[3][0]])
	file.write(f'l[{i}] = {fx};')
	file.write('\r\n')
file.write(f"plhs[0] = mxCreateDoubleMatrix({len(realm)}, 1, mxREAL);")
file.write(f"double* l0 = mxGetDoubles(plhs[0]);\nfor(int i =0; i < {len(realm)}; i++)\n\tl0[i] = l[i];\n")
file.write("}")
file.close()





user_choice = sys.argv[3]
#user_choice = 'no'
if user_choice == 'yes':

	file_name_script = 'treescore.m'
	file = open(file_name_script, "w+")
	file.write('function [F] = treescore()\nclear all\n\nmex -silent lambda_formula.cpp -R2018a\n objective = @(x) myObjective(x);\n')
	file.write(f'x0 = ones(1,{n_edge*2+len(gm)+n_leaf});\nlb = [zeros(1,{n_edge+len(gm)+n_leaf}),ones(1,{n_edge})*-10];\nub = [ones(1,{n_edge+len(gm)+n_leaf}),ones(1,{n_edge})*10];\n\n')
	file.write(f"options = optimoptions('fmincon', 'Display', 'off','MaxFunctionEvaluations',{(n_edge*2+len(gm)+n_leaf)*5000},'MaxIterations',5000);\n")
	file.write("[x_optimal,fval,exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);\n\n")
	file.write(f"disp('Optimal solution t:');\ndisp(-log(x_optimal(1:{n_edge+len(gm)})));\ndisp('Optimal solution q:');\ndisp(x_optimal({1+n_edge+len(gm)}:{n_edge+n_leaf+len(gm)}));\ndisp('Optimal solution c:');\ndisp(exp(x_optimal({n_edge+n_leaf+len(gm)+1}:{n_edge*2+len(gm)+n_leaf})));\ndisp('Optimal log likelihood:');\ndisp(-fval);\n")
	file.write("F= -fval;\n\nend\n\n")
	file.write("function result = myObjective(x)\n\tl=lambda_formula(x);\n\t")
	file.write('k='+str(kvector)+';\n\t')
	file.write("F = k*log(l)-sum(k)*log(sum(l));\n\tresult = -F;\nend\n")
	file.close()

else:
	file_name_script = 'treescore.m'
	file = open(file_name_script, "w+")
	file.write('function [F] = treescore()\nclear all\n\nmex -silent lambda_formula.cpp -R2018a\n objective = @(x) myObjective(x);\n')
	file.write(f'x0 = ones(1,{2*n_edge+len(gm)+n_leaf});\nlb = zeros(1,{2*n_edge+len(gm)+n_leaf});\nub = ones(1,{2*n_edge+len(gm)+n_leaf});\n')
	file.write(f'Aeq = zeros({n_edge}, length(x0));\nfor i = 1:{n_edge}\n\tAeq(i, i + {n_edge+len(gmloc)+n_leaf}) = 1;\nend\nBeq = zeros({n_edge}, 1);\n\n')
	file.write(f"options = optimoptions('fmincon', 'Display', 'off','MaxFunctionEvaluations',{(n_edge*2+len(gm)+n_leaf)*5000},'MaxIterations',5000);\n")
	file.write("[x_optimal,fval,exitflag] = fmincon(objective, x0, [], [], Aeq, Beq, lb, ub, [], options);\n\n")
	file.write(f"disp('Optimal solution t:');\ndisp(-log(x_optimal(1:{n_edge+len(gm)})));\ndisp('Optimal solution q:');\ndisp(x_optimal({n_edge+len(gm)+1}:{n_edge+len(gm)+n_leaf}));\ndisp('Optimal log likelihood:');\ndisp(-fval);\n")
	file.write("F= -fval;\n\nend\n\n")
	file.write("function result = myObjective(x)\n\tl=lambda_formula(x);\n\t")
	file.write('k='+str(kvector)+';\n\t')
	file.write("F = k*log(l)-sum(k)*log(sum(l));\n\tresult = -F;\nend\n")
	file.close()



