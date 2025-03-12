import subprocess
import matlab.engine
import io
import math
import re
import os


def find_matching_parentheses(s):
	'''
	input: a string, eg: '((a,b),(c,d),e,f)'
	return: a list of tuples, each tuple contains the index of matching parentheses
	'''
	stack = []
	result = []

	for i, char in enumerate(s):
		if char == '(':
			stack.append(i)
		elif char == ')':
			if stack:
				result.append([stack.pop(), i])
	
	#order the result by the first element of each tuple
	result.sort(key=lambda x: x[0])

	return result


def findsplit(s):
	'''
	input: a newick tree string, eg: '(a,b),(c,d),e,f'
	return: a list of strings, each string is a part of the tree
	'''
	result = []
	stack = []
	current_part = ""

	for char in s:
		if char == ',' and not stack:
			result.append(current_part.strip())
			current_part = ""
		else:
			current_part += char
			if char == '(':
				stack.append(char)
			elif char == ')':
				if stack:
					stack.pop()

	result.append(current_part.strip())
	return result


def find_minimum_parentheses(s, target):
	# 用一个栈来追踪括号的开闭
	stack = []
	min_paren = None
	start_idx = None
	
	# 遍历字符串
	for i, char in enumerate(s):
		if char == '(':
			# 遇到左括号，开始一个新的层次
			stack.append(i)
		elif char == ')':
			# 遇到右括号，结束一个层次
			start = stack.pop()
			# 如果这个括号是我们关注的括号，判断X是否在其中
			if min_paren is None and target in s[start + 1:i]:
				min_paren = (start, i)
				start_idx = start + 1
		
	return s[min_paren[0]:min_paren[1] + 1] if min_paren else None


def treescore(treetopo,markerfile,usec,useq):
	print('tree:',treetopo)
	if useq =='yes':
		subprocess.run(["python3", "newlambda_q.py", treetopo, markerfile, usec])
		eng = matlab.engine.start_matlab()
		stdout_capture = io.StringIO()
		LogL,t,c,q=eng.treescore(nargout =4, stdout=stdout_capture)
		print("LogL = ", LogL)
		eng.exit()
	else:
		subprocess.run(["python3", "newlambda_cpp.py", treetopo, markerfile, usec])
		eng = matlab.engine.start_matlab()
		stdout_capture = io.StringIO()
		LogL,t,c =eng.treescore(nargout =3, stdout=stdout_capture)
		print("LogL = ", LogL)
		eng.exit()

	return LogL

def nni(s,s_cluster):

	subtree = find_minimum_parentheses(s,s_cluster)
	two_clades = findsplit(s_cluster[1:-1])
	lst = findsplit(subtree[1:-1])
	third_clade = lst[0] if lst[0] != s_cluster else lst[1]
	 
	new_nni1 = '((' + two_clades[0]+','+third_clade+'),'+ two_clades[1]+')'
	new_nni2 = '((' + two_clades[1]+','+third_clade+'),'+ two_clades[0]+')'
	
	newtree1 = s.replace(subtree,new_nni1)
	newtree2 = s.replace(subtree,new_nni2)
	
	edge_index1 = newtree1.split(new_nni1)[0].count('(')+1
	edge_index2 = newtree2.split(new_nni2)[0].count('(')+1
	
	return newtree1,newtree2,edge_index1,edge_index2


def is_valid_path(path):
    return len(path) < 255 and os.path.exists(path)

user_input = input('Please input the tree(Newick string or filename): ')
if is_valid_path(user_input):
    with open(user_input, 'r') as file:
        newick_str = file.read().strip()
else:
    newick_str = user_input
newick_str = newick_str.replace(" ", "").replace(";", "")
newick_str = re.sub(r':\d+(\.\d+)?([eE][+-]?\d+)?', '', newick_str)
treetopo = re.sub(r'\)(\d+(\.\d+)?([eE][+-]?\d+)?)', ')', newick_str)


markerfile = input("Please input the marker filename: ")
usec = input("use parameter c (yes/no): ").lower()
useq = input("use parameter q (yes/no): ").lower()


F0 = treescore(treetopo,markerfile,usec,useq)
treecluster = [treetopo[start:end+1] for start, end in find_matching_parentheses(treetopo)[1:]]

NNIscore =[]
support =[]
for edge in range(len(treecluster)):
	new_tree_1, new_tree_2, edge_index1, edge_index2 = nni(treetopo, treecluster[edge])
	
	# Calculate the F value for both new trees
	F1 = treescore(new_tree_1,markerfile,usec,useq)
	F2 = treescore(new_tree_2,markerfile,usec,useq)

	NNIscore.append([F1,F2])
	result = 1 / (1 + math.exp(F1 - F0) + math.exp(F2 - F0))
	support.append(result)

print('support value:',support)

