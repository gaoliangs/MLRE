import subprocess
import matlab.engine
import io
import re
import os
from io import StringIO
from Bio import Phylo


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
			content_taxa = s[start + 1:i].replace("(","").replace(')',"").split(",")
			target_taxa = target.replace("(","").replace(')',"").split(",")
			if (min_paren is None) and (set(target_taxa) < set(content_taxa)):
				min_paren = (start, i)
				start_idx = start + 1
		
	return s[min_paren[0]:min_paren[1] + 1] if min_paren else None



def find_target_bracket(s, target_cluster):
	# 定义一个栈来追踪括号位置
	stack = []
	target_element = target_cluster.replace('(','').replace(')','').split(',')
	
	# 遍历字符串，检查每个括号的内容
	for i in range(len(s)):
		if s[i] == '(':
			stack.append(i)  # 记录左括号位置
		elif s[i] == ')':
			start = stack.pop()  # 获取最近的左括号位置
			# 获取当前括号内的子串
			content = s[start+1:i]
			content_elements = content.replace('(','').replace(')','').split(',')
			# 如果目标子串在括号内，返回该括号的范围
			if set(target_element) == set(content_elements):
				return s[0:start].count('(') 
	return None  # 如果没有找到，返回None

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

#print(nni('((((((newdataMJLFD,newdataK),BALRE),((CHLUN,TAUER),CUCCA)),newdataH),COLLI),newdataI)','((CHLUN,TAUER),CUCCA)'))



def spr_x(s,s_subtree):
	#prune the subtree
	cut_subtree = find_minimum_parentheses(s, s_subtree)
	lst = findsplit(cut_subtree[1:-1])
	sister_subtree = lst[0] if lst[0] != s_subtree else lst[1]
	#print(sister_subtree)
	pruned_tree = s.replace(cut_subtree,sister_subtree)

	#regraft the subtree
	sprtree = []
	bracket_in_pruned = find_matching_parentheses(pruned_tree)
	cluster_in_pruned = [pruned_tree[i[0]:i[1]+1] for i in bracket_in_pruned]
	taxa = pruned_tree.replace('(','').replace(')','').split(',') 
	subtree = taxa+cluster_in_pruned
	subtree.remove(sister_subtree)

	#remove nni tree
	if '(' in sister_subtree:
		sister_children = findsplit(sister_subtree[1:-1])
		subtree = [item for item in subtree if item not in sister_children]

	if pruned_tree != sister_subtree:
		subtree_parent = find_minimum_parentheses(pruned_tree, sister_subtree)
		subtree.remove(subtree_parent)
		lst = findsplit(subtree_parent[1:-1])
		parent_neighbor = lst[0] if lst[0] != sister_subtree else lst[1]
		subtree.remove(parent_neighbor)

	for i in subtree:
		if "(" not in i:
			pattern = rf'\b{i}\b'  
		else:
			pattern = re.escape(i)
		new_spr = re.sub(pattern, f"({i},{s_subtree})", pruned_tree)
		sprtree.append(new_spr)

	return sprtree

def spr(s):
	cluster = [s[start:end+1] for start, end in find_matching_parentheses(s)[1:]]
	taxa = s.replace('(','').replace(')','').split(',')
	
	subtree = cluster+taxa
	spr = []
	for i in subtree:
		new_spr = spr_x(s,i)
		spr.append(new_spr)

	sprtree = [item for sublist in spr for item in sublist]
	return sprtree


def edge_distance(s,target_edge):
	#find the target edge
	bracket = find_matching_parentheses(s)
	target = bracket[target_edge]

	#find the distance between the target edge and the other edges
	edge_distance = []
	for i in bracket[1:]:
		if (i[0] < target[0]) & (i[1] > target[1]):
			edge_distance.append(s[i[0]:target[0]].count('(')-s[i[0]:target[0]].count(')'))
		elif (i[0] > target[0]) & (i[1] < target[1]):
			edge_distance.append(s[target[0]:i[0]].count('(')-s[target[0]:i[0]].count(')'))
		elif i == target:
			edge_distance.append('max')
		elif (i[0] < target[0]) & (i[1] < target[0]):
			edge_distance.append(s[i[1]:target[0]].count('(')+s[i[1]:target[0]].count(')')-2*len(find_matching_parentheses(s[i[1]:target[0]])))
		elif (i[0] > target[1]) & (i[1] > target[1]):
			edge_distance.append(s[target[1]:i[0]].count('(')+s[target[1]:i[0]].count(')')-2*len(find_matching_parentheses(s[target[1]:i[0]])))
		
	return edge_distance


#search_tree = '(((((((((APAVI,(LEPDI,PICPU_MERNU)),BUCRH),((CATAU,HALLE_HALAL),TYTAL)),(CARCR,((FALPE,NESNO_MELUN),TAEGU_MANVI_GEOFO_CORBR_ACACH))),COLST),(((((BALRE,((((EGRGA,(NIPNI,PELCR)),PHACA),(FULGL,PYGAD_APTFO)),GAVST)),CHAVO),OPHHO),(CAPCA,CHAPE_CALAN)),(CHLUN,TAUER))),(EURHE,PHALE)),(COLLI,(CUCCA,PODCR_PHORU))),(MESUN,PTEGU))'
#constraint_tree = '((CAPCA,CHAPE_CALAN),(EURHE,PHALE),(MESUN,PTEGU),(((FULGL,PYGAD_APTFO),((EGRGA,NIPNI,PELCR),PHACA)),GAVST),((((NESNO_MELUN,TAEGU_MANVI_GEOFO_CORBR_ACACH),FALPE),CARCR),((CATAU,HALLE_HALAL),((BUCRH,PICPU_MERNU),APAVI,LEPDI),TYTAL),COLST),BALRE,CUCCA,CHLUN,COLLI,TAUER,OPHHO,CHAVO,PODCR_PHORU)'


user_input = input('Please input the initial tree(enter Newick tree or press enter to use default): ')
constraint = input('Please input the constraint tree(press enter to skip, buneman for default, or enter Newick tree): ')
markerfile = input('Please input the marker filename (CSV format):')
usec = input("use parameter c (yes/no): ").lower()
useq = input("use parameter q (yes/no): ").lower()


def remove_information(tree):
    for clade in tree.find_clades():
        clade.branch_length = None  # 去掉枝长
        clade.comment = None     # 去掉注释
        clade.confidence = None  # 去掉支持值
    return tree


if os.path.exists(constraint):
	with open(constraint, 'r') as file:
		newick_str = file.read().strip()
		tree = Phylo.read(StringIO(newick_str), "newick")
		topology_str = remove_information(tree).format("newick")
		str_noedge = re.sub(r":\d+(\.\d+)?", "", topology_str)
		constraint_tree = str_noedge.strip().replace(";", "")

else:
	if constraint == "":
		with open(markerfile, 'r') as csv_file:
			csv_content = csv_file.readlines()
			constraint_tree = '('+csv_content[0].strip()+')'
	elif constraint.lower()  == "buneman":
		threshold = input('Please input the constraint threshold(default:0.95): ')
		if threshold == '':
			threshold = str(0.95)
		result = subprocess.run(["python3", "bunemantree.py", markerfile, threshold], stdout=subprocess.PIPE, text=True)
		constraint_tree = result.stdout.strip()
	else:
		tree = Phylo.read(StringIO(constraint), "newick")
		topology_str = remove_information(tree).format("newick")
		str_noedge = re.sub(r":\d+(\.\d+)?", "", topology_str)
		constraint_tree = str_noedge.strip().replace(";", "")
print(constraint_tree)

if os.path.exists(user_input):
	with open(user_input, 'r') as file:
		newick_str = file.read().strip()
		tree = Phylo.read(StringIO(newick_str), "newick")
		topology_str = remove_information(tree).format("newick")
		str_noedge = re.sub(r":\d+(\.\d+)?", "", topology_str)
		search_tree = str_noedge.strip().replace(";", "")
else:
	if user_input == "":
		result = subprocess.run(["python3", "initialtree.py", markerfile,constraint_tree], stdout=subprocess.PIPE, text=True)
		search_tree = result.stdout.strip()
	else:
		tree = Phylo.read(StringIO(user_input), "newick")
		topology_str = remove_information(tree).format("newick")
		str_noedge = re.sub(r":\d+(\.\d+)?", "", topology_str)
		search_tree = str_noedge.strip().replace(";", "")
print(search_tree)


def check_inclusion(search_tree, con_tree):
	clusters_search_tree = [search_tree[start:end+1] for start, end in find_matching_parentheses(search_tree)]
	clusters_constraint_tree = [con_tree[start:end+1] for start, end in find_matching_parentheses(con_tree)]
	delete = []
	for cluster in clusters_constraint_tree:
		cluster_constraint = cluster.replace(')','').replace('(','').split(',')
		for taxon in clusters_search_tree:
			clusters_search = taxon.replace(')','').replace('(','').split(',')
			if set(cluster_constraint) == set(clusters_search):
				#print(cluster)
				break
		else:
			delete.append(cluster)
	constrainttree = con_tree
	for cluster in delete:
		constrainttree = constrainttree.replace(cluster, cluster[1:-1])
		
	return constrainttree

print('The search tree is:', search_tree)
constraint_tree = check_inclusion(search_tree,constraint_tree)
print('The constraint tree is:', constraint_tree)
constraint_cluster = [constraint_tree[start:end+1] for start, end in find_matching_parentheses(constraint_tree)[1:]]


print("Beginning the search...")


#################NNI search#################
# Initialize the best tree and its F value

best_tree = search_tree
best_F_value = treescore(search_tree,markerfile,usec,useq)


current_tree = ''
while current_tree != best_tree:
	current_tree = best_tree

	#################NNI search#################
	print('start NNI')
	# Initialize the edge index to 0
	edge_num = [0]*(current_tree.count(',')-1)
	while edge_num != ['max']*(current_tree.count(',')-1):
		current_tree = best_tree

		constraint_edge = [find_target_bracket(current_tree,cc) for cc in constraint_cluster]
		for idx in constraint_edge:
			edge_num[idx-1] = 'max'
		search_cluster = [current_tree[start:end+1] for start, end in find_matching_parentheses(current_tree)[1:]]
		numbers = [(i, value) for i, value in enumerate(edge_num) if isinstance(value, (int, float))]
		sorted_numbers = sorted(numbers, key=lambda x: x[1])
		sorted_indices = [i for i, _ in sorted_numbers]
		#print(edge_num)

		for edge in sorted_indices:
			new_tree_1, new_tree_2, edge_index1, edge_index2 = nni(current_tree, search_cluster[edge])
			
			# Calculate the F value for both new trees
			F_value_1 = treescore(new_tree_1,markerfile,usec,useq)
			F_value_2 = treescore(new_tree_2,markerfile,usec,useq)
			
			# If a better tree is found, update the best tree
			if F_value_1 > max(best_F_value,F_value_2):
				best_tree = new_tree_1
				best_F_value = F_value_1
				edge_num = edge_distance(best_tree,edge_index1)
				print('update',edge_num)
				break  # Exit the loop and re-evaluate the tree
			
			if F_value_2 > max(best_F_value,F_value_1):
				best_tree = new_tree_2
				best_F_value = F_value_2
				edge_num = edge_distance(best_tree,edge_index2)
				print('update',edge_num)
				break  # Exit the loop and re-evaluate the tree

			if best_F_value > max(F_value_1,F_value_2):
				edge_num[edge] = 'max'
				print(edge_num)

	print('The best NNI tree is:', best_tree)
	print('The best NNI score is:', best_F_value)


	#################SPR search#################
	print('start SPR')
	current_tree = best_tree

	result =[]
	def recursive_split(s):
		parts = findsplit(s[1:-1])
		if len(parts)>2:
			result.append(s)
	    
		for part in parts:
			if '(' in part and ')' in part:  # part contains nested brackets
				layer = recursive_split(part)

		return result

	unresolved_node = recursive_split(constraint_tree)
	#print('unresolved_node',unresolved_node)


	def find_corresponding_elements(listA, listB):
		string = []
		node = []
		# 对listA中的每个元素，查找其对应的listB元素
		for a in listA:
			taxon_a = a.replace(')','').replace('(','').split(',')
			for b in listB:
				taxon_b = b.replace(')','').replace('(','').split(',')
				if set(taxon_a) == set(taxon_b):
					string.append(a)
					node.append(b)

		return string,node


	def replace_with_split(s, replacement_list):
		replacements = {replacement: f"split{i+1}" for i, replacement in enumerate(replacement_list)}

		new_s = s
		for old, new in replacements.items():
			new_s = new_s.replace(old, new)

		return new_s, replacements  # 返回替换后的字符串和替换字典


	def restore_replacements(result_list, replacements):
		restored_list = []
		for item in result_list:
			restored_item = item
			for split, original in replacements.items():
				restored_item = restored_item.replace(original, split)
			restored_list.append(restored_item)
		return restored_list


	tree_cluster = [current_tree[start:end+1] for start, end in find_matching_parentheses(current_tree)]
	string, node = find_corresponding_elements(tree_cluster,unresolved_node)
	#print('aa',string,node)

	result = []
	for i in range(len(string)):
		split_list = findsplit(node[i][1:-1])
		st, n = find_corresponding_elements(tree_cluster,split_list)
		#print('mm',st,n)
		new_s, replacements = replace_with_split(string[i],st)
		#print('new_s',new_s,replacements)
		spr_subtree = spr(new_s)
		#print(spr_subtree)
		spr_tree = restore_replacements(spr_subtree,replacements)
		if len(spr_tree) > 0:
			sprtree = [current_tree.replace(string[i],spr) for spr in spr_tree]
			result.append(sprtree)

	flattened_result = [item for sublist in result for item in sublist]
	#print(flattened_result)

	for tree in flattened_result:
		F = treescore(tree,markerfile,usec,useq)
		if F > best_F_value:
			best_tree = tree
			best_F_value = F
			break

	print('SPR tree:',best_tree)
	print('SPR score:',best_F_value)
