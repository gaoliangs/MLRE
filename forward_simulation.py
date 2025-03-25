import random
import re
import math
import csv
import os
from io import StringIO
from Bio import Phylo

num_samples = input("Please input the number of samples: ")
num_samples = int(num_samples)

num_marker = input("Please input the number of markers: ")
num_marker = int(num_marker)

input_markerfile = input("Please input the marker filename: ")

user_input = input('Please input the tree(Newick string or filename): ')

def remove_information(tree):
    for clade in tree.find_clades():
        clade.branch_length = None 
        clade.comment = None
        clade.confidence = None 
    return tree

if os.path.exists(user_input):
	with open(user_input, 'r') as file:
		newick_str = file.read().strip()
		tree = Phylo.read(StringIO(newick_str), "newick")
else:
	tree = Phylo.read(StringIO(user_input), "newick")


t = [clade.branch_length for clade in tree.get_nonterminals()]
print("internal_edge: ",t)
leaf_length = [clade.branch_length for clade in tree.get_terminals()]
print("leaf_edge: ",leaf_length)

c = [clade.comment for clade in tree.get_nonterminals()]
c = [1.0 if x is None else x for x in c]
print("insertion_rate: ",c)


topology_str = remove_information(tree).format("newick")
str_noedge = re.sub(r":\d+(\.\d+)?", "", topology_str)
treetopo = str_noedge.strip().replace(";", "")
taxa_list = treetopo.replace('(','').replace(')','').split(',')
print("taxa: ",taxa_list)

mapping = dict(zip(taxa_list,leaf_length ))

def find_matching_parentheses(s):
	stack = []
	result = []

	for i, char in enumerate(s):
		if char == '(':
			stack.append(i)
		elif char == ')':
			if stack:
				result.append([stack.pop(), i])

	return result


def findsplit(s):
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


def generate_marker(clu,tvalue):
	left = findsplit(clu[1:-1])[0]
	right = findsplit(clu[1:-1])[1]

	if '(' in left:
		left_p = [(1-math.exp(-tvalue[1]))/2,math.exp(-tvalue[1]),(1-math.exp(-tvalue[1]))/2]
	else:
		if mapping[left] is None:
			left_p = [1/2, 0, 1/2]
		else:
			tlp = math.exp(-mapping[left])
			left_p = [(1-tlp)/2,tlp,(1-tlp)/2]
	if '(' in right:
		right_p = [(1-math.exp(-tvalue[1+left.count('(')]))/2,math.exp(-tvalue[1+left.count('(')]),(1-math.exp(-tvalue[1+left.count('(')]))/2]
	else:
		if mapping[right] is None:
			right_p = [1/2, 0, 1/2]
		else:
			tlp = math.exp(-mapping[right])
			right_p = [(1-tlp)/2,tlp,(1-tlp)/2]

	left_node = random.choices(range(len(left_p)), weights=left_p, k=1)[0]
	right_node = random.choices(range(len(right_p)), weights=right_p, k=1)[0]

	if left_node == 0:
		left_marker = [0]*(left.count('(')+1)
	if left_node == 2:
		left_marker = [1]*(left.count('(')+1)
	if right_node == 0:
		right_marker = [0]*(right.count('(')+1)
	if right_node == 2:
		right_marker = [1]*(right.count('(')+1)
	if (left_node == 1) & ('(' in left):
		left_marker = generate_marker(left,tvalue[1:left.count('(')+1])
	elif (left_node == 1) & (mapping[left] is not None):
		left_marker = ["01"]*(left.count('(')+1)
	if (right_node == 1) & ('(' in right):
		right_marker = generate_marker(right,tvalue[left.count('(')+1:])
	elif (right_node == 1) & (mapping[right] is not None):
		right_marker = ["01"]*(right.count('(')+1)


	marker = left_marker+right_marker

	return marker


if input_markerfile == "":
	pass
else:
	with open(input_markerfile, mode='r', newline='') as infile:
		reader = csv.reader(infile)

		header = next(reader)
		new_order_indices = [header.index(col) for col in taxa_list]
		reordered_data = [taxa_list]  # 首先加入重新排列后的列名

		# 重新排列每一行的数据并加入到结果列表中
		for row in reader:
			reordered_row = [row[i] for i in new_order_indices]
			reordered_data.append(reordered_row)

	question_location=[]
	for i in range(1,len(reordered_data)):
		marker = reordered_data[i]
		positions = [index for index, value in enumerate(marker) if value == '?']
		question_location.append(positions)


bracketloc = sorted(find_matching_parentheses(treetopo))
cluster = [treetopo[i[0]:i[1]+1] for i in bracketloc]
probabilities = [a * b for a, b in zip(t, c)]



for example_n in range(num_samples):

	markerlist = []
	while len(markerlist) < num_marker:
		selected_index = random.choices(range(len(probabilities)), weights=probabilities, k=1)[0]
		#print(selected_index)

		insert_edge = cluster[selected_index]
		#selected_taxa = re.findall(r"taxa\d+", insert_edge)
		selected_taxa = insert_edge.replace('(','').replace(')','').split(',')  # Extract taxa within the cluster

		insert_rate = [c[selected_index]*(t[selected_index]/2-(1-math.exp(-t[selected_index]))/2),c[selected_index]*(1-math.exp(-t[selected_index])),c[selected_index]*(t[selected_index]/2-(1-math.exp(-t[selected_index]))/2)]
		insert = random.choices(range(len(insert_rate)), weights=insert_rate, k=1)[0]

		marker = ["?" if taxa in selected_taxa else 0 for taxa in taxa_list]
		#print("a",marker)

		if (insert == 0) & (selected_index>0):
			last_common_ancestor = [0]*len(selected_taxa)
		elif (insert == 2) & (selected_index>0):
			last_common_ancestor = [1]*len(selected_taxa)
		else:
			last_common_ancestor = generate_marker(insert_edge,t[selected_index:selected_index+insert_edge.count('(')])

		start = marker.index("?")  # 第一个 `?` 的位置
		end = len(marker) - marker[::-1].index("?")  # 最后一个 `?` 的位置
		result = marker[:start] + last_common_ancestor + marker[end:]  # 替换中间的 `?`
		#print(result,result.count(1),result.count(0))
		
		if input_markerfile != "":
			qm = question_location[len(markerlist)]
			for pos in qm:
				result[pos] = '?'

		if (result.count(1)>1) & (result.count(0)>0):
			markerlist.append(result)


	file_path = f'simulate_marker_{example_n+1}.csv'
	with open(file_path, "w", newline="") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerows([taxa_list])
		writer.writerows(markerlist)

print(f"Successfully generated {num_samples} marker files.")

