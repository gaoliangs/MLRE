import matlab
import matlab.engine
import subprocess
import re
import io
import time
import os
from io import StringIO
from Bio import Phylo


user_input = input('Please input the tree(Newick string or filename): ')

def remove_information(tree):
    for clade in tree.find_clades():
        clade.branch_length = None  # 去掉枝长
        clade.comment = None     # 去掉注释
        clade.confidence = None  # 去掉支持值
    return tree


if os.path.exists(user_input):
	with open(user_input, 'r') as file:
		alltreetopo = [line.strip() for line in file]
else:
	alltreetopo = [user_input]

alltree = []
for newick_str in alltreetopo:
	tree = Phylo.read(StringIO(newick_str), "newick")
	topology_str = remove_information(tree).format("newick")
	str_noedge = re.sub(r":\d+(\.\d+)?", "", topology_str)
	alltree.append(str_noedge.strip().replace(";", ""))


markerfile = input("Please input the marker filename: ")
usec = input("use parameter c (yes/no): ").lower()
useq = input("use parameter q (yes/no): ").lower()


start_time = time.time()

for tree in alltree:
	print('tree:',tree)
	if useq =='yes':
		subprocess.run(["python3", "newlambda_q.py", tree, markerfile, usec])
		eng = matlab.engine.start_matlab()
		stdout_capture = io.StringIO()
		LogL,t,c,q=eng.treescore(nargout =4, stdout=stdout_capture)
		print("LogL = ", LogL)
		print("t = ", t[0])
		print("c = ", c[0])
		print("q = ", q[0])
		eng.exit()
	else:
		subprocess.run(["python3", "newlambda_cpp.py", tree, markerfile, usec])
		eng = matlab.engine.start_matlab()
		stdout_capture = io.StringIO()
		LogL, t, c =eng.treescore(nargout =3, stdout=stdout_capture)
		print("LogL = ", LogL)
		print("t = ", t[0])
		print("c = ", c[0])
		eng.exit()



end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time:.4f} seconds")
