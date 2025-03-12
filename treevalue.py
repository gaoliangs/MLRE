import matlab
import matlab.engine
import subprocess
import re
import io
import time
import os

def is_valid_path(path):
    return len(path) < 255 and os.path.exists(path)

user_input = input('Please input the tree(Newick string or filename): ')

if is_valid_path(user_input):
    with open(user_input, 'r') as file:
        alltreetopo = [line.strip() for line in file]
else:
    alltreetopo = [user_input]
    
alltree =[]
for newick_str in alltreetopo:
    newick_str = newick_str.replace(" ", "").replace(";", "")
    # remove edge length
    newick_str = re.sub(r':\d+(\.\d+)?([eE][+-]?\d+)?', '', newick_str)
    # remove bootstrap
    newick_str = re.sub(r'\)(\d+(\.\d+)?([eE][+-]?\d+)?)', ')', newick_str)
    alltree.append(newick_str)

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
