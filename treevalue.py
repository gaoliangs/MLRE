import matlab
import matlab.engine
import subprocess
import re


user_input = input('Please input the tree(Newick string or filename): ')
# 检查用户输入是否是文件名，如果是，则读取文件内容
try:
    with open(user_input, 'r') as file:
        alltreetopo = [line.strip() for line in file]
except FileNotFoundError:
    # 如果输入不是文件名，则直接使用用户输入的字符串
    alltreetopo = [user_input]
    
alltree =[]
for newick_str in alltreetopo:
    newick_str = newick_str.replace(" ", "").replace(";", "")
    newick_str = re.sub(r':\d+(\.\d+)?', '', newick_str)
    newick_str = re.sub(r'\)(\d+(\.\d+)?)', ')', newick_str)
    alltree.append(newick_str)

markerfile = input("Please input the marker filename: ")
usec = input("use parameter c (yes/no): ").lower()
useq = input("use parameter q (yes/no): ").lower()



for tree in alltree:
    print('tree:',tree)
    if useq =='yes':
        subprocess.run(["python3", "newlambda_q.py", tree, markerfile, usec])
    else:
        subprocess.run(["python3", "newlambda_cpp.py", tree, markerfile, usec])
    eng = matlab.engine.start_matlab()
    F0=eng.treescore(nargout =1)
    eng.exit()

