import matlab
import matlab.engine
import subprocess


markerfile = input("Please input the marker file: ")
newicktree = input("Please input the tree file: ")
usec = input("use parameter c (yes/no): ").lower()
useq = input("use parameter q (yes/no): ").lower()


alltreetopo = []
with open(newicktree, 'r') as mark:
    for line in mark:
        alltreetopo.append(line.strip())
#print(alltreetopo)


for treetopo in alltreetopo:
    print('tree:',treetopo)
    if useq =='yes':
        subprocess.run(["python3", "newlambda_q.py", treetopo, markerfile, usec])
    else:
        subprocess.run(["python3", "newlambda_cpp.py", treetopo, markerfile, usec])
    eng = matlab.engine.start_matlab()
    F0=eng.treescore(nargout =1)
    eng.exit()

