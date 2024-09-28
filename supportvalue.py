import matlab
import matlab.engine
import re
from operator import itemgetter
import subprocess
import math


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



def findunresolvedsplit(tree):
	brloc = find_matching_parentheses(tree)
	for br in brloc:
		unsp = tree[br[0]:br[1]+1]
		csp1 = tree[:br[0]]
		csp2 = tree[br[1]+1:]
		cn = unsp.count(',')
		bn = unsp.count('(')
		#print(unsp,cn,bn)
		if cn-bn>0:
			return([unsp,csp1,csp2])
			break

def isthesame(str1,str2):
	b1 = find_matching_parentheses(str1)
	b2 = find_matching_parentheses(str2)
	clu1 = [sorted(str1[bb1[0]:bb1[1]+1].replace('(', '').replace(')', '').split(',')) for bb1 in b1]
	clu2 = [sorted(str2[bb2[0]:bb2[1]+1].replace('(', '').replace(')', '').split(',')) for bb2 in b2]
	return sorted(clu1) == sorted(clu2)



def find_last_common_ancestor(bloc,pair1,pair2):
	for start,end in sorted(bloc,key =itemgetter(1)):
		if (start<=pair1[0]) & (start<=pair2[0]) & (end>=pair1[1]) & (end>=pair2[1]):
			return [start,end]
			break


def cluster_level(s,pair):
	return s[:pair[0]].count('(')-s[:pair[0]].count(')')
#print(cluster_level(treetopo,[6, 18]))




def spr_x(s):
    #get all clusters in tree string
    bracket = find_matching_parentheses(s)
    clu = [s[i[0]:i[1]+1] for i in bracket]
    taxon = s.replace("(", "").replace(")", "").split(',')
    fullclu = clu+taxon
    
    for c in clu:
        two_split = findsplit(c[1:-1])
        if ('X' == two_split[0])|('X' == two_split[1]):
            original_c = c
            newc = [two_split[0] if two_split[1]== 'X' else two_split[1]]
            fullclu.remove(newc[0])
            fullclu.remove('X')
            fullclu.remove(c)
    #print(fullclu)

    delete_x = s.replace(original_c,newc[0])
    clu_delete_x = [i.replace(original_c,newc[0]) for i in fullclu]
    #print(clu_delete_x)

    new = [f'({i},X)' for i in clu_delete_x]
    #print('sdfgsdgf',clu_delete_x,new)
    #de_replaced = re.sub(r'\bta1\b', , delete_x)

    sprtree = []
    for i in range(len(new)):
        if '(' in clu_delete_x[i]:
            sprtree.append(delete_x.replace(clu_delete_x[i],new[i]))
        else:
            sprtree.append(re.sub(r'\b' + re.escape(clu_delete_x[i]) + r'\b',new[i], delete_x))
        #print(sprtree)

    return sprtree

def remove_duplicates(trees):
    unique_trees = []

    for i, tree1 in enumerate(trees):
        is_duplicate = False

        for j, tree2 in enumerate(trees[i + 1:]):
            if isthesame(tree1, tree2):
                is_duplicate = True
                break

        if not is_duplicate:
            unique_trees.append(tree1)

    return unique_trees



treetopo = input("starting tree: ")
markerfile = input("input file name: ")
option = input("use parameter c (yes/no): ").lower()
useq = input("use parameter q (yes/no): ").lower()


eng = matlab.engine.start_matlab()

bracketloc = sorted(find_matching_parentheses(treetopo))
#print(bracketloc)

if useq =='yes':
	subprocess.run(["python3", "newlambda_q.py", treetopo, markerfile, option])
else:
	subprocess.run(["python3", "newlambda.py", treetopo, markerfile, option])
#eng = matlab.engine.start_matlab()
F0=eng.treescore(nargout =1)



NNI = []
for edge in range(len(bracketloc)-1):
	loc = bracketloc[edge+1]

	newtopo = treetopo[:loc[0]]+treetopo[loc[0]+1:loc[1]]+treetopo[loc[1]+1:]
	unresolvedtree = findunresolvedsplit(newtopo)
	newtree = findsplit(unresolvedtree[0][1:-1])

	newtreetopo1 = unresolvedtree[1]+'(('+newtree[0]+','+newtree[1]+'),'+newtree[2]+')'+unresolvedtree[2]
	newtreetopo2 = unresolvedtree[1]+'(('+newtree[0]+','+newtree[2]+'),'+newtree[1]+')'+unresolvedtree[2]
	newtreetopo3 = unresolvedtree[1]+'(('+newtree[1]+','+newtree[2]+'),'+newtree[0]+')'+unresolvedtree[2]

	newtreetopo =[]
	if not isthesame(newtreetopo1,treetopo):
		newtreetopo.append(newtreetopo1)
	if not isthesame(newtreetopo2,treetopo):
		newtreetopo.append(newtreetopo2)
	if not isthesame(newtreetopo3,treetopo):
		newtreetopo.append(newtreetopo3)

	NNI.append(newtreetopo)

print(NNI)

NNIscore=[]
support=[]
for edge in range(len(NNI)):
	NNItopo = NNI[edge][0]
	if useq =='yes':
		subprocess.run(["python3", "newlambda_q.py", NNItopo, markerfile, option])
	else:
		subprocess.run(["python3", "newlambda.py", NNItopo, markerfile, option])
	#eng = matlab.engine.start_matlab()
	F1=eng.treescore(nargout =1)
	#eng.exit()

	NNItopo = NNI[edge][1]
	if useq =='yes':
		subprocess.run(["python3", "newlambda_q.py", NNItopo, markerfile, option])
	else:
		subprocess.run(["python3", "newlambda.py", NNItopo, markerfile, option])
	#eng = matlab.engine.start_matlab()
	F2=eng.treescore(nargout =1)

	NNIscore.append([F1,F2])
	result = 1 / (1 + math.exp(F1 - F0) + math.exp(F2 - F0))
	support.append(result)

print(support)

