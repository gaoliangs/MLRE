import matlab
import matlab.engine
import re
from operator import itemgetter
import subprocess
import os
import itertools
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
#print(findsplit('a,b,(c,d)'))



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
#print(findunresolvedsplit('((a,b,(c,d)),(e,f,g))'))



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


######### split in tree

besttreetopo = input("starting tree: ")
markerfile = input("input file name: ")
option = input("use parameter c (yes/no): ").lower()
useq = input("use parameter q (yes/no): ").lower()
searchtree = input("constraint tree:")
'''
besttreetopo = '((((((((((((FULGL,PYGAD_APTFO),(((EGRGA,PELCR),NIPNI),PHACA)),GAVST),(CAPCA,CHAPE_CALAN)),((OPHHO,(EURHE,PHALE)),CHAVO)),((COLST,(((CATAU,HALLE_HALAL),TYTAL),(LEPDI,(APAVI,(BUCRH,PICPU_MERNU))))),(((TAEGU_MANVI_GEOFO_CORBR_ACACH,NESNO_MELUN),FALPE),CARCR))),BALRE),(CHLUN,TAUER)),CUCCA),PODCR_PHORU),COLLI),(PTEGU,MESUN))'
markerfile = 'newseq1.csv'
option='yes'
useq = 'no'
searchtree='((CAPCA,CHAPE_CALAN),(EURHE,PHALE),(MESUN,PTEGU),(((FULGL,PYGAD_APTFO),(((EGRGA,PELCR),NIPNI),PHACA)),GAVST),((((NESNO_MELUN,TAEGU_MANVI_GEOFO_CORBR_ACACH),FALPE),CARCR),((((CATAU,HALLE_HALAL),TYTAL),(((BUCRH,PICPU_MERNU),APAVI),LEPDI)),COLST)),TAUER,COLLI,CHLUN,CHAVO,BALRE,PODCR_PHORU,CUCCA,OPHHO)'
'''

bracketloc = find_matching_parentheses(searchtree)
cluster_unchange = [searchtree[i[0]:i[1]+1].replace(')','').replace('(','').split(',') for i in bracketloc]


treetopo=''
eng = matlab.engine.start_matlab()

while treetopo != besttreetopo:
	treetopo = besttreetopo


	##NNI
	print('start NNI')
	print('topo:',treetopo)

	if useq =='yes':
		subprocess.run(["python3", "newlambda_q.py", treetopo, markerfile, option])
	else:
		subprocess.run(["python3", "newlambda_cpp.py", treetopo, markerfile, option])
	#eng = matlab.engine.start_matlab()
	F0=eng.treescore(nargout =1)
	#eng.exit()

	bracketloc = sorted(find_matching_parentheses(treetopo))

	edgescore = [1]*(len(bracketloc)-1)
	while edgescore != [len(bracketloc)]*(len(bracketloc)-1):
		#找到score最小的那条边及括号位置
		#edge_index = edgescore.index(min(edgescore))+1

		minindex = min(edgescore)
		edge_index = [i for i in range(len(edgescore)) if edgescore[i] == minindex]
		for edge in edge_index:
			loc = bracketloc[edge+1]
			newsplit = treetopo[loc[0]:loc[1]+1]
			newsplit = newsplit.replace('(','').replace(')','').split(',')
			edge_unchange=False
			for bs in cluster_unchange:
				if set(bs) == set(newsplit):
					edgescore[edge]=len(bracketloc)
					print(edgescore)
					edge_unchange=True
					break
				
			if edge_unchange:
				continue
			else:
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


				F2 = []
				for tt in newtreetopo:
					print('topo:',tt)
					if useq =='yes':
						subprocess.run(["python3", "newlambda_q.py", tt, markerfile, option])
					else:
						subprocess.run(["python3", "newlambda_cpp.py", tt, markerfile, option])
					#eng = matlab.engine.start_matlab()
					F=eng.treescore(nargout =1)
					#eng.exit()

					F2.append(F)
			    

				if max(F2) > F0:
					F0 = max(F2)
					print('update tree',F0)
					treetopo = newtreetopo[0] if F2[0] >= F2[1] else newtreetopo[1]
					bracketloc = sorted(find_matching_parentheses(treetopo))

					for j in range(len(edgescore)):
						ancestor = find_last_common_ancestor(bracketloc,bracketloc[j+1],bracketloc[edge+1])
						pathnum = cluster_level(treetopo,bracketloc[j+1])+cluster_level(treetopo,bracketloc[edge+1])-2*cluster_level(treetopo,ancestor)
						edgescore[j] = pathnum
					edgescore[edge] = len(bracketloc)
				else:
					edgescore[edge]= len(bracketloc)

			print(edgescore)


	print('NNItree:',treetopo)
	print('NNI_F:',F0)
	besttreetopo = treetopo
	bestF =F0


	##spr
	print('start SPR')
	topospr = besttreetopo
	taxa = topospr.replace("(", "").replace(")", "").split(',')

	bracketloc = sorted(find_matching_parentheses(topospr))
	clustera = [topospr[i[0]:i[1]+1].replace(')','').replace('(','').split(',') for i in bracketloc]
	clustera_sorted = [sorted(sublist) for sublist in clustera]
	cluster_unchange_sorted = [sorted(sublist) for sublist in cluster_unchange]
	complement_indices = [i for i, sublist in enumerate(clustera_sorted) if sublist not in cluster_unchange_sorted]

	cluster_topo = [topospr[i[0]:i[1]+1] for i in bracketloc]
	non_bunemacluster = [cluster_topo[i] for i in complement_indices]
	#print(non_bunemacluster)

	buneman_split=findsplit(searchtree[1:-1])
	buneman_taxa_sorted = [sorted(bs.replace(')','').replace('(','').split(',')) for bs in buneman_split]
	complement_indices = [i for i, sublist in enumerate(clustera_sorted) if sublist in buneman_taxa_sorted]
	buneman_cluster = [cluster_topo[i] for i in complement_indices]+[buneman_taxa for buneman_taxa in buneman_split if buneman_taxa.count('(')==0]
	#print(buneman_cluster)

	replacements = [f"split{i+1}" for i in range(len(buneman_cluster))]

	non_bc_split=[]
	for non_b in non_bunemacluster:
		split =False
		for i, item in enumerate(buneman_cluster):
			if item in non_b:
				split = True
			non_b = non_b.replace(item, replacements[i])
		if split: 		
			non_bc_split.append(non_b)
	#print(non_bc_split)

	for i, item in enumerate(buneman_cluster):
		topospr = topospr.replace(item, replacements[i])

	allspr = []
	for i in non_bc_split+replacements:
		xtopo = topospr.replace(i,'X')
		xspr = spr_x(xtopo)
		spr = [j.replace('X',i) for j in xspr]
		allspr.append(spr)

	allspr = [element for i in allspr for element in i]
	#print(allspr)

	allsprtree = []
	for ast in allspr:
		for i, item in enumerate(buneman_cluster):
			ast = re.sub(rf"\bsplit{i+1}\b", item, ast)
		allsprtree.append(ast)

	#print(len(allsprtree))


	for tree in allsprtree:
		print('topo:',tree)
		if useq =='yes':
			subprocess.run(["python3", "newlambda_q.py", tree, markerfile, option])
		else:
			subprocess.run(["python3", "newlambda_cpp.py", tree, markerfile, option])
		#eng = matlab.engine.start_matlab()
		F=eng.treescore(nargout =1)
		#eng.exit()

		if F>bestF:
			besttreetopo = tree
			bestF =F
			break

	print('spr_tree:',besttreetopo)
	print('spr_F:',bestF)


eng.exit()

