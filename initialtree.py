import sys
import itertools
import math


file_name = sys.argv[1]
with open(file_name, 'r') as csv_file:
    csv_content = csv_file.readlines()

bunemantree = sys.argv[2]

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

bunemancluster = [bunemantree[start:end+1] for start, end in find_matching_parentheses(bunemantree)[1:]]
taxa_in_buneman = [bc.replace('(','').replace(')','').split(',') for bc in bunemancluster]


taxaname = csv_content[0].strip().split(',')
taxanum = len(taxaname)

data = [line.strip().split(',') for line in csv_content[1:]]
taxa_data=list(map(list, zip(*data)))



triplet = []
for j in itertools.combinations(list(range(taxanum)), 3):
	na =0
	nb =0
	nc =0
	for i in range(len(taxa_data[0])):
		if (taxa_data[j[0]][i] == '1')&(taxa_data[j[1]][i] == '1')&(taxa_data[j[2]][i] == '0'):
			nc+=1
		if (taxa_data[j[0]][i] == '1')&(taxa_data[j[1]][i] == '0')&(taxa_data[j[2]][i] == '1'):
			nb+=1
		if((taxa_data[j[0]][i] == '0')&(taxa_data[j[1]][i] == '1')&(taxa_data[j[2]][i] == '1')):
			na+=1
	triplet.append([na,nb,nc])
#print(triplet)

#marker = list(itertools.combinations(taxaname, 3))
#m_qw = list(zip(marker,triplet))
#print(m_qw)


def findindex(u,v,w):
	qu = [u,v,w]
	qu.sort()
	qindex = math.comb(taxanum, 3) - math.comb(taxanum - qu[0] - 1, 3) - math.comb(taxanum - qu[1] - 1, 2) - math.comb(taxanum - qu[2] - 1, 1)
	if (set([u,v]) == set([qu[0],qu[1]])):
		orindex = 2
	if (set([u,v]) == set([qu[0],qu[2]])):
		orindex = 1
	if (set([u,v]) == set([qu[1],qu[2]])):
		orindex = 0
	return([qindex,orindex])
#print('i',findindex(2,1,0))



def find_valid_newpair(score, a, b, taxa_in_buneman):
	while True:
		newpairloc = score.index(max(score))
		newpair = a[newpairloc] + b[newpairloc]

		for tb in taxa_in_buneman:
			if (len(set(tb)&set(newpair)) >0) & (len(set(tb)&set(newpair))<len(newpair)) & (len(set(tb)&set(newpair))<len(tb)): 
				score[newpairloc] = -float('inf') 
				break
		else:
			return newpairloc


index_map = {name: idx for idx, name in enumerate(taxaname)}

cluster = taxaname

newpairlist =[]
while len(cluster)>2:

	pair = list(itertools.combinations(cluster, 2))
	a = [t[0].split(',') for t in pair]
	b = [t[1].split(',') for t in pair]
	rest = [sorted(list(set(cluster).difference(set(j)))) for j in pair]

	score = []
	for i in range(len(pair)):
		aa = a[i]
		bb = b[i]
		qw = 0
		d = 0
		
		for cc in rest[i]:
			#print(cc)
			c1 = cc.split(',')
			d = len(aa)*len(bb)*len(c1)
			tw = 0
			for ci in c1:
				for ai in aa:
					for bi in bb:
						qn = findindex(index_map[ai], index_map[bi], index_map[ci])
						#tw = tw+float(m_qw[qn[0]-1][1][qn[1]])
						tw = tw+float(triplet[qn[0]-1][qn[1]])
			qw = qw +tw/d
		
		score.append(qw)


	#newpairloc = score.index(max(score))  # 找到最大得分的位置
	newpairloc = find_valid_newpair(score,a,b,taxa_in_buneman)
	newpair = a[newpairloc] + b[newpairloc]  # 获取对应的newpair
	#print(newpair)

	
	newpairc = ','.join(newpair)
	cluster = rest[newpairloc]
	cluster.append(newpairc)
	cluster.sort()
	newpairlist.append(newpairc)

	
fullcluster = ','.join(cluster)
newpairlist = newpairlist+[fullcluster]
newpairlist = list(reversed(newpairlist))


for i in range(len(newpairlist)):
	for j in range(i+1,len(newpairlist)):
		if newpairlist[j] in newpairlist[i]:
			newpairlist[i] = newpairlist[i].replace(newpairlist[j],'('+newpairlist[j]+')')


tjtree = '('+newpairlist[0]+')'

print(tjtree)


