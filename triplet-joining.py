import os
import itertools
import math

current_directory = os.getcwd()

file_name = input("input file name: ")
#file_name = 'newseq1.csv'
input_file_path = os.path.join(current_directory, file_name)
with open(input_file_path, 'r') as csv_file:
    csv_content = csv_file.readlines()

data =[]
for i in csv_content[1:]:
	data.append(i[:-1].split(','))
taxa_data=list(map(list, zip(*data)))



taxaname = csv_content[0][:-1].split(',')
taxanum = len(taxaname)
print(taxaname)


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

print(triplet)


tjtree =[]



marker = []
for j in itertools.combinations(taxaname, 3):
   #mm = list(j)+['R']
   marker.append(list(j))

m_qw = list(zip(marker,triplet))
#print(m_qw)


def c(n,m):
	if n >=m:
		return(math.factorial(n)/(math.factorial(m)*math.factorial(n-m)))
	else:
		return(0)



def findindex(u,v,w):
	qu = [u,v,w]
	qu.sort()
	qindex = c(taxanum,3) - c(taxanum-qu[0]-1,3)-c(taxanum-qu[1]-1,2)-c(taxanum-qu[2]-1,1)
	if (set([u,v]) == set([qu[0],qu[1]])):
		orindex = 2
	if (set([u,v]) == set([qu[0],qu[2]])):
		orindex = 1
	if (set([u,v]) == set([qu[1],qu[2]])):
		orindex = 0
	return([int(qindex),orindex])

#print('i',findindex(2,1,0))



cluster = taxaname

#cluster= ['A','B','C','D','E','F','G','H','I','J','K','L','M']
#cluster = ['APAVI', 'BALRE', 'BUCRH', 'CAPCA', 'CARCR', 'CATAU', 'CHAPE_CALAN', 'CHAVO', 'CHLUN', 'COLST', 'COLLI', 'CUCCA', 'EGRGA', 'EURHE', 'FALPE', 'FULGL', 'GAVST', 'HALLE_HALAL', 'LEPDI', 'MESUN', 'NESNO_MELUN', 'NIPNI', 'OPHHO', 'PELCR', 'PHALE', 'PHACA', 'PICPU_MERNU', 'PODCR_PHORU', 'PTEGU', 'PYGAD_APTFO', 'TAEGU_MANVI_GEOFO_CORBR_ACACH', 'TAUER', 'TYTAL']
#cluster.sort()
#print(cluster)

newpairlist =[]
clusterlist =[]
while len(cluster)>2:

	pair = []
	a = []
	b = []
	rest = []
	resttaxa = []
	for j in itertools.combinations(cluster, 2):
	    pair.append(j)
	for t in pair:
		a.append(t[0].split(','))
		b.append(t[1].split(','))



	for j in pair:
		rest.append(sorted(list(set(cluster).difference(set(j)))))


	score = []
	for i in range(len(pair)):
		aa = a[i]
		bb = b[i]
		qw = 0
		d = 0
		#if (bb == ['MELUN_NESNO']) & (aa == ['FULGL']):
		#	print('sdddd')
		
		for cc in rest[i]:
			#print(cc)
			c1 = cc.split(',')
			#if aa == ['APTFO','PYGAD']:
		
			#d = d+len(aa)*len(bb)*len(c1)
			d = len(aa)*len(bb)*len(c1)
			tw = 0
			for ci in c1:
				for ai in aa:
					for bi in bb:
						qn = findindex(taxaname.index(ai),taxaname.index(bi),taxaname.index(ci))
						tw = tw+float(m_qw[qn[0]-1][1][qn[1]])
			qw = qw +tw/d
		#if (aa == ['NESNO_MELUN']) & (bb == ['TAEGU_MANVI_GEOFO_CORBR_ACACH']):
		#	print('sddd',qw)
		
		score.append(qw)


	#print(score)
	newpairloc = score.index(max(score))
	#print(score)
	#print(max(score))
	#print(newpairloc)
	newpair = a[newpairloc]+b[newpairloc]
	#print(newpair)
	#print(rest[newpairloc])
	
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



tjtree.append('('+newpairlist[0]+')')

print(tjtree)

file_name = 'tjtree_newick.txt'
output_file_path = os.path.join(current_directory, file_name)
file = open(output_file_path, "w+")
for i in tjtree:
    file.write(i)
    file.write('\r\n')
file.close()


