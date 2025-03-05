#preprocessing
import os
import itertools 
import math
import numpy as np
from scipy.stats import binomtest
from collections import OrderedDict


# 读取NEXUS文件内容
file_name = input("input file name: ")
#file_name = 'Neoaves_2118_44.nex'
with open(file_name, 'r') as nexus_file:
    nexus_content = nexus_file.read()

# 查找数据矩阵的开始和结束位置
matrix_start = nexus_content.lower().find("matrix") + len("matrix")
matrix_end = nexus_content.find(";", matrix_start)

# 提取数据矩阵部分
matrix_data = nexus_content[matrix_start:matrix_end].strip()


# 读取每一行的名字和内容
rows = matrix_data.split('\n')
taxanames = []
marker = []
for row in rows:
    taxanames.append(row.split(' ')[0])
    marker.append(row.split(' ')[-1])
taxanum = len(taxanames)


###conflict free marker
def check_conflict(listA, listB):
    for a, b in zip(listA, listB):
        if (a == '1' and b not in ['1', "?"]) or (a == '0' and b not in ['0', "?"]):
            return True  # conflict
    return False  # no conflict


def replace_question_mark(listA, listB):
    listC = [a if a != "?" else b for a, b in zip(listA, listB)]
    return "".join(listC)


    
# merge non-conflict marker
new_data = []
new_name = []
i = 0
while i < len(marker):
    k = 0
    for j in range(i+1,len(marker)) :
        if check_conflict(marker[i], marker[j]):
            k +=1
        else:
            merged = replace_question_mark(marker[i], marker[j])
            marker[j]=merged
            taxanames[j] = str(taxanames[j]+'_'+taxanames[i])
            break
    if k == len(marker)-i-1:
        new_data.append(marker[i])
        new_name.append(taxanames[i])
    i += 1



###find outgroup taxa
def filter_lists(names, binaries):
    filtered_names = []
    filtered_binaries = []

    for name, binary in zip(names, binaries):
        if '1' in binary:
            filtered_names.append(name)
            filtered_binaries.append(binary)

    return filtered_names, filtered_binaries


taxaname, taxa_data = filter_lists(new_name, new_data)
taxanum = len(taxaname)
print('non-cinflict taxa:',taxaname)



# triplet weight
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



# buneman weight
weight = triplet
buneman_weight = [[0,0,0] for i in range(len(triplet))]
for i in range(len(triplet)):
    for j in range(3):
        if weight[i][j] == max(weight[i]):
            #buneman_weight[i][j] = weight[i][j] - sum(weight[i])+max(weight[i])+min(weight[i])
            result = binomtest(max(weight[i]), sum(weight[i])-min(weight[i]), 0.5, alternative='less')
            buneman_weight[i][j] = result.pvalue
        else:
            buneman_weight[i][j] = 0
    

# 寻找特定triplet在buneman weight矩阵中的位置
def find(a, b, c):
    s = sorted([a, b, c])
    row = math.comb(taxanum, 3) - math.comb(taxanum - s[0], 3) - math.comb(taxanum - s[1], 2) - math.comb(taxanum - s[2], 1)
    if c == min(a, b, c):
        return [row-1, 0]
    if c == max(a, b, c):
        return [row-1, 2]
    if c == sum([a, b, c]) - max(a, b, c) - min(a, b, c):
        return [row-1, 1]

# find first buneman cluster
indexnum = [max(buneman_weight[0])]
clusterB = [buneman_weight[0].index(indexnum[0])+1]
clusterA = list(set([1,2,3])-set(clusterB))
cluster=[[clusterA,clusterB]]



# find buneman cluster
for i in range(4,taxanum+1):
    newcluster =[]
    newindex =[]

    for s in range(len(indexnum)):
        num = indexnum[s]
        clusterA = cluster[s][0]
        clusterB = cluster[s][1]
      
        # i in clusterA
        for x in range(len(clusterA)):
            a = clusterA[x]
            b = i
            for y in range(len(clusterB)):
                c = clusterB[y]
                weight_aib = buneman_weight[find(a, b, c)[0]][find(a, b, c)[1]]
                num = min(num, weight_aib)
        if num>0:
            newindex.append(num)
            clustera = clusterA + [i]
            newcluster.append([clustera,clusterB])

        # i in clusterB
        num = indexnum[s]
        for x in range(len(clusterA) - 1):
            for y in range(x + 1, len(clusterA)):
                a = clusterA[x]
                b = clusterA[y]
                c = i
                weight_aai = buneman_weight[find(a, b, c)[0]][find(a, b, c)[1]]
                num = min(num, weight_aai)
        if num>0:
            newindex.append(num)
            clusterb = clusterB + [i]
            newcluster.append([clusterA,clusterb])

    # i & 1...i-1
    for x in range(1, i):
        a = x
        b = i
        setc = list(set(range(1, i)) - set([x]))
        num = []
        for y in range(len(setc)):
            c = setc[y]
            num.append(buneman_weight[find(a, b, c)[0]][find(a, b, c)[1]])
        if 0 not in num:
            newindex.append(min(num))
            clustera = [a,i]
            clusterb = setc
            newcluster.append([clustera,clusterb])

    # 1..i-1|i
    num = []
    for x in range(1,i):
        for y in range(x + 1, i):
            a = x
            b = y
            c = i
            num.append(buneman_weight[find(a, b, c)[0]][find(a, b, c)[1]])      
    if 0 not in num:
        newindex.append(min(num))
        clustera = list(range(1, i))
        clusterb = [i]
        newcluster.append([clustera,clusterb])

    cluster = newcluster
    indexnum = newindex


# print buneman cluster, buneman score
for i,c in enumerate(cluster):
    bc= []
    for j in c[0]:
        bc.append(taxaname[j-1])
    print(bc,indexnum[i])

#is number
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		pass
	return False



#threshold
user_input = input("choose threshold: ")
bunemancluster = []
if is_number(user_input):
	for i in range(len(indexnum)):
		if indexnum[i]> float(user_input):
			bunemancluster.append(cluster[i][0])
else:
	user_input = max(indexnum)+1
	for i in range(len(indexnum)):
		if indexnum[i]> float(user_input):
			bunemancluster.append(cluster[i][0])


bunemancluster = [[taxaname[index-1] for index in sub_list] for sub_list in bunemancluster]


# full set
X = taxaname
# find child
def find_child(full,B):
    B.sort(key=len)
    subset = []
    restB =[]

    for i, a in enumerate(B):
        is_subset = False
        for j, b in enumerate(B[i + 1:]):
            if set(a).issubset(set(b)):
                restB.append(a)
                is_subset = True
                break
        if not is_subset:
            full = [item for item in full if item not in a]
            subset.append(a)


    child = []
    for j in subset:
        childj = []
        for jj in restB:
            if set(jj).issubset(set(j)):
                childj.append(jj)
        child.append(childj)

    return [subset+full,child]


#print(X,bunemancluster)
root = find_child(X,bunemancluster)[0]
restbuneman = find_child(X,bunemancluster)[1]


unresolved = []
if len(root)>2:
	unresolved.append(root)


####
while len(restbuneman) >0:
    newroot=[]
    newbuneman=[]
    for i,a in enumerate(restbuneman):
        if len(root[i])>2:

            clu = find_child(root[i],a)[0]
            rest = find_child(root[i],a)[1]

            if len(clu)>2:
                unresolved.append(clu)

        
            for j,b in enumerate(rest):
                if b != []:
                    newroot.append(clu[j])
                    newbuneman.append(b)
                if (len(clu[j])>2) & (len(b)==0):
                    unresolved.append(clu[j])

    root = newroot
    restbuneman=newbuneman



#生成group marker
def merge_marker(group,taxaname,data):
	seq =[]
	for g in group:
		loc = taxaname.index(g)
		seq.append(data[loc])


	merged_sequence = [
    '?' if set(chars) == {'?'} else
    '1' if '0' not in chars else
    '0' if '1' not in chars else
    'B' if set(chars) in ({'0', '1'}, {'0', '1', '?'}) else
    '?' for chars in zip(*seq)
	]

	return ''.join(merged_sequence)


#######################################################################

# 输出markers
for j in range(len(unresolved)):
    sequences =[]
    for gi in unresolved[j]:
        if type(gi) == list:
            sequences.append(merge_marker(gi,taxaname,taxa_data))
        if type(gi) == str:
            sequences.append(taxa_data[taxaname.index(gi)])


    transposed_sequences = list(zip(*sequences))

    # 保留条件：不满足删除条件的索引
    indices_to_keep = OrderedDict()
    for i, chars in enumerate(transposed_sequences):
        if not (chars.count('1') + chars.count('?') == len(chars)) and \
           not (chars.count('0') + chars.count('?') == len(chars)) and \
           not (chars.count('1') + chars.count('B') < 2) and \
           not (chars.count('1') + chars.count('B') + chars.count('0') < 3):
            indices_to_keep[i] = None


    # 删除不满足条件的字符
    filtered_sequences = [
        ''.join(chars[i] for i in indices_to_keep.keys())
        for chars in sequences
    ]
    marker = list(zip(*filtered_sequences))


    #unresolved data
    output_name = 'newseq'+str(j+1)+'.csv'
    with open(output_name, 'w') as output_file:
        names = []
        for gi in unresolved[j]:
            if type(gi) == list:
                names.append(str('_'.join(gi)))
            if type(gi) == str:
                names.append(gi)
        output_file.write(','.join(names) + "\n")
        for item in marker:
            output_file.write(','.join(item)+ "\n")









