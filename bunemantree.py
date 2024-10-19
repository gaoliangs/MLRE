import os
import itertools
import math

#buneman tree
current_directory = os.getcwd()
file_name=input("marker file: ")
current_directory = os.getcwd()
input_file_path = os.path.join(current_directory, file_name)
data = []
with open(input_file_path, 'r') as mark:
   line = mark.readlines()
data.append(line[0][:-1].split(","))
for i in range(1,len(line)):
   line1 = line[i][:-1].replace('B','01')
   marker = line1.split(",")
   data.append(marker)


taxaname= data[0]
taxa_data = data[1:]
taxanum=len(taxaname)

triplet = []
for j in itertools.combinations(list(range(len(taxa_data[0]))), 3):
	na =0
	nb =0
	nc =0
	for i in range(len(taxa_data)):        
		if (taxa_data[i][j[0]] == '1')&(taxa_data[i][j[1]]== '1')&(taxa_data[i][j[2]] == '0'):
			nc+=1
		if (taxa_data[i][j[0]] == '1')&(taxa_data[i][j[1]] == '0')&(taxa_data[i][j[2]] == '1'):
			nb+=1
		if((taxa_data[i][j[0]] == '0')&(taxa_data[i][j[1]] == '1')&(taxa_data[i][j[2]] == '1')):
			na+=1
	triplet.append([na,nb,nc])
     
#print(triplet)
# buneman weight
weight = triplet
buneman_weight = [[0,0,0] for i in range(len(triplet))]
for i in range(len(triplet)):
    for j in range(3):
        buneman_weight[i][j] = weight[i][j] - sum(weight[i])+max(weight[i])+min(weight[i])
        if buneman_weight[i][j] < 0:
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

#print(cluster)

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



user_input = input("choose threshold: ")
#user_input=3000
#threshold
bunemancluster = []
for i in range(len(indexnum)):
    if indexnum[i]> float(user_input):
        bunemancluster.append(cluster[i][0])


bunemancluster = [[taxaname[index-1] for index in sub_list] for sub_list in bunemancluster]
bunemancluster = sorted(bunemancluster,key = len)


def flatten(nested_list):
    result = []
    for item in nested_list:
        if isinstance(item, list):
            result.extend(flatten(item))  # 递归处理嵌套的列表
        else:
            result.append(item)  # 非列表元素直接添加到结果中
    return result

bunemantree=[]
for i in range(len(bunemancluster)):
	found= False
	newcluster = []
	for j in bunemantree:
		bj = flatten(j)
		if set(bj).issubset(set(bunemancluster[i])):
			bunemantree.remove(j)
			newcluster.append(j)
			found = True
	if found:	
		rest = list(set(bunemancluster[i])-set(flatten(newcluster)))
		bunemantree.append(newcluster+rest)	
	else: 	
		bunemantree.append(bunemancluster[i])
	



fulltree = bunemantree+list(set(taxaname)-set(flatten(bunemantree)))

fulltree = str(fulltree).replace("'","").replace("[","(").replace("]",")").replace(" ","")

print(fulltree)


file_name = 'bunemantree.txt'
output_file_path = os.path.join(current_directory, file_name)
file = open(output_file_path, "w+")
file.write(fulltree+';')
file.close()

