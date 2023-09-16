from Bio import Phylo
from itertools import combinations
import re
import os
import pandas as pd

def ProcessOutputHmmscan(input_path,output_path):
    files = os.listdir(input_path)
    for file in files:
        flag = []
        try:
            rm_line = []
            with open (input_path+file,'r') as f:
                for i in f.readlines():
                    flag.append(i)
                    if '#' in i:
                        continue
                    else:
                        rm_line.append(i.strip())
            if '#' not in flag[-1]:
                print(file+" is an incomplete file")  
                continue
            domain = []
            gene = []
            for i in rm_line:
                i = i.split()
                domain.append(str(i[0]))
                gene.append(str(i[2]))
            c = {"domain":domain,"gene":gene}
            data = pd.DataFrame(c)
            try:
                data.to_csv(output_path+file+'_out', header=False, index=False)
            except:
                print(file+" has a fault(on to_csv)")
        except:
            print(file+" has a fault")
            continue
            
def CreateSet(input_path,output_path):
    files = os.listdir(input_path)
    for file in files:
        name = re.search('(.*)\.',file).group(1)
        set1 = set()
        list1 = []
        with open(input_path+file,'r') as f:
            for i in f.readlines():
                i = i.strip().split(',')
                set1.add(i[0])
        for i in set1:
            list1.append(i)
        list1.sort()
        with open(output_path+name+'_set','w') as f:
            f.write(name+'\n')
            for domain in list1:
                f.write(domain+'\n')
                
def CreateCombineSet(input_path,output_path):
    files = os.listdir(input_path)
    for file in files:
        name = re.search('(.*)\.',file).group(1)
        set1 = set()
        list1 = []
        with open(input_path+file,'r') as f:
            combination = {}
            for i in f.readlines():
                domain = i.strip().split(',')[0]
                gene = i.strip().split(',')[1]
                try:
                    combination[gene].append(domain)
                except:
                    combination[gene] = []
                    combination[gene].append(domain)
        for a,b in combination.items(): 
            c= ''
            b.sort()
            for domain in b:
                c = c +domain+' '
            set1.add(c.strip())
        for i in set1:
            list1.append(i)
        list1.sort()
        with open(output_path+name+'_combine_set','w') as f:
            f.write(name+'\n')
            for domain in list1:
                f.write(domain+'\n')
                
def IntTree(file,format1):
    AllLeafs = {}
    AllNodes = {}
    Parents = {}
    Childs = {}
    tree = Phylo.read(file,format1)
    Phylo.draw(tree)
    for i in tree.get_terminals():
        AllLeafs[str(i)] = i
        
    for i in tree.get_nonterminals():
        AllNodes[str(i)] = i
        
    for clade in tree.find_clades(order="level"):
        for child in clade:
            Parents[str(child)] = clade
            
    for node_str,node_clade in AllNodes.items():
        for leaf_str in AllLeafs.keys():
            if node_clade.is_parent_of(leaf_str):
                try:
                    Childs[node_str].add(leaf_str)
                except:
                    Childs[node_str] = set()
                    Childs[node_str].add(leaf_str)
    return(AllLeafs,AllNodes,Parents,Childs,tree)

def AncestralContent(species,tree):
    AncestralNode = {}
    for node in tree.get_nonterminals():
        AncestralNode[str(node)] = set()    
    LeafCombination = list(combinations(list(species.keys()),2))
    for i in LeafCombination:
        Trace = tree.trace(i[0],i[1])
        if str(Trace.pop()).startswith('M'):
            pass
        else:
            print('pop failed')
            break     
        A_content = species[i[0]]&species[i[1]]
        for A in Trace:
            AncestralNode[str(A)] = AncestralNode[str(A)]|A_content
    return AncestralNode

def FormLineageSet(dataset):
    all_domain = set()
    lineage = {}
    files = os.listdir(dataset)
    for file in files:
        lineage[file] = set()
        with open(dataset+file,'r') as f:
            for i in f.readlines():
                lineage[file].add(i.strip())
                all_domain.add(i.strip())   
    return (lineage,all_domain)
                
def DmtbloutInit(dmtblout_path,name=0):
    with open(dmtblout_path,'r') as f:
        flag = []
        rm_line = []
        for i in f.readlines():
            flag.append(i)
            if '#' in i:
                continue
            else:
                rm_line.append(i.strip())
        if '#' not in flag[-1]:
            print(name+" is an incomplete file")  
    combination = {}
    for i in rm_line[0:30]:
        gene =  i.split()[3]
        domain = i.split()[0]
        position1 = i.split()[19]
        position2 = i.split()[20]
        score = i.split()[21]
        length = int(position2)-int(position1)
        try:
            combination[gene][domain].append((position1,position2,score))
        except:
            try:
                combination[gene][domain] = []
                combination[gene][domain].append((position1,position2,score))
            except:
                combination[gene] = {}
                combination[gene][domain] = []
                combination[gene][domain].append((position1,position2,score))
    return combination            
                
def LineageSet(file_lists,path):
    set_list = []
    for i in file_lists:
        locals()[str(i)] = set()
        with open(path+i,'r') as f:
            num = 0
            for a in f.readlines():
                if num ==0:
                    num = num+1
                    continue
                if num > 0:
                    a = a.strip()
                    locals()[str(i)].add(a)
        set_list.append(locals()[str(i)])
    set1 = set_list[0]
    for b in set_list:
        set1 = set1.union(b)
    return set1

def EverySet(file_lists,path):
    set_list = []
    for i in file_lists:
        locals()[str(i)] = set()
        with open(path+i,'r') as f:
            num = 0
            for a in f.readlines():
                if num ==0:
                    num = num+1
                    continue
                if num > 0:
                    a = a.strip()
                    locals()[str(i)].add(a)
        set_list.append(locals()[str(i)])
    return set_list

def CoreSet(set1,everyset):
    core = set()
    for domain in set1:
        num = 0
        for s in everyset:
            if domain in s:
                num = num+1
        if num/len(everyset)>=0.8:
            core.add(domain)
    return core

if __name__ == '__main__':
    hmmscan_files = os.listdir('./examples/input/hmmscan')
    for i in hmmscan_files:
        os.mkdir('./examples/input/hmmscan_output/'+i)
        ProcessOutputHmmscan('./examples/input/hmmscan/'+i+'/','./examples/input/hmmscan_output/'+i+'/')
        os.mkdir('./examples/input/domain_set/'+i)
        CreateSet('./examples/input/hmmscan_output/'+i+'/','./examples/input/domain_set/'+i+'/')
        lineageset = LineageSet(os.listdir('./examples/input/domain_set/'+i+'/'),'./examples/input/domain_set/'+i+'/')
        everyset = EverySet(os.listdir('./examples/input/domain_set/'+i+'/'),'./examples/input/domain_set/'+i+'/')
        coreset = CoreSet(lineageset,everyset)
        with open('./examples/input/domain_coreset/'+i,'w') as f:
            for i in coreset:
                f.write(i+'\n')
        
    lineage = FormLineageSet('./examples/input/domain_coreset/')[0]
    returns = IntTree("Species.nwk","newick")
    
    AllLeafs = returns[0]
    AllNodes = returns[1]
    Parents = returns[2]
    Childs = returns[3]
    tree = returns[-1]

    AncestralNode = AncestralContent(lineage,tree)
    AncestralNode_gain_loss = {}
    for a,b in AncestralNode.items():
        if a =='A5':
            with open('./examples/output/'+a,'w') as f:
                f.write(a+'\n')
                f.write('Ancestral: '+str(len(b))+'\n')
            continue
        loss = AncestralNode[str(Parents[a])] - b
        gain = b - AncestralNode[str(Parents[a])]
        with open('./examples/output/'+a,'w') as f:
            f.write(a+'\n')
            f.write('Ancestral: '+str(len(b))+'\n')
            f.write('gain: '+str(len(gain))+'\n')
            f.write('loss: '+str(len(loss))+'\n')
