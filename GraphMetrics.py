__author__ = 'Rafael'

import itertools
import re

def read_input(filename):
    f = open(filename,'r')
    #First, we read the node labels...
    l = f.readline()
    pattern = re.compile("{.*}")
    l = pattern.match(l).group()
    l = l[1:-1]
    l = l.split(',') #Now, we have our labels...
    #print(l)
    seq = f.readline()
    if pattern.match(seq) is None:
        print("Error in adjacency list formation!")
        raise(Exception)
    seq = seq[1:-1]
    pattern = re.compile("((\((.*-.*)\),)*?\(.*-.*\))?")
    seq = pattern.match(seq).group()
    v = []
    if seq is not '':
        for k in seq.split(','):
            p = k[1:-1].split('-')
            v.append(p)
	#v now has our adjacency list, so let's build up the resulting adjacency matrix...
    m = [[] for i in l]
    for i in range(0,len(m)):
        m[i] = [0 for i in l]
    for i in v:
        idx_first = l.index(i[0])
        idx_second = l.index(i[1])
        m[idx_first][idx_second] += 1
        m[idx_second][idx_first] += 1
    #print(m)
    #print(l)
    return m,l


"""
Older methods for calculating alpha and omega - brute force over the combinations
"""
#Idea: calculate this for every connected component then aggregate results
#If is_alpha is set to true, calculates the alpha and maximum independent set for the input graph. Otherwise, finds maximum clique.

def calc_param(adj_mtx,is_alpha, name_lst = None):
    rs_lst = None
    rs_siz = len(adj_mtx)

    __found = False
    __own_names = [i for i in range(0, rs_siz)]

    #While the solution still wasn't found...:
    while not __found and rs_siz > 0:
        lst_comb = list(itertools.combinations(__own_names, rs_siz))

        for i in lst_comb:
            t_el = True
            for u in i:
                for v in i[i.index(u) + 1:]:
                    if (is_alpha and adj_mtx[u][v] != 0) or (not is_alpha and adj_mtx[u][v] == 0):
                        t_el = False
                        break
                if t_el == False:
                    break
            if t_el == True:
                rs_lst = list(i)
                __found = True
                break
        if not __found:
            rs_siz -= 1
    if rs_lst is not None:
        if name_lst is not None:
            rs_lst = [name_lst[i] for i in rs_lst]
    else:
        rs_lst = []
    return (rs_lst,rs_siz)

def calc_alpha_old(adj_mtx, name_lst = None):
    return calc_param(adj_mtx, True, name_lst)
	

def calc_omega_old(adj_mtx, name_lst = None):
    return calc_param(adj_mtx, False, name_lst)


"""
Newer methods - using maghout's method for calculating link coverages
"""

def p_m(mat):
    for i in mat:
        print(i)

#Multiplies "or" statements together,
def mul_and_verify(l_1,l_2):

    #Multiplication of two conglomerate elements: union of two sets
    #e.g.: ABCD * BCDE = ABCDE
    def mul(el_1,el_2):
        __res = []
        for k in el_1:
            if k not in __res:
                __res.append(k)
        for k in el_2:
            if k not in __res:
                __res.append(k)
        return __res

    #Verifies whether we can simplify two expressions linked by an "or" operation
    #e.g.: ABCD or BC => BC
    #returns: (is_subset(set_1,set_2))

    def verify_if_subset(set_1,set_2):
        s1 = set_1 #Pointer to the smallest of the sets
        s2 = set_2
        if len(set_1) > len(set_2):
            s1 = set_2
            s2 = set_1
        for i in s1:
            if i not in s2:
                return False, s2
        #We managed to find every element of s1 in s2...
        return True, s2

    #Returns the final, verified list
    def verify_whole_lst(lst):
        #print(lst)
        __res = list(lst)
        for i in range(0,len(lst)):
            for k in range(i+1,len(lst)):
                is_subset,larger_el = verify_if_subset(lst[i],lst[k])
                if is_subset:
                    if larger_el in __res:
                        __res.remove(larger_el)
        return __res

    __final_res = []

    for i in l_1:
        for k in l_2:
            __final_res.append(mul(i,k))
    return verify_whole_lst(__final_res)



#Function that calculates graph link coverages
def calc_maghout(mat,name_lst = None):
    rs_lst = None
    rs_siz = len(mat)
    __own_names = [i for i in range(0, rs_siz)]
    link_lst = []

    for i in range(0,len(mat)-1):
        __scnd_term = []
        for k in range(i+1,len(mat[i])):
            if mat[i][k] is not 0:
                __scnd_term.append(__own_names[k])
        if len(__scnd_term) is not 0:
            link_lst.append(([__own_names[i]], __scnd_term))

    #We should check if the resulting list is empty...
    if len(link_lst) is 0:
        return [[]]
    #Now, we calculate the product of the pairs:
    __final_result = link_lst.pop(0)
    for i in link_lst:
        __final_result = mul_and_verify(__final_result, i)

    if name_lst is not None:
        __translated = []
        for i in __final_result:
            v = []
            for k in i:
                v.append(name_lst[k])
            __translated.append(v)
        return __translated
    else:
        return __final_result

def calc_complimentary_set(name_lst,cvrg_lst):
    __res_lst = []
    for i in cvrg_lst:
        __new_i_set = []
        for k in name_lst:
            if k not in i:
                __new_i_set.append(k)
        __res_lst.append(__new_i_set)
    return __res_lst

def calc_set_sizes(lst_sets):
    return [len(i) for i in lst_sets]

def calc_max_size(lst_sizes):
    max = lst_sizes.pop(0)
    for i in lst_sizes:
        if i > max:
            max = i
    return max

def calc_complimentary_mat(mat):
    res_mat = [[] for i in mat]
    for i in res_mat:
        for k in mat:
            i.append(0)
    for i in range(0,len(mat)):
        for k in range(0,len(mat[i])):
            if i != k:
                if mat[i][k] == 0:
                    res_mat[i][k] = 1
                else:
                    res_mat[i][k] = 0
    return res_mat

def calc_cliques(l_name,mat):
    m_comp = calc_complimentary_mat(mat)
    cov = calc_maghout(m_comp,l_name)
    indep = calc_complimentary_set(l_name,cov)
    return (calc_max_size(calc_set_sizes(indep)),indep)
	
def calc_indep_sets(l_name,mat):
    cov = calc_maghout(mat,l_name)
    indep = calc_complimentary_set(l_name,cov)
    return (calc_max_size(calc_set_sizes(indep)),indep)
	

#Reading the graph input file...
m,l = read_input("input_g2.txt")

#Calculating alpha and the maximal independent sets - that simple!
alpha,i_sets = calc_indep_sets(l,m)

i_sets.sort(key=len,reverse=True)
print("\n\nAlpha: %d\n\nIndependent Sets:\n\t%s"%(alpha,str(i_sets)))
maxim_i = []
for i in i_sets:
    if len(i) == alpha:
        maxim_i.append(i)
print("Largest (maximum) independent set(s):\n\t%s"%maxim_i)

#Calculating omega and the maximal cliques...
omega,cliques = calc_cliques(l,m)
cliques.sort(key=len,reverse=True)
print("\n\nOmega: %d\n\nCliques:\n\t%s"%(omega,str(cliques)))
maxim_c = []
for i in cliques:
    if len(i) == omega:
        maxim_c.append(i)
print("Largest (maximum) clique(s):\n\t%s"%maxim_c)


