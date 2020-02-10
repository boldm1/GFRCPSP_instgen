import numpy as np
import os
import random
import math
import datetime

random.seed(1996)

#fixed parameters
min_sources = 1
max_sources = 10
min_sinks = 1
max_sinks = 10
max_pred = 3 #max number of non-redundant predecessors per task
max_succ = 3 #max number of non-redundant successors per task
min_ncs = 1 #min number of cycle structures
max_ncs = 5
min_cs = 4 #min size of a cycle structure (counting event nodes)
max_cs = 20
min_back = 0.05 #min proportion of max time-lags
max_back = 0.2
epsilon = 2 #max relative deviation of min time-lag from d_max
rho = 0.02 #degree of redundancy
delta = 0.5 #cycle densification parameter
tau = 0.3 #cycle structure tightness
RF = 0.75 #resource factor

n_res = 5
min_qmin = 1
max_qmin = 3
min_qmax = 3
max_qmax = 6
min_w = 12 #should be at least l_min*max_qmin
max_w = 40

#variable parameters
n = 30
OS = [0.4] #order strength
RS = [0.05,0.15,0.25] #resource strength

class Task():
    def __init__(self, id):
        self.id = id
        self.q_min = [0 for i in range(n_res)]
        self.q_max = [0 for i in range(n_res)]
        self.w = 0
        self.d_min = 0
        self.d_max = 0
        self.SS_succ = []
        self.SF_succ = []
        self.FS_succ = []
        self.FF_succ = []
        self.ES = 0
        self.LS = 0
        self.EF = 0
        self.LF = 0

class Node(): #node in project network corresponding to start or end of task
    def __init__(self, num):
        self.num = num
        self.succ = [] #list of forward successors
        self.pred = [] #list of forward predecessors

try:
    os.mkdir('p{}'.format(n))
except FileExistsError:
    pass
os.chdir('p{}'.format(n))
p_file = open('params.txt', 'w+')
p_file.write('{}\n\nProblem type: GFRCPSP\n----------------------------------\nPrecedence network parameters:\n----------------------------------\nNumber of activities: {}\nMin number of source activities: {}\nMax number of source activities: {}\nMin number of sink activities: {}\nMax number of sink activities: {}\nMax number of non-redundant predecessor activities: {}\nMax number of non-redundant successors activities: {}\nOrder strength: {}\nDegree of redundancy: {} \nMax relative deviation of max time-lags from max duration: {}\nMin proportion of max time-lags: {}\nMax proportion of max time-lags: {}\nMin number of cycle structures: {}\nMax number of cycle structures: {}\nMin number of nodes per cycle structure: {}\nMax number of nodes per cycle structure: {}\nCycle structure density parameter: {}\nCycle structure tightness parameter: {}\n----------------------------------\nResource parameters:\n----------------------------------\nNumber of resources: {}\nResource factor (proportion of non-zero resource requirements): {}\nMin resource allocation LB: {}\nMax resource allocation LB: {}\nMin resource allocation UB: {}\nMax resource allocation UB: {}\nMin principle resource requirement: {}\nMax principle resource requirement: {}\nResource strength: {}'.format(datetime.datetime.now().replace(microsecond=0), n, min_sources, max_sources, min_sinks, max_sinks, max_pred, max_succ, OS, rho, epsilon, min_back, max_back, min_ncs, max_ncs, min_cs, max_cs, delta, tau, n_res, RF, min_qmin, max_qmin, min_qmax, max_qmax, min_w, max_w, RS))
p_file.close()
stat_file = open('stat.txt', 'w+')
stat_file.write('Name:\tTime of generation: \t Number of activities: \t Number of source activities:\tNumber of sink activities:\tMean number of non-redundant predecessor nodes:\t Number of non-redundant arcs:\tOrder strength: \t  Number of redundant arcs:\t Number of max time-lags: \t Number of cycle structures:\t Mean number of nodes per cycle structure: \tMin number of nodes in a cycle structure:\tMax number of nodes in a cycle structure:\tNumber of resources \tMean principle resource requirement: \tMean minimal activity duration: \t Mean maximal activity duration: \t Resource factor: \t Resource strength: \t Network-based LB: \t UB on optimal makespan:\n')

inst = 0 #instance number
for o in OS:
    for rs in RS:
        for rep in range(30):
            
            tasks = {i:Task(i) for i in range(1,n+1)}
            V = {i:Node(i) for i in range(2*n+2)}
            N = [i for i in range(1,2*n+1)] #node numbers
            
            inst += 1
            print('###############################################################################')
            print('generating instance', inst)
            
            #dict for containing instance stats
            stats = {'name':'psp{}'.format(inst), 'n_act':n, 'n_res':n_res, 'RS':rs}

            l_min = random.randint(2,4)

            ### Resource requirements and bounds
            #principle resource
            for i in tasks:
                tasks[i].w = random.randint(min_w, max_w)
                tasks[i].q_min[0] = random.randint(min_qmin, max_qmin)
                tasks[i].q_max[0] = min(random.randint(max(tasks[i].q_min[0], min_qmax), max_qmax), math.floor(tasks[i].w/l_min))
                tasks[i].d_min = math.ceil(tasks[i].w/tasks[i].q_max[0])
                tasks[i].d_max = math.ceil(tasks[i].w/tasks[i].q_min[0])
            #non-principle resources
            rf = sum(tasks[i].q_min[r] != 0 for i in tasks for r in range(n_res))/(n_res*n)
            while rf < RF:
                i = random.choice([i for i in tasks if len([r for r in range(n_res) if tasks[i].q_min[r] == 0]) > 0])
                r = random.choice([r for r in range(n_res) if tasks[i].q_min[r] == 0])
                tasks[i].q_min[r] = random.randint(min_qmin, max_qmin)
                tasks[i].q_max[r] = max(tasks[i].q_min[0], random.randint(min_qmax, max_qmax))
                rf = sum(tasks[i].q_min[r] != 0 for i in tasks for r in range(n_res))/(n_res*n)

            stats['rf'] = rf
            avg_w = sum([tasks[i].w for i in tasks])/n
            stats['avg_w'] = avg_w
            avg_dmin = sum([tasks[i].d_min for i in tasks])/n
            stats['avg_dmin'] = avg_dmin
            avg_dmax = sum([tasks[i].d_max for i in tasks])/n
            stats['avg_dmax'] = avg_dmax
                
            ###Project network
            G = {i:[] for i in range(2*n+2)}
            A = np.zeros((2*n+2, 2*n+2)) #adjacency matrix
            R = np.identity(2*n+2) #reachability matrix
            R_forw = np.identity(2*n+2) #reachability matrix for forward arcs only. (A full reachability matrix with backwards max duration arcs is required for preventing infeasible network. This reachability matrix is solely used for computing order strength)
            R2 = np.identity(2*n+2) #squared reachability matrix - # nodes in all paths connecting i to j

            n_sources = random.randint(min_sources, max_sources)
            n_sinks = random.randint(min_sinks, max_sinks)
            sources = [i for i in range(1,n_sources+1)]
            sinks = [i for i in range(n-n_sinks+1,n+1)]
#            print('sources', sources)
#            print('sinks', sinks)

            stats['n_sources'] = n_sources
            stats['n_sinks'] =n_sinks
            
            forw = [] #min time-lag arcs
            back = [] #max time-lag arcs
            SS = []
            SF = []
            FS = []
            FF = []

            #parameters to track stats
            n_nred = 0 #number of redundant arcs
            n_red = 0 #number of redundant arcs
            n_mtl = 0 #number of max time-lags
            
            #adding min duration arcs
            for i in tasks:
                forw.append((2*i-1,2*i,tasks[i].d_min))
                SF.append((2*i-1,2*i,tasks[i].d_min))
                V[2*i-1].succ.append((2*i,tasks[i].d_min))
                V[2*i].pred.append((2*i-1,tasks[i].d_min))
                G[2*i-1].append(2*i)
                A[2*i-1,2*i] = 1
                R[2*i-1,2*i] = 1
                R_forw[2*i-1,2*i] = 1
                R2[2*i-1,2*i] = 2
            #adding max duration arcs
            for i in tasks:
                back.append((2*i,2*i-1,-tasks[i].d_max))
                FS.append((2*i,2*i-1,-tasks[i].d_max))
                G[2*i].append(2*i-1)
                A[2*i,2*i-1] = 1
                R[2*i,2*i-1] = 1
                R2[2*i,2*i-1] = 2
                R2[2*i-1,2*i-1] = 2
                R2[2*i,2*i] = 2

            source_nodes = [i for i in range(1,2*n_sources+1)]
            sink_nodes = [i for i in range(2*n-2*n_sinks+1,2*n+1)]
            
            #adding predecessors for every activity
            print('adding predecessors')
            non_sources = [i for i in tasks if i not in sources]
            while non_sources != []:
                j = random.choice(non_sources)
                non_sources.remove(j)
                j = random.choice([2*j-1,2*j])
                P = [i for i in N if i not in sink_nodes if R[i,j] + sum(R[g,i]*R[g,j] for g in N if R[g,i]==1) == 0 if R[j,i]==0]
#                P = [i for i in N if i not in sink_nodes if R[i,j]+sum(R[e[0],i]*R[j,e[1]] for e in forw)==0 if R[j,i]==0 if len(V[i].succ)<max_succ]
                if P == []:
                    continue
                else:
                    i = random.choice(P)
                    #create arc i->j
                    n_nred += 1
                    weight = random.randint(0,math.ceil(epsilon*tasks[math.ceil(i/2)].d_min))
                    print((i,j), weight)
                    V[i].succ.append((j,weight))
                    V[j].pred.append((i,weight))
                    forw.append((i,j,weight))
                    if i%2 == 1:
                        if j%2 == 1:
                            SS.append((i,j,weight))
                            tasks[math.ceil(i/2)].SS_succ.append((math.ceil(j/2),weight)) 
                        else:
                            SF.append((i,j,weight))
                            tasks[math.ceil(i/2)].SF_succ.append((math.ceil(j/2),weight)) 
                    else:
                        if j%2 == 1:
                            FS.append((i,j,weight))
                            tasks[math.ceil(i/2)].FS_succ.append((math.ceil(j/2),weight)) 
                        else:
                            FF.append((i,j,weight))
                            tasks[math.ceil(i/2)].FF_succ.append((math.ceil(j/2),weight)) 
                    G[i].append(j)
                    A[i,j] = 1
                    R[i,j] = 1
                    R_forw[i,j] = 1
                    for g in V:
                        if R[g,i] == 1:
                            for h in V:
                                if R[j,h] == 1:
                                    R[g,h] = 1
                                    R_forw[g,h] = 1
                    for g in V:
                        if R[g,i] == 1:
                            for h in V:
                                if R[j,h] == 1:
                                    R2[g,h] = sum([R[g,x]*R[x,h] for x in V])

            #adding successors for every activity
            print('adding successors')
            non_sinks = [i for i in tasks if i not in sinks]
            while non_sinks != []:
                i = random.choice(non_sinks)
                non_sinks.remove(i)
                i = random.choice([2*i-1,2*i])
                S = [j for j in N if j not in source_nodes if R[i,j] + sum(R[g,i]*R[g,j] for g in N if R[g,i]==1) == 0 if R[j,i]==0 if len(V[i].pred)<max_pred]
#                S = [j for j in N if j not in source_nodes if R[i,j]+sum(R[e[0],i]*R[j,e[1]] for e in forw)==0 if R[j,i]==0 if len(V[j].pred)<max_pred]
                if S == []:
                    continue
                else:
                    j = random.choice(S)
                    #create arc i->j
                    n_nred += 1
                    weight = random.randint(0,math.ceil(epsilon*tasks[math.ceil(i/2)].d_min))
                    print((i,j), weight)
                    V[i].succ.append((j,weight))
                    V[j].pred.append((i,weight))
                    forw.append((i,j,weight))
                    if i%2 == 1:
                        if j%2 == 1:
                            SS.append((i,j,weight))
                            tasks[math.ceil(i/2)].SS_succ.append((math.ceil(j/2),weight)) 
                        else:
                            SF.append((i,j,weight))
                            tasks[math.ceil(i/2)].SF_succ.append((math.ceil(j/2),weight)) 
                    else:
                        if j%2 == 1:
                            FS.append((i,j,weight))
                            tasks[math.ceil(i/2)].FS_succ.append((math.ceil(j/2),weight)) 
                        else:
                            FF.append((i,j,weight))
                            tasks[math.ceil(i/2)].FF_succ.append((math.ceil(j/2),weight)) 
                    G[i].append(j)
                    A[i,j] = 1
                    R[i,j] = 1
                    R_forw[i,j] = 1
                    for g in V:
                        if R[g,i] == 1:
                            for h in V:
                                if R[j,h] == 1:
                                    R[g,h] = 1
                                    R_forw[g,h] = 1
                    for g in V:
                        if R[g,i] == 1:
                            for h in V:
                                if R[j,h] == 1:
                                    R2[g,h] = sum([R[g,x]*R[x,h] for x in V])

            #adding non-redundant arcs until OS is satisfied
            print('Adding non-redundant arcs')
            P = [i for i in N if i not in sink_nodes if len(V[i].succ) < max_succ]
            current_os = (sum(R_forw[i,j] for i in N for j in N) - 3*n)/(2*n*(n+1))
            print('current_os', current_os)
            print('OS', o)
            while current_os < o: #o - requires OS value
                if P == []:
                    break
                else:
                    i = random.choice(P)
                    S = [j for j in N if j not in source_nodes if R[i,j] + sum(R[g,i]*R[g,j] for g in N if R[g,i]==1) == 0 if R[j,i]==0 if len(V[i].pred)<max_pred]
#                    S = [j for j in N if j not in source_nodes if R[i,j]+sum(R[e[0],i]*R[j,e[1]] for e in forw)==0 if R[j,i]==0 if len(V[j].pred)<max_pred]
                    if S == []:
                        P.remove(i)
                    else:
                        j = random.choice(S)
                        #create arc i->j
                        n_nred += 1
                        weight = random.randint(0,math.ceil(epsilon*tasks[math.ceil(i/2)].d_min))
                        print((i,j), weight)
                        V[i].succ.append((j,weight))
                        if len(V[i].succ) == max_succ:
                            P.remove(i)
                        V[j].pred.append((i,weight))
                        forw.append((i,j,weight))
                        if i%2 == 1:
                            if j%2 == 1:
                                SS.append((i,j,weight))
                                tasks[math.ceil(i/2)].SS_succ.append((math.ceil(j/2),weight)) 
                            else:
                                SF.append((i,j,weight))
                                tasks[math.ceil(i/2)].SF_succ.append((math.ceil(j/2),weight)) 
                        else:
                            if j%2 == 1:
                                FS.append((i,j,weight))
                                tasks[math.ceil(i/2)].FS_succ.append((math.ceil(j/2),weight)) 
                            else:
                                FF.append((i,j,weight))
                                tasks[math.ceil(i/2)].FF_succ.append((math.ceil(j/2),weight)) 
                        G[i].append(j)
                        A[i,j] = 1
                        R[i,j] = 1
                        R_forw[i,j] = 1
                        for g in V:
                            if R[g,i] == 1:
                                for h in V:
                                    if R[j,h] == 1:
                                        R[g,h] = 1
                                        R_forw[g,h] = 1
                        for g in V:
                            if R[g,i] == 1:
                                for h in V:
                                    if R[j,h] == 1:
                                        R2[g,h] = sum([R[g,x]*R[x,h] for x in V])
                        current_os = (sum(R_forw[i,j] for i in N for j in N) - 3*n)/(2*n*(n+1))
                        print('current_os', current_os)
            stats['n_nred'] = n_nred 
            stats['OS'] = current_os

            #adding redundant arcs
            print('adding redundant arcs')
            max_red = sum(R_forw[i,j] for i in N for j in N) - 2*n - len(forw)
            n_red = math.floor(rho*max_red)
            P = [i for i in N if i not in sink_nodes if len(V[i].succ) < 2*n - len(source_nodes)]
            for iteration in range(n_red):
                if P == []:
                    break
                else:
                    i = random.choice(P)
                    S = [j for j in N if j not in source_nodes if R[i,j] == 1 if A[i,j] == 0 if R[j,i] == 0]
                    if S == []:
                        P.remove(i)
                    else:
                        j = random.choice(S)
                        #create arc i->j
                        n_red += 1
                        weight = random.randint(0,math.ceil(epsilon*tasks[math.ceil(i/2)].d_min))
                        print((i,j), weight)
                        V[i].succ.append((j,weight))
                        if len(V[i].succ) == 2*n - len(source_nodes):
                            P.remove(i)
                        V[j].pred.append((i,weight))
                        forw.append((i,j,weight))
                        if i%2 == 1:
                            if j%2 == 1:
                                SS.append((i,j,weight))
                                tasks[math.ceil(i/2)].SS_succ.append((math.ceil(j/2),weight)) 
                            else:
                                SF.append((i,j,weight))
                                tasks[math.ceil(i/2)].SF_succ.append((math.ceil(j/2),weight)) 
                        else:
                            if j%2 == 1:
                                FS.append((i,j,weight))
                                tasks[math.ceil(i/2)].FS_succ.append((math.ceil(j/2),weight)) 
                            else:
                                FF.append((i,j,weight))
                                tasks[math.ceil(i/2)].FF_succ.append((math.ceil(j/2),weight)) 
                        G[i].append(j)
                        A[i,j] = 1

            stats['n_red'] = n_red
            
            #computing distance matrix 
            dgraph = np.full((2*n+2,2*n+2),-10000)
            for i in range(2*n+2):
                dgraph[i,i] = 0
            for arc in forw + back:
                dgraph[arc[0],arc[1]] = arc[2]
#            print(dgraph)
            #floyd-warshall algorithm
            for k in range(2*n+2):
                for i in range(2*n+2):
                    for j in range(2*n+2):
                        dgraph[i,j] = max(dgraph[i,j], dgraph[i,k]+dgraph[k,j])

            #Creating new cycle structures
            print('creating new cycles structures')
            n_back = random.randint(math.ceil(len(forw)*min_back), math.floor(len(forw)*max_back))
            n_cs = random.randint(min_ncs, min(n_back, max_ncs)) #number of cycle structures
            P = [i for i in N]
            n_cycles = 0
            cycles = {}
            while n_cycles < n_cs:
                if P == []:
                    break
                else:
                    i = random.choice(P)
                    S = [j for j in N if R[i,j] == 1 if R[j,i] == 0 if sum([R[i,h]*R[h,j] for h in N if R2[h,h] > 2]) == 0 if min_cs <= R2[i,j] if R2[i,j] <= max_cs]
                    if S == []:
                        P.remove(i)
                    else:
                        j = random.choice(S)
                        #create arc j->i
#                        print('i,j', i,j)
                        min_weight = dgraph[i,j] #min max time-lag from i to j
                        #calculating max weight
                        max_weight = 0
                        for g in [g for g in N if g != j if R[i,g]*R[g,j] == 1]:
                            if g%2 == 1:
                                max_SS = max([arc[2] for arc in SS if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_SF = max([max(tasks[math.ceil(g/2)].d_max+tasks[math.ceil(arc[1]/2)].d_max, arc[2]) for arc in SF if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_weight += max(max_SS, max_SF)
                            else:
                                max_FS = max([arc[2] for arc in FS if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_FF = max([max(tasks[math.ceil(arc[1]/2)].d_max, arc[2]) for arc in FF if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_weight += max(max_FS, max_FF)
                        weight = random.randint(math.floor(max_weight-2*(max_weight-min_weight)*tau+(max_weight-min_weight)*tau*tau),math.ceil(max_weight-(max_weight-min_weight)*tau*tau))
                        weight = -weight
                        print((j,i), weight)
                        n_mtl += 1
                        if weight >= 0:
                            V[j].succ.append((i,weight))
                            V[i].pred.append((j,weight))
                            forw.append((j,i,weight))
                        else:
                            back.append((j,i,weight))

                        if j%2 == 1:
                            if i%2 == 1:
                                SS.append((j,i,weight))
                                tasks[math.ceil(j/2)].SS_succ.append((math.ceil(i/2),weight)) 
                            else:
                                SF.append((j,i,weight))
                                tasks[math.ceil(j/2)].SF_succ.append((math.ceil(i/2),weight)) 
                        else:
                            if i%2 == 1:
                                FS.append((j,i,weight))
                                tasks[math.ceil(j/2)].FS_succ.append((math.ceil(i/2),weight)) 
                            else:
                                FF.append((j,i,weight))
                                tasks[math.ceil(j/2)].FF_succ.append((math.ceil(i/2),weight)) 
                        G[j].append(i)
                        for x in P:
                            if R[x,j] == 1 and R[i,x] == 1:
                                P.remove(x) #in this case, x is part of the cycle structure containing i and j
                        A[j,i] = 1
                        R[j,i] = 1
                        for h in N:
                            if R[h,j] == 1:
                                for g in N:
                                    if R[i,g] == 1:
                                        R[h,g] = 1
                        for h in N:
                            if R[h,j] == 1:
                                for g in N:
                                    if R[i,g] == 1:
                                        R2[h,g] = sum([R[h,x]*R[x,g] for x in V])
                        #updating dgraph
                        dgraph[j,i] = weight
                        for l in range(2*n+2):
                            for m in range(2*n+2):
                                dgraph[l,m] = max(dgraph[l,m], dgraph[l,j]+dgraph[j,i]+dgraph[i,m])
                        n_cycles += 1
                        n_back -= 1
                        cycles[n_cycles] = i #a representative activity from cycle from which full cycle is recovered later


            stats['n_cycles'] = n_cycles

            #Extending cycle structures
            print('extending cycle structures')
            n_ext = math.floor((1-delta)*n_back)
            P = [i for i in N]
            new_arcs = 0
            while new_arcs < n_ext:
                if P == []:
                    break
                else:
                    i = random.choice(P)
                    S = [j for j in N if R[i,j] == 1 if R[j,i] == 0 if sum([R[i,h]*R[h,j] for h in N if R2[h,h] > 2]) >= 1 if min_cs <= R2[i,j] if R2[i,j] <= max_cs]
#                    S = [j for j in N if R[i,j] == 1 if R[j,i] == 0 if min_cs <= R2[i,j] if R2[i,j] <= max_cs]
                    if S == []:
                        P.remove(i)
                    else:
                        j = random.choice(S)
                        #create arc j->i
                        n_mtl += 1
                        new_arcs += 1
                        min_weight = dgraph[i,j] #min weight of mtl from i to j
                        #calculating max weight
                        max_weight = 0
                        for g in [g for g in N if g != j if R[i,g]*R[g,j] == 1]:
                            if g%2 == 1:
                                max_SS = max([arc[2] for arc in SS if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_SF = max([max(tasks[math.ceil(g/2)].d_max+tasks[math.ceil(arc[1]/2)].d_max, arc[2]) for arc in SF if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_weight += max(max_SS, max_SF)
                            else:
                                max_FS = max([arc[2] for arc in FS if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_FF = max([max(tasks[math.ceil(arc[1]/2)].d_max, arc[2]) for arc in FF if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_weight += max(max_FS, max_FF)
                        weight = random.randint(math.floor(max_weight-2*(max_weight-min_weight)*tau+(max_weight-min_weight)*tau*tau),math.ceil(max_weight-(max_weight-min_weight)*tau*tau))
                        weight = -weight
                        print((j,i), weight)
                        if weight >= 0:
                            V[j].succ.append((i,weight))
                            V[i].pred.append((j,weight))
                            forw.append((j,i,weight))
                        else:
                            back.append((j,i,weight))
                        if j%2 == 1:
                            if i%2 == 1:
                                SS.append((j,i,weight))
                                tasks[math.ceil(j/2)].SS_succ.append((math.ceil(i/2),weight)) 
                            else:
                                SF.append((j,i,weight))
                                tasks[math.ceil(j/2)].SF_succ.append((math.ceil(i/2),weight)) 
                        else:
                            if i%2 == 1:
                                FS.append((j,i,weight))
                                tasks[math.ceil(j/2)].FS_succ.append((math.ceil(i/2),weight)) 
                            else:
                                FF.append((j,i,weight))
                                tasks[math.ceil(j/2)].FF_succ.append((math.ceil(i/2),weight)) 
                        G[j].append(i)
                        for x in P:
                            if R[i,x] == 1 and R[x,j] == 1:
                                P.remove(x) #in this case, x is part of cycle structure containing i and j
                        A[j,i] = 1
                        R[j,i] = 1
                        for h in N:
                            if R[h,j] == 1:
                                for g in N:
                                    if R[i,g] == 1:
                                        R[h,g] = 1
                        for h in N:
                            if R[h,j] == 1:
                                for g in N:
                                    if R[i,g] == 1:
                                        R2[h,g] = sum([R[h,x]*R[x,g] for x in V])
                        #updating dgraph with floyd-warshall algorithm
                        dgraph[j,i] = weight
                        for l in range(2*n+2):
                            for m in range(2*n+2):
                                dgraph[l,m] = max(dgraph[l,m], dgraph[l,j]+dgraph[j,i]+dgraph[i,m])
                        n_back -= 1
            
            #Densification of cycle structures
            print('Densification of cycle structures')
            P = [i for i in N]
            while n_back > 0:
                if P == []:
                    break
                else:
                    i = random.choice(P)
                    S = [j for j in N if j != i if A[j,i] == 0 if R[i,j] == 1 if R[j,i] == 1]
                    if S == []:
                        P.remove(i)
                    else:
                        j = random.choice(S)
                        #create arc j->i
                        min_weight = dgraph[i,j] #min max time-lag from i to j
                        #calculating max weight
                        max_weight = 0
                        for g in [g for g in N if g != j if R[i,g]*R[g,j] == 1]:
                            if g%2 == 1:
                                max_SS = max([arc[2] for arc in SS if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_SF = max([max(tasks[math.ceil(g/2)].d_max+tasks[math.ceil(arc[1]/2)].d_max, arc[2]) for arc in SF if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_weight += max(max_SS, max_SF)
                            else:
                                max_FS = max([arc[2] for arc in FS if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_FF = max([max(tasks[math.ceil(arc[1]/2)].d_max, arc[2]) for arc in FF if g==arc[0] and R[arc[1],j]==1] + [0])
                                max_weight += max(max_FS, max_FF)
#                        print('min_weight', min_weight)
#                        print('max_weight', max_weight)
                        weight = random.randint(math.floor(max_weight-2*(max_weight-min_weight)*tau+(max_weight-min_weight)*tau*tau),math.ceil(max_weight-(max_weight-min_weight)*tau*tau))
                        weight = -weight
                        print((j,i), weight)
                        if weight >= 0:
                            V[j].succ.append((i,weight))
                            V[i].pred.append((j,weight))
                            forw.append((j,i,weight))
                        else:
                            back.append((j,i,weight))
                            n_mtl += 1
                        if j%2 == 1:
                            if i%2 == 1:
                                SS.append((j,i,weight))
                                tasks[math.ceil(j/2)].SS_succ.append((math.ceil(i/2),weight)) 
                            else:
                                SF.append((j,i,weight))
                                tasks[math.ceil(j/2)].SF_succ.append((math.ceil(i/2),weight)) 
                        else:
                            if i%2 == 1:
                                FS.append((j,i,weight))
                                tasks[math.ceil(j/2)].FS_succ.append((math.ceil(i/2),weight)) 
                            else:
                                FF.append((j,i,weight))
                                tasks[math.ceil(j/2)].FF_succ.append((math.ceil(i/2),weight)) 
                        G[j].append(i)
                        A[j,i] = 1
                        #updating dgraph with floyd-warshall algorithm
                        dgraph[j,i] = weight
                        for l in range(2*n+2):
                            for m in range(2*n+2):
                                dgraph[l,m] = max(dgraph[l,m], dgraph[l,j]+dgraph[j,i]+dgraph[i,m])
                        n_back -= 1
            
            for c in cycles:
                cycles[c] = [j for j in N if R[cycles[c],j] == 1 if R[j,cycles[c]] == 1]
            stats['min_cs'] = min([len(cycles[c]) for c in cycles])
            stats['max_cs'] = max([len(cycles[c]) for c in cycles])
            stats['avg_cs'] = sum([len(cycles[c]) for c in cycles])/len(cycles)
            stats['n_mtl'] = n_mtl
            stats['avg_pred'] = sum([len(V[i].pred) for i in N])/(2*n)

            #adding supersource and supersink arcs
            without_pred = [i for i in tasks if V[2*i-1].pred == []]
            for i in without_pred:
                forw.append((0,2*i-1,0)) #(i,j,weight)
                FS.append((0,i,0))
                V[0].succ.append((2*i-1,0))
                V[2*i-1].pred.append((0,0))
                G[0].append(2*i-1)
                A[0,2*i-1] = 1
                R[0,2*i-1] = 1
                R2[0,2*i-1] = 2
                for j in V:
                    dgraph[0,j] = max(dgraph[0,j], dgraph[2*i-1,j])
            without_succ = [i for i in tasks if V[2*i].succ == []]
            for i in without_succ:
                forw.append((2*i,2*n+1,0))
                FS.append((i,n+1,0))
                V[2*i].succ.append((2*n+1,0))
                V[2*n+1].pred.append((2*i,0))
                G[2*i].append(2*n+1)
                A[2*i,2*n+1] = 1
                R[2*i,2*n+1] = 1
                R2[2*i,2*n+1] = 2
                tasks[i].FS_succ.append((n+1,0))
                for j in V:
                    dgraph[j,2*n+1] = max(dgraph[j,2*n+1], dgraph[j,2*i])

           
            #UB on min project duration
            T_max = 0
            for i in tasks:
                T_max += max([tasks[i].d_max, max([arc[2] for arc in SS if 2*i-1 == arc[0]] + [0]), max([arc[2]-tasks[math.ceil(arc[1]/2)].d_min for arc in SF if 2*i-1 == arc[0]] + [0]), max([tasks[i].d_max+arc[2] for arc in FS if 2*i == arc[0]] + [0]), max([tasks[i].d_max+arc[2]-tasks[math.ceil(arc[1]/2)].d_min for arc in FF if 2*i == arc[0]] + [0])])
            print('T_max', T_max)
            
            for i in tasks:
                tasks[i].ES = dgraph[0][2*i-1]
                tasks[i].LS = T_max-dgraph[2*i-1][2*n+1]
                tasks[i].EF = dgraph[0][2*i]
                tasks[i].LF = T_max-dgraph[2*i][2*n+1]
                print('task {} EF'.format(i), tasks[i].EF)
    
            stats['net_LB'] = max([tasks[i].EF for i in tasks]) 
            stats['T_max'] = T_max
        
            #Resource availability
            min_Rmax = [0 for r in range(n_res)]
            for r in range(n_res):
                min_Rmax[r] = max(tasks[i].q_max[r] for i in tasks)
            max_Rmax = [0 for r in range(n_res)]
            for r in range(n_res):
                max_Rmax[r] = max([sum(tasks[i].q_max[r] for i in tasks if tasks[i].ES <= t if t <= tasks[i].LF) for t in range(T_max)])
            
            R_max = [0 for r in range(n_res)]
            for r in range(n_res):
                R_max[r] = min_Rmax[r] + math.ceil(rs*(max_Rmax[r]-min_Rmax[r]))


            
            ### write folder of test instances
            f = open('psp{}.sch'.format(inst), 'w+')
            raw_line1 = '{}\t{}\n'.format(n, n_res)
            f.write(raw_line1)
            b1_line = ['0','1','0','0',['{}'.format(len(V[0].succ))],'0', ['{}'.format(math.ceil(succ[0]/2)) for succ in V[0].succ], ['{}'.format([succ[1]]) for succ in V[0].succ]]
            b1_line.append('\n')
            b1_line = [x for sublist in b1_line for x in sublist]
            raw_b1_line = '\t'.join(b1_line)
            f.write(raw_b1_line)
            for i in tasks:
                b1_line = [['{}'.format(i)], '1', ['{}'.format(len(tasks[i].SS_succ))], ['{}'.format(len(tasks[i].SF_succ))], ['{}'.format(len(tasks[i].FS_succ))], ['{}'.format(len(tasks[i].FF_succ))], ['{}'.format(succ[0]) for succ in tasks[i].SS_succ], ['{}'.format([succ[1]]) for succ in tasks[i].SS_succ], ['{}'.format(succ[0]) for succ in tasks[i].SF_succ], ['{}'.format([succ[1]]) for succ in tasks[i].SF_succ], ['{}'.format(succ[0]) for succ in tasks[i].FS_succ], ['{}'.format([succ[1]]) for succ in tasks[i].FS_succ], ['{}'.format(succ[0]) for succ in tasks[i].FF_succ], ['{}'.format([succ[1]]) for succ in tasks[i].FF_succ]]
                b1_line.append('\n')
                b1_line = [x for sublist in b1_line for x in sublist]
                raw_b1_line = '\t'.join(b1_line)
                f.write(raw_b1_line)
            b1_line = [['{}'.format(n+1)],'1','0','0','0','0']
            b1_line.append('\n')
            b1_line = [x for sublist in b1_line for x in sublist]
            raw_b1_line = '\t'.join(b1_line)
            f.write(raw_b1_line)
            b2_line1 = ['0', '1', '0', '0', ['0' for r in range(2*n_res)]]
            b2_line1.append('\n')
            b2_line1 = [x for sublist in b2_line1 for x in sublist]
            raw_b2_line1 = '\t'.join(b2_line1)
            f.write(raw_b2_line1)
            for i in tasks:
                q = []
                for r in range(n_res):
                    q.append('{}'.format(tasks[i].q_min[r]))
                    q.append('{}'.format(tasks[i].q_max[r]))
                b2_line = [['{}'.format(i)], '1', '0', ['{}'.format(tasks[i].w)], q]
                b2_line.append('\n')
                b2_line = [x for sublist in b2_line for x in sublist]
                raw_b2_line = '\t'.join(b2_line)
                f.write(raw_b2_line)
            b2_last_line = [['{}'.format(n+1)], '1', '0', '0', ['0' for r in range(2*n_res)]]
            b2_last_line.append('\n')
            b2_last_line = [x for sublist in b2_last_line for x in sublist]
            raw_b2_last_line = '\t'.join(b2_last_line)
            f.write(raw_b2_last_line)
            last_line = []
            for r in range(n_res):
                last_line.append('{}'.format(R_max[r]))
            last_line.append('{}'.format(l_min))
            raw_last_line = '\t'.join(last_line)
            f.write(raw_last_line)
            f.close()
            
            stats['gen_t'] = datetime.datetime.now().replace(microsecond=0)
            
            stat_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(stats['name'],stats['gen_t'],stats['n_act'],stats['n_sources'],stats['n_sinks'],stats['avg_pred'],stats['n_nred'],stats['OS'],stats['n_red'],stats['n_mtl'],stats['n_cycles'],stats['avg_cs'],stats['min_cs'],stats['max_cs'],stats['n_res'],stats['avg_w'],stats['avg_dmin'],stats['avg_dmax'],stats['rf'],stats['RS'],stats['net_LB'], stats['T_max']))

stat_file.close()

