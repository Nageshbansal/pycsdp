
# import pycsdp.utilites.pca
from pycsdp.utilites import pca
import numpy as np
from joblib import Parallel,delayed
import multiprocessing
import autograd.numpy as anp
from pymoo.model.problem import Problem
import math
from pymoo.algorithms.hfidea_d_non_dom import HFIDEAD
from pymoo.optimize import minimize
from pymoo.factory import get_performance_indicator
import pandas as pd
from pymoo.util.normalization import normalize as norm
import hvwfg
from pymoo.util.function_loader import load_function
from pymoo.factory import get_problem, get_crossover, get_mutation, get_reference_directions, get_selection
seeds = 1 + np.arange(31)
print(seeds)
def get_ref(M):
    if M==2:
        return get_reference_directions('das-dennis',2,n_partitions=99)
    elif M==3:
        return get_reference_directions('das-dennis',3,n_partitions=13)
    elif M==4:
        return get_reference_directions('das-dennis',4,n_partitions=8)
    elif M==5:
        return get_reference_directions('das-dennis',5,n_partitions=6)
    elif M==6:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 6, n_partitions=4, scaling=1.0), get_reference_directions("das-dennis", 6, n_partitions=4, scaling=0.5))
    elif M==7:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 7, n_partitions=4, scaling=1.0), get_reference_directions("das-dennis", 7, n_partitions=3, scaling=0.5))
    elif M==8:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 8, n_partitions=3, scaling=1.0), get_reference_directions("das-dennis", 8, n_partitions=3, scaling=0.5))
    elif M==9:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 9, n_partitions=3, scaling=1.0), get_reference_directions("das-dennis", 9, n_partitions=3, scaling=0.5))
    elif M==10:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 10, n_partitions=3, scaling=1.0), get_reference_directions("das-dennis", 10, n_partitions=3, scaling=0.5))
    elif M==11:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 11, n_partitions=2, scaling=1.0), get_reference_directions("das-dennis", 11, n_partitions=2, scaling=0.5))
    elif M==12:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 12, n_partitions=2, scaling=1.0), get_reference_directions("das-dennis", 12, n_partitions=2, scaling=0.5))
    elif M==13:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 13, n_partitions=2, scaling=1.0), get_reference_directions("das-dennis", 13, n_partitions=1, scaling=0.5))
    elif M==14:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 14, n_partitions=2, scaling=1.0), get_reference_directions("das-dennis", 14, n_partitions=1, scaling=0.5))
    elif M==15:
        return get_reference_directions("multi-layer", get_reference_directions("das-dennis", 15, n_partitions=2, scaling=1.0), get_reference_directions("das-dennis", 15, n_partitions=1, scaling=0.5))
for i in range(2,16):
    print(i,len(get_ref(i)))
print("worked")

from pymoo.model.callback import Callback
class MyCallback(Callback):

    def __init__(self) -> None:
        super().__init__()
        self.data["pop"] = []
#         self.data["nadir-update"] = []
#         self.data["insights"] = []
#         self.data["term-gen"] = []

    def notify(self, algorithm):
        self.data["pop"].append(algorithm.pop)
#         self.data["nadir-update"].append(algorithm.term_gen)
#         self.data["insights"].append(algorithm.data)
#         self.data["term-gen"] = algorithm.termination_suggestion

class MyCallback1(Callback):

    def __init__(self) -> None:
        super().__init__()
        self.data["pop"] = []
#         self.data["nadir-update"] = []
#         self.data["insights"] = []
        self.data["term-gen"] = []

    def notify(self, algorithm):
        self.data["pop"].append(algorithm.pop)
#         self.data["nadir-update"].append(algorithm.term_gen)
#         self.data["insights"].append(algorithm.data)
        self.data["term-gen"] = algorithm.termination_suggestion

import autograd.numpy as anp
import math

from pymoo.model.problem import Problem

import autograd.numpy as anp
from pymoo.model.problem import Problem
class DTLZ5_3_5(Problem):

    def __init__(self, ess):

        # define lower and upper bounds -  1d array with length equal to number of variable
        xl = 0 * anp.ones(14)
        xu = 1 * anp.ones(14)

        super().__init__(n_var=14, n_obj=len(ess), n_constr=3, xl=xl, xu=xu, elementwise_evaluation=True)
        self.ess = ess

    def _evaluate(self, x, out, *args, **kwargs):
        M = 5
        I = 3
        noc = M-I+1
        k = 1
        p = np.zeros(M-I+1)
        p = [(M-I+2)-(i+1) for i in range(len(p))]
        p[0] = M-1
        g = sum([(x[i]-0.5)**2 for i in range(M-1,M+k-1)])
        theta = np.zeros(M-1)
        for i in range(I-1):
            theta[i] = math.pi*x[i]/2
        for i in range(I-1,M-1):
            theta[i] = math.pi*(1+2*g*x[i])/(4*(1+g))
        F = np.zeros(M)
#         val = 0.5*(1+g)
#         for i in range(M-1):
#             val = val*math.cos(theta[i])
#         F[0] = val
#         F[M-1] = (1+g)*math.sin(theta[0])
#         for i in range(1,M-1):
#             F[i] = (1+g)*math.sin(theta[M-(i+1)])
#             for j in range(M-i):
#                 F[i] = F[i]*math.cos(theta[j])
        for i in range(M):
            if i==0:
                sine = 1
            else:
                sine = math.sin(theta[M-1-i])
            coss = 1
            for j in range(M-1-i):
                coss = coss*math.cos(theta[j])
            F[i] = (1+g)*coss*sine
                
                
                
        G = np.zeros(len(p))
#         for i in range(len(p)):
#             val = 0
#             for j in range(I-2+1):
#                 val += F[M-1-(j)]**2
#             G[i] = 1 - val - (2**p[i])*(F[i]**2)
            
        ind = 0
        for j in range(I-1):
            ind += F[M-j-1]**2
        for i in range(M-I+1):
            if i==0:
                p=M-I
            else:
                p=M-I+1-i
            G[i] = -ind +1 -(2**p)*(F[i]**2)

        out["F"] = anp.column_stack(F)[:,self.ess]
        out["G"] = anp.column_stack(G)

class GEARBOX(Problem):

    def __init__(self, ess):

        # define lower and upper bounds -  1d array with length equal to number of variable
        xl = np.array([18]*24)
        xu = np.array([40]*24)
        
        xl[23] = 0
        xu[23] = 1
        xl[9] = 0
        xu[9] = 50
        
        for i in range(9):
            xl[i], xu[i] = 0.5,2.5

        super().__init__(n_var=24, n_obj=3, n_constr=52, xl=xl, xu=xu,  elementwise_evaluation=True)
        self.ess = ess

    def _evaluate(self, x, out, *args, **kwargs):
        Kc = 1.5
        Kd = 1.1
        r_max = 3.5
        beta = math.pi/9.0
        Sb = 2500.0
        Sw = 17500.0 
        YM_E = 2100000.0
        G = int(9)
        
        r_i = np.zeros(G)
        t_i = np.zeros(G)
        a_i = np.zeros(G)
        w_i = np.zeros(G)
        nw_i = np.zeros(G)
        y_i = np.zeros(G)
        n = np.zeros(2*G)
        omega = np.zeros(40)
        out_rpm = np.zeros(2*G)
        com_sb_sw = np.zeros(G)
        true_Sb = np.zeros(G)
        true_Sw = np.zeros(G)
        
        power = x[9]
        module = x[23]
        
        n[0] = x[10]
        n[1] = x[11]
        n[2] = x[12]
        n[3] = n[0] + n[1] - n[2]
        n[4] = x[13]
        n[5] = n[0] + n[1] - n[4]
        n[6] = x[14]
        n[7] = x[15]
        n[8] = x[16]
        n[9] = n[6] + n[7] - n[8]
        n[10] = x[17]
        n[11] = n[6] + n[7] - n[10]
        n[12] = x[18]
        n[13] = x[19]
        n[14] = x[20]
        n[15] = x[21]
        n[16] = x[22]
        n[17] = n[14] + n[15] - n[16]
        
        omega[0]= 1400.0
        omega[1]= 1400.0 * (n[0]/n[1])
        omega[2]= 1400.0 * (n[2]/n[3])
        omega[3]= 1400.0 * (n[4]/n[5])
        
        for i in range(3):
            omega[4+(3*i)]= omega[1+i] * (n[6]/n[7])
            omega[5+(3*i)]= omega[1+i] * (n[8]/n[9])
            omega[6+(3*i)]= omega[1+i] * (n[10]/n[11])
            
        for i in range(9):
            omega[13+i]= omega[4+i] * (n[12]/n[13])
            
        for i in range(9):
            omega[22+(2*i)]= omega[13+i] * (n[14]/n[15])
            omega[23+(2*i)]= omega[13+i] * (n[16]/n[17])
            
        for i in range(3):
            w_i[i]= omega[1+i] if omega[0]>omega[1+i] else omega[0]

        for i in range(3):
            w_i[3+i]= omega[4+i] if omega[1]>omega[4+i] else omega[1]

        w_i[6]= omega[13] if omega[4]>omega[13] else omega[4]

        for i in range(2):
            w_i[7+i]= omega[22+i] if omega[13]>omega[22+i] else omega[13]

        sum_over_G = 0.0
        
        for i in range(9):
            sum_over_G += x[i] + math.pow(n[2*i],2) + math.pow(n[2*i+1],2)
            r_i[i] = n[2*i]/n[2*i+1] if n[2*i] > n[2*i+1] else n[2*i+1]/n[2*i]
            t_i[i] = x[i]
            a_i[i] = module * (n[2*i] + n[2*i+1])/2.0
            com_sb_sw[i] = (97500.0 * power * Kc * Kd * (r_i[i]+1.0))/(w_i[i] * t_i[i] * math.cos(beta))
            nw_i[i] = (n[2*i]) if n[2*i] > n[2*i+1] else (n[2*i+1])
            y_i[i] = 0.52 * (1 + 20/nw_i[i])
            
        for i in range(9):
            true_Sb[i] = com_sb_sw[i]/(a_i[i] * module * r_i[i] * y_i[i])
            true_Sw[i] = ((0.59 * (r_i[i]+1))/(r_i[i] * a_i[i])) * math.pow(((com_sb_sw[i]* YM_E)/(2.0 * math.sin(beta))), 0.5)
            
        f1 = (((math.pi*module*module)/4.0)*sum_over_G)
        f2 = -power
        f3 = (module/2.0)*((n[0]+n[1]) + (n[6]+n[7]) + (n[12]+n[13]) + (n[14]+n[15]))
        
        g1 = n[0] + n[1] -  n[10]
        g2 = n[6] + n[7] -  n[1]
        g3 = n[6] + n[7] -  n[12]
        g4 = n[12] + n[13] -  n[7]
        g5 = n[12] + n[13] -  n[16]
        g6 = n[14] + n[15] -  n[13]
        
        temp = []
        for i in range(9):
            temp.append(3.5 - r_i[i])
        [g7, g8, g9, g10, g11, g12, g13, g14, g15] = temp
        
        modulus_diff = (omega[39] - 1400.0) if omega[39] > 1400.0 else (1400.0 - omega[39])
        g16 = 1.0 - (modulus_diff/(1400.0*0.1))
        
        for i in range(3):
            out_rpm[0+i] = omega[22+(2*i)]

        for i in range(3):
            out_rpm[3+i] = omega[23+(2*i)]

        for i in range(3):
            out_rpm[6+i] = omega[28+(2*i)]

        for i in range(3):
            out_rpm[9+i] = omega[29+(2*i)]

        for i in range(3):
            out_rpm[12+i] = omega[34+(2*i)]

        for i in range(3):
            out_rpm[15+i] = omega[35+(2*i)]
            
        temp = [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16]
        for i in range(17):
            out_rpm_ratio = out_rpm[1+i]/out_rpm[i]
            modulus_diff = (out_rpm_ratio - 1.14) if out_rpm_ratio > 1.14 else (1.14 - out_rpm_ratio)
            temp.append(1.0 - (modulus_diff/(0.1*1.14)))
            
        modulus_diff = (omega[22] - 150.0) if omega[22] > 150.0 else (150.0 - omega[22])
        temp.append(1.0 - (modulus_diff/(150.0 * 0.1)))
        
        for i in range(9):
            temp.append(Sb - true_Sb[i])
            
        for i in range(9):
            temp.append(Sw - true_Sw[i])
            
        temp = -1 * np.array(temp)

        out["F"] = np.array([f1,f2,f3])[self.ess]
        out["G"] = anp.column_stack(temp)

class WATER(Problem):

    def __init__(self, ess):

        # define lower and upper bounds -  1d array with length equal to number of variable
        xl = 0.01 * anp.ones(3)
        xu = 0.1 * anp.ones(3)
        xu[0] = 0.45

        super().__init__(n_var=3, n_obj=len(ess), n_constr=7, xl=xl, xu=xu, elementwise_evaluation=True)
        self.ess = ess

    def _evaluate(self, x, out, *args, **kwargs):
        x1 = x[0]
        x2 = x[1]
        x3 = x[2]
        f1 = 106780.37*(x2+x3)+61704.67
        f2 = 3000*x1
        f3 = (305700*2289*x2)/((0.06*2289)**0.65)
        f4 = 250*2289*math.exp(-39.75*x2+9.9*x3+2.74)
        f5 = 25*(1.39/(x1*x2)+4940*x3-80)
        g1 = 0.00139/(x1*x2)+4.94*x3-0.08-1
        g2 = 0.000306/(x1*x2)+1.082*x3-0.0986-1
        g3 = 12.307/(x1*x2)+49408.24*x3+4051.02-50000
        g4 = 2.098/(x1*x2)+8046.33*x3-696.71-16000
        g5 = 2.138/(x1*x2)+7883.39*x3-705.04-10000
        g6 = 0.417*(x1*x2) + 1721.26*x3-136.54-2000
        g7 = 0.164/(x1*x2) + 631.13*x3 - 54.48 - 550

        out["F"] = np.array([f1,f2,f3,f4,f5])[self.ess]
        out["G"] = anp.column_stack([g1,g2,g3,g4,g5,g6,g7])

class CARSIDE(Problem): #GM problem

    def __init__(self, ess):

        super().__init__(n_var=7, n_obj=len(ess), n_constr=0, type_var=anp.double)
        self.xl = anp.array([0.5, 0.45, 0.5, 0.5, 0.875, 0.4, 0.4])
        self.xu = anp.array([1.5, 1.35, 1.5, 1.5, 2.625, 1.2, 1.2])
        self.ess = ess

    def _evaluate(self, x, out, *args, **kwargs):
        g1 = 1.16 - 0.3717 * x[:,1] * x[:,3] - 0.0092928 * x[:,2]
        g2 = 0.261 - 0.0159 * x[:,0] * x[:,1] - 0.188 * x[:,0] * 0.345 - 0.019 * x[:,1] * x[:,6] + 0.0144 * x[:,2] * x[:,4] + 0.08045 * x[:,5] * 0.192
        g3 = 0.214 + 0.00817 * x[:,4] - 0.131 * x[:,0] * 0.345 - 0.0704 * x[:,0] * 0.192 + 0.03099 * x[:,1] * x[:,5] - 0.018 * x[:,1] * x[:,6] + 0.0208 * x[:,2] * 0.345 + 0.121 * x[:,2] * 0.192 - 0.00364 * x[:,4] * x[:,5] - 0.018 * x[:,1] ** 2
        g4 = 0.74 - 0.61 * x[:,1] - 0.031296 * x[:,2] - 0.166 * x[:,6] * 0.192 + 0.227 * x[:,1] ** 2
        g5 = 28.98 + 3.818 * x[:,2] - 4.2 * x[:,0] * x[:,1] + 6.63 * x[:,5] * 0.192 - 7.77 * x[:,6] * 0.345
        g6 = 33.86 + 2.95 * x[:,2] - 5.057 * x[:,0] * x[:,1] - 11 * x[:,1] * 0.345 - 9.98 * x[:,6] * 0.345 + 22 * 0.345 * 0.192
        g7 = 46.36 - 9.9 * x[:,1] - 12.9 * x[:,0] * 0.345
        g8 = 4.72 - 0.5 * x[:,3] - 0.19 * x[:,1] * x[:,2]
        g9 = 10.58 - 0.674 * x[:,0] * x[:,1] - 1.95 * x[:,1] * 0.345
        g10 = 16.45 - 0.489 * x[:,2] * x[:,6] - 0.843 * x[:,4] * x[:,5]

        f1 = 1.98 + 4.9 * x[:,0] + 6.67 * x[:,1] + 6.98 * x[:,2] + 4.01 * x[:,3] + 1.78 * x[:,4] + 0.00001 * x[:,5] + 2.73 * x[:,6]

        F = anp.column_stack([f1, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10])
        out["F"] = np.array([row[self.ess] for row in F])


ngen= 1

names = ['GEARBOX'] # change problem name

for i in range(1):
    # building the pca libraries (nl_pca or l_pca)
    nl_fn, free_fs = pca.pca_reduction('nl_pca','Release')
    #######
    algo = 'HFIDEAD'
    name = names[i]
    n_obj = 3 # change number of objectives
    print(name, ngen)
    crossover = get_crossover('real_sbx',prob=0.9,eta=20)
    mutation = get_mutation('real_pm',prob=1/15,eta=20)
    ref_dirs = get_ref(n_obj)
    algorithm = HFIDEAD(pop_size=len(ref_dirs),ref_dirs=ref_dirs,crossover=crossover,mutation=mutation,eliminate_duplicates=True)
    def optimize(seed_value):
        global problem,algorithm,ngen,n_obj, name
        ess_init = [i for i in range(n_obj)]
        pfs = []
        term_gens = []
        for interation in range(n_obj):
            problem = GEARBOX(ess_init) #change problem name
            M = len(ess_init)
            ref_dirs = get_ref(M)
            algorithm = HFIDEAD(pop_size=len(ref_dirs),ref_dirs=ref_dirs,crossover=crossover,mutation=mutation,eliminate_duplicates=True)
            res = minimize(problem,algorithm,('n_gen', ngen),callback=MyCallback1(),seed=seed_value)
            df = np.array(res.algorithm.callback.data["pop"])
            pf = df[-1]
            pfs.append(pf)
            term_gens.append(len(df))
            pf = np.array([ind.F for ind in pf])

            # filter only the natural cluster representatives - check the number of clusters and take the initial solutions
            pfn = norm(pf)
            dist_matrix = load_function('calc_perpendicular_distance')(pfn, ref_dirs) # i,j where i is solution and j is RV
            assoc_index = np.array([np.argmin(row) for row in dist_matrix])
            count = len(np.unique(assoc_index))
            print(count)
            pf = pf[:count]

            # np.savetxt('pf_data_nl.txt', pf)
            # out = !/home/sukrit/Documents/HFiDEA/NL_MVU_PCA/Release/nl_mvu_pca_dmoss_kemoss ./pf_data_nl.txt
            # out = out[0][7:]
            # out = out.split(' ')
            # calling the nl_pca function
            out = pca.nl_pca_reduction(pf, nl_fn, free_fs)
            # print(out)
            ess = np.array([int(row)-1 for row in out])
            if len(ess_init)==len(ess):
                break
            else:
                ess_init = np.array(ess_init)
                ess_init = ess_init[ess]
        np.save('data/'+name+'-fronts-'+str(seed_value)+'.npy',pfs)
        np.save('data/'+name+'-gens-'+str(seed_value)+'.npy',term_gens)
    for seed in seeds:
        optimize(seed)
    pca._clean_pca('nl_pca')
#     Parallel(n_jobs=7)(delayed(optimize)(i) for i in seeds[:])



