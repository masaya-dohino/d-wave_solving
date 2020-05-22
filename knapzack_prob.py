#d-wave_solving
import numpy as np
from pyqubo import Array, Constraint, Placeholder, solve_qubo
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import time

np.random.seed(2)
# BINARY変数
NUM_VER = 10
x = Array.create('x', shape=NUM_VER, vartype='BINARY')

limit = 2*NUM_VER
y = Array.create('y',shape=limit + 1, vartype='BINARY')

value = np.random.randint(1, NUM_VER, NUM_VER)
print("value")
print(value)

weight = np.random.randint(1, NUM_VER, NUM_VER)
print("weight")
print(weight)

detail = False

#使いたい計算
def func1(limit,a):
    sum1 = 0
    for i in range(1,limit+1):
        sum1 += i*a[i]
    
    return sum1

def func2(N,p,x):
    sum2 = 0
    for i in range(N):
        sum2 += p[i]*x[i]
    
    return sum2


# プレースホルダー
param_cover1 = Placeholder("cover1")
param_cover2 = Placeholder("cover2")

#ナップザック問題におけるハミルトニアン
H_A = param_cover1 * Constraint(((1 - sum(y)))**2,"cover1") + param_cover1 * Constraint(((func1(limit,y) - func2(NUM_VER,weight,x)))**2,"cover2")
H_B = param_cover2 * func2(NUM_VER,value,x)
H = H_A - H_B


# モデルをコンパイル
model = H.compile()

# 制約項のペナルティウェイト
param_cover1 = 3.0
param_cover2 = 1.0

# プレースホルダーと合わせてQUBOを作成
feed_dict = {"cover1": param_cover1,"cover2": param_cover2}
qubo, offset = model.to_qubo(feed_dict=feed_dict)

#dwave 計算

sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_qubo(qubo,num_reads = 100)

if detail == True:



    df_result = pd.DataFrame()

    k = 0

    for sample, energy, num_occurrences, chain_break_fraction in list(response.data()):

        #print(sample, "Energy: ", energy, "Occurrences: ", num_occurrences)

        df_tmp = pd.DataFrame(dict(sample), index=[k])

        df_tmp['Energy'] = -energy

        df_tmp['Occurrences'] = num_occurrences

        df_result = df_result.append(df_tmp)

        k+=1

    

    result = df_result.pivot_table(index=df_result.columns[:n].tolist()+['Energy'], values=['Occurrences'], aggfunc='sum').sort_values('Energy', ascending=False)

        

    print(result)

    

else:

    print(list(response.data())[0])
