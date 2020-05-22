import numpy as np
from pyqubo import Array, Constraint, Placeholder, solve_qubo
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import time


goal = 0
np.random.seed(56)

NUM_VER = 8
x = Array.create('x', shape = NUM_VER, vartype='SPIN')


member = np.random.randint(1, NUM_VER, NUM_VER)
print("member")
print(member)
print("goal")
goal = sum(member) / 2
print(goal)

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

start = time.time()
#数の分割問題
H = param_cover1*(func2(NUM_VER,member,x))**2


# モデルをコンパイル
model = H.compile()

# 制約項のペナルティウェイト
param_cover1 = 1.0
param_cover2 = 1.0

# プレースホルダーと合わせてQUBOを作成
feed_dict = {'cover1': param_cover1}
qubo, offset = model.to_qubo(feed_dict=feed_dict)

#dwave 計算
sampler = EmbeddingComposite(DWaveSampler())

print("Start anealing", time.time() - start)
response = sampler.sample_qubo(qubo,num_reads = 10)
print("End anealing", time.time() - start)

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
    
print("Total_real_time ", response.info["timing"]["total_real_time"], "us")

print(time.time() - start)
print(response)
print(" ")
print(next(response.samples()))
p = dict(next(response.samples()))
ans = 0
for i in range(NUM_VER):
      ans += p["x[{}]".format(i)] * member[i]
        
print('answer')
print(ans)
