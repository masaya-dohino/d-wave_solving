# d-wave_solving
import numpy as np
from pyqubo import Array, Constraint, Placeholder, solve_qubo
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import time
import networkx as nx

#バイナリ変数作成
num = 5
x = Array.create('x', shape = (num,num), vartype = 'BINARY')


#無向グラフの作成
g = np.array([
    [0, 1, 1, 0, 0,],
    [1, 0, 1, 1, 1,],
    [1, 1, 0, 0, 1,],
    [0, 1, 0, 0, 1,],
    [0, 1, 1, 1, 0,],
   
])
nodes = np.array(['a', 'b', 'c', 'd', 'e'])
G = nx.Graph()
G.add_nodes_from(nodes)
edges = []
for hi, hv  in enumerate(g):
    for wi, wv in enumerate(hv):
        if(wv): edges.append((nodes[hi], nodes[wi]))
G.add_edges_from(edges)
pos = nx.spring_layout(G)
nx.draw_networkx(G, pos, with_labels=True)
plt.axis("off")
plt.show()

detail = False
# プレースホルダー
param_cover1 = Placeholder("cover1")

#必要な計算
def func1(N,bi):
    sum1 = 0
    ans1 = 1
    for v in range(N):
        ans1 += (1 - sum1)**2
        for j in range(N):
            sum1 += bi[v][j]
    return ans1

def func2(N,bi):
    sum2 = 0
    ans2 = 0
    for j in range(N):
        ans2 += (1 - sum2)**2
        for v in range(N):
            sum2 += bi[v][j]
    return ans2

def func3(N,bi,G):
    ans3 = 0
    for i in range(N):
        for j in range(i,N):
            for k in range(N - 1):
                if(not(G[i][j])):
                    ans3 += bi[i][k]*bi[j][k+1]
    
    return ans3

# エネルギー関数の定義（ハミルトン経路問題のイジングモデル）

H1 = func1(num,x)
H2 = func2(num,x)
H3 = func3(num,x,g)

H = param_cover1*(H1 + H2 +H3)


# モデルをコンパイル
model = H.compile()

# 制約項のペナルティウェイト
param_cover1 = 1.0

# プレースホルダーと合わせてQUBOを作成
feed_dict = {'cover1': param_cover1}
qubo, offset = model.to_qubo(feed_dict=feed_dict)

#dwave 計算

sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_qubo(qubo,num_reads =5000)　#アニーリング試行回数

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

    print(list(response.data())[0])　#出力
