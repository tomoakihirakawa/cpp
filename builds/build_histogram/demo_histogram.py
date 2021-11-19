import math
from fundamental import *
import json

f = open("data.json")
data = json.load(f)["data"]
f.close()

# class Histogram():
#     def __init__(self,dataIN):        
#         self.data = sorted(dataIN)
#         self.lendata = len(self.data)
#         self.sturges = (int)(math.log2(self.lendata)+1)

#         self.bins = [[] for i in range(self.sturges)]
#         s = self.sturges
#         min_data = min(self.data)
#         self.bin_width = (max(self.data)-min_data)/s
#         self.interval = [i*self.bin_width+min_data for i in range(s+1)]
#         for i in range(len(self.data)-1):
#             d = self.data[i]
#             self.bins[round((d-min_data)/self.bin_width-1/2.)].append(d)
#         self.bins[-1].append(self.data[-1])

#         self.der=[]
#         for i in range(0, s-1):
#             self.der.append((len(self.bins[i+1])-len(self.bins[i]))/self.lendata)

h = Histogram(data)

data = h.data
bin = h.bins
der = h.diff
s = len(bin)
print(data)
print(bin)
print("長さ：",[len(v) for v in bin])
print(h.interval)
print(der)
    
import matplotlib.pyplot as plt 

fig, ax = plt.subplots()
cum_bar = ax.bar(h.mid_interval, h.cumulative_count)
bar = ax.bar(h.mid_interval, h.count)
plot = ax.plot([h.interval[i+1] for i in range(len(bin)-1)], [i for i in der],'r-o')
ax.set_xlabel("x")
ax.set_ylabel("y")
# ax.set_xticklabels(["{:.3f}".format(h.interval[i]) for i in range(len(bin)-1)])
ax.set_xticks([x for x in h.interval])
plt.show()

# plt.bar([i for i in range(len(bin))], [len(i) for i in bin])
# plt.xlabel("x")
# plt.ylabel("y")
# plt.title("total:"+str(len(data)))
# plt.pause(10.)


# https://fraserlab.com/2014/08/20/Figures-with-Python/

