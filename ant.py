#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import *
from pylab import *
import random
import sys
import G2D
import copy
import time

def random_int_list(start, stop, length):
    start, stop = (int(start), int(stop)) if start <= stop else (int(stop), int(start))
    length = int(abs(length)) if length else 0
    random_list = []
    for i in range(length):
        random_list.append(random.randint(start, stop))
    return random_list

NNUM = 25                                          #地图边长

G = matrix([[0,1,0,0,0,0,0],            #自定义地图 1为障碍，边长要对应地图边长
						[0,1,0,0,0,1,1],
						[0,0,0,0,0,1,1],
						[0,1,0,1,0,1,1],
						[0,1,0,1,0,1,0],
						[1,0,0,0,1,0,1],
						[1,1,1,0,0,0,0]])

#```                 以下为随机地图			若选择自定义地图，去掉前面的#号和后面的#号			
prar = random_int_list(1,NNUM*NNUM - 2,NNUM * 8)          #最后一个参数为障碍数量，地图边长变短时可能需要适当减少，不然可能生成不恰当的地图
G = zeros((NNUM,NNUM));
for pp in prar:
	x = pp//NNUM
	y = pp%NNUM
	G[x][y] = 1
print(prar)               #障碍的位置
#```        随机地图生成结束

time_start=time.time()
MM = G.shape[0]   #G 地形图为01矩阵，如果为1表示障碍物
Tau = ones((MM*MM,MM*MM)) # Tau 初始信息素矩阵（认为前面的觅食活动中有残留的信息素）
Tau = 8.*Tau
K = 50   #K 迭代次数（指蚂蚁出动多少波）
M = 30   #M 蚂蚁个数（每一波蚂蚁有多少个）
S = 0   #S 起始点（最短路径的起始点）
E = MM*MM - 1   #E 终止点（最短路径的目的点）
Alpha = 1    # Alpha 表征信息素重要程度的参数
Beta = 7     #Beta 表征启发式因子重要程度的参数
Rho = 0.3      # Rho 信息素蒸发系数
Q = 1         # Q 信息素增加强度系数
minkl = inf
mink = 0
minl = 0
load = ''             #进度条显示
D = G2D.G2D(G)

N = D.shape[0] 
a = 1
Ex = a * (mod(E,MM - 1) - 0.5)   #终止点横坐标
if(Ex == -0.5):
	Ex = MM - 0.5
Ey = a * (MM + 0.5 - ceil(E/MM))  #终止点纵坐标
Eta = zeros((N,N))

#构造启发式信息矩阵
i = 0
while(i != N):
	ix = a*(mod(i,MM - 1) - 0.5)
	if(ix == -0.5):
		ix = MM - 0.5
	iy = a * (MM + 0.5 - ceil(i/MM))
	if (i != E):
		Eta[i][0]=1/(((ix - Ex)**2+(iy-Ey)**2)**0.5)
	else:
		Eta[i][0] = 100
	i += 1

ROUTES = []                     #路径储存 包含每一代每一只蚂蚁的路径
PL = zeros((K,M))          #路径长度 同上


#正式开始 有K代蚂蚁 每代有M只蚂蚁
k = 0
while(k != K):
	m = 0
	while(m != M):
		#初始化
		W = S      #当前节点
		Path = []
		Path = [S]       #初始化路径
		PLkm = 0
		TABUkm = ones(N)
		TABUkm[S] = 0
		DD = copy.deepcopy(D)
		#检测下一步可以走的节点
		DW = DD[W,:]
		DW1 = nonzero(DW)
		j = 0
		while(j != len(DW1[0])):
			non0 = DW1[0][j]
			if(TABUkm[non0] == 0):
				DW[non0] = 0
			j += 1
		LJD = nonzero(DW)
		Len_LJD = len(LJD[0])    #计算出可以走的节点数量
		while((W != E) & (Len_LJD >= 1) ):              #如果无路可走或走到食物就退出
		#用转轮赌法选择下一步
			PP = zeros((Len_LJD,Len_LJD))
			i = 0
			while(i !=  Len_LJD):
				PP[i][0] = (Tau[W,LJD[0][i]]**Alpha)*((Eta[LJD[0][i]][0])**Beta)
				i += 1
			sumPP = sum(PP)
			PP = PP/sumPP
			Pcum = zeros(Len_LJD)
			Pcum[0] = PP[0][0]
			i = 1
			while(i != Len_LJD):
				Pcum[i] = Pcum[i - 1] + PP[i][0]
				i += 1
			rand = random.random()
			Select = nonzero(Pcum >= rand)
			to_visit = LJD[0][Select[0][0]]
			
			
			# 更新状态
			Path.append(to_visit) #更新路径
			PLkm = PLkm + DD[W,to_visit] #更新长度
			W = to_visit 
			kk = 0
			while(kk != N):
				if(TABUkm[kk] == 0):
					DD[W,kk] = 0
					DD[kk,W] = 0
				kk += 1
			TABUkm[W] = 0
			DW = DD[W,:]
			DW1 = nonzero(DW)
			j = 0
			while(j != len(DW1[0])):
				if(TABUkm[DW1[0][j]] == 0):
					DW[j] = 0
				j += 1
			LJD = nonzero(DW)
			Len_LJD = len(LJD[0])
			
		#记录每代的蚂蚁路径和长度
		ROUTES.append(Path)
		if (Path[-1] == E):
			PL[k][m] = PLkm
			if(PLkm < minkl):
				mink = k
				minl = m 
				minkl = PLkm
		else:
			PL[k][m] = 0
		m += 1
		
	#更新信息素
	Delta_Tau = zeros((N,N))
	m = 0
	while(m != M):
		if(PL[k][m] != 0):
			ROUT = copy.deepcopy(ROUTES[k*M + m])
			TS = len(ROUT) - 1
			PL_km = PL[k][m]
			s = 0
			while(s != TS):
				x = ROUT[s]
				y = ROUT[s + 1]
				Delta_Tau[x][y] = Delta_Tau[x][y] + Q/PL_km
				Delta_Tau[y][x] = Delta_Tau[y][x] + Q/PL_km
				s += 1
		m += 1
	Tau = (1 - Rho) * Tau + Delta_Tau #信息素矩阵更新
	load += '-' 
	print(load)   #进度显示
	k += 1
#算法结束，下面是图像显示
print(ROUTES[mink*M + minl])		#输出最短的一条路径
t = 0
xian = ''
while(t != NNUM):
	xian += '---'
	t += 1
Pic = [xian]
i = 0
while(i != NNUM):
	j = 0
	Pstr = '|'
	while(j != NNUM):
		P = False
		for q in ROUTES[mink*M + minl]:
			if((i*NNUM + j) == q ):
				P = True
		if (P):
			Pstr += ' *|'     #路径
		elif(G[i,j] == 0):
			Pstr += '  |'
		else:
			Pstr += '▇|'    #障碍
		j += 1	
	Pic.append(Pstr)
	Pic.append(xian)
	i += 1
for x in Pic:
	print(x)
time_end=time.time()
print (str(time_end-time_start) , 's') #输出程序运行时间