#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import *

def G2D(G):
	L = G.shape[0]
	D = zeros((L * L,L * L),float)
	i = 0
	while(i != L):	
		j = 0
		while(j != L):		
			if(G[i,j] ==  0):
				m = 0
				while(m != L):		
					n =  0
					while(n != L):
						if(G[m,n] == 0):
							im = float(abs(i - m))
							jn = float(abs(j - n))
							if(((im + jn) == 1) | ((im == 1) & (jn == 1))):
								D[(i)*L+j,(m)*L+n] = (im+jn)**0.5
						n += 1
					m += 1
			j += 1
		i += 1
	return D
	
	