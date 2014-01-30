# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 22:31:12 2014

@author: martin
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 10:44:42 2013
@author: martin + julian
"""
#Using the magic encoding
#-*- coding: utf-8 -*-
from scipy import * 
import matplotlib.pyplot as p


from uncertainties import *
import math


def make_LaTeX_table(data,header, flip= 'false', onedim = 'false'):
    output = '\\begin{table}\n\\centering\n\\begin{tabular}{'
    #Get dimensions
    if(onedim == 'true'):
        if(flip == 'false'):
        
            data = array([[i] for i in data])
        
        else:
            data = array([data])
    
    row_cnt, col_cnt = data.shape
    header_cnt = len(header)
    
    if(header_cnt == col_cnt and flip== 'false'):
        #Make Format
        
        for i in range(col_cnt):
            output += 'l'
        output += '}\n\\toprule\n{'+ header[0]
        for i in range (1,col_cnt):
            output += '} &{ ' + header[i]
        output += ' }\\\\\n\\midrule\n'
        for i in data:
            if(isinstance(i[0],(int,float,int32))):
                output +=  str( i[0] ) 
            else:
                output += ' ${:L}$ '.format(i[0])
            for j in range(1,col_cnt):
                if(isinstance(i[j],(int,float,int32))):
                    output += ' & ' + str( i[j])   
                else:          
                    output += ' & ${:L}$ '.format(i[j])                
                
            output += '\\\\\n'
        output += '\\bottomrule\n\\end{tabular}\n\\label{}\n\\caption{}\n\\end{table}\n'
                            
        return output

    else:
        return 'ERROR'



    
def err(data):
    mean = data.mean()
    N = len(data)
    err = 0
    for i in data:
        err += (i - mean)**2
    err = sqrt(err/((N-1)*N))
    return ufloat(mean,err)


def lin_reg(x,y):
    N = len(x)
    sumx = x.sum()
    sumy = y.sum()
    sumxx = (x*x).sum()
    sumxy = (x*y).sum()
    m = (sumxy -  sumx*sumy/N)/(sumxx- sumx**2/N)
    b = sumy/N - m*sumx/N
    
    sy = sqrt(((y - m*x - b)**2).sum()/(N-1))
    m_err = sy *sqrt(N/(N*sumxx - sumx**2))
    b_err= m_err * sqrt(sumxx/N)
    return ufloat(m,m_err),ufloat(b,b_err)
    

import sympy as sym

def error(f, err_vars=None):
    from sympy import Symbol, latex
    s = 0
    latex_names = dict()    
    if err_vars == None:
        err_vars = f.free_symbols      
    for v in err_vars:
        err = Symbol('latex_std_' + v.name)
        s += f.diff(v)**2 * err**2
        latex_names[err] = '\\sigma_{' + latex(v) + '}'        
    return latex(sqrt(s), symbol_names=latex_names)

###################################################

# Daten erste Messung

m = ufloat(0.5122,0.0002)
r_k = 0.5*ufloat(50.76,0.0035532)*1e-3
print 'R_K %s'% r_k
theta_K = 22.5e-7
r = 0.5*ufloat(0.18e-3,0.01e-3)
L=0.63


theta = (2*m*r_k**2)/5+theta_K

print theta

x, y = sym.var('m_k r_k')

#print error(2*x*y**2/5)

T1= array([18.474,18.471,18.473,18.473,18.463,18.466,18.466,18.462,18.469,18.460,18.460,18.454,18.459,18.453,18.457,18.435,18.457,18.453,18.457,18.447,18.450])

T_quer = err(T1)
print T_quer

G = 8* pi * theta * L /(r**4*T_quer**2)

print G




# T bzgl versch. StromstÃ¤rken in A

T_b =array([[10.402,10.391,10.381,10.380,10.364],[8.701,8.701,8.682,8.679,8.673],[7.371,7.362,7.355,7.355,7.349],[6.674,6.671,6.669,6.665,6.662],[6.123,6.124,6.120,6.113,6.120],[5.653,5.643,5.643,5.648,5.640],[5.286,5.256,5.272,5.253,5.247],[5.103,5.093,5.085,5.077,5.071],[4.795,4.772,4.784,4.772,4.777],[4.526,4.527,4.517,4.516,4.522]])

T_b = array([err(i) for i in T_b ])

I = array([.1,.2,.3,.4,.5,.6,.7,.8,.9,1 ])

B=1.25663706147e-6 * 8*I*390/(sqrt(125)*0.078)
print B

data = array([I,[round(i,6) for i in B],T_b])
m,b = lin_reg(B,1/array([i.n for i in T_b])**2)
print make_LaTeX_table(data.T, [r'$\frac{I}{\si{\ampere}}$',r'$\frac{B}{\si{\tesla}}$',r'$\frac{T}{\si{\second}}$'])
print "b: %s und m: %s"% (b,m)

mag = 4*pi**2*theta* m

print "Das magnetische Moment beträgt : %s " % mag
p.close()
x = linspace(0,.005)

p.plot(B,1/array([i.n for i in T_b])**2, 'x', label = 'Messwerte')
p.plot(x,m.n*x+b.n, label = "Lineare Regression")
p.ylim(0,0.05)
p.xlabel('B in T')
p.ylabel(r' $T^{-2}$ in $s^{-2}$',rotation=0)
p.legend(loc="bootom")
p.subplots_adjust(left = 0.2)
p.savefig('./Abb4.png')


# Daten zweite Messung

T2=array([18.876,18.906,18.908,18.913,18.913,18.913,18.914,18.916,18.912,18.921,18.923,18.921,18.926,18.918,18.934,18.925,18.949,18.928,18.944,18.29])

T2 = err(T2)
print T2
D = pi*G*r**4/(2*L)
B_erd = -(theta/(4*pi**2*T2**2)-D)/mag
print "Dasn Erdmagnetfeld beträgt %s Tesla" % B_erd

