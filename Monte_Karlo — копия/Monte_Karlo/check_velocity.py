from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import datetime
import math
from math import fabs

def read_fiel(name):
        arr =[]
        with open(name,"r") as f:
                reader = f.read().splitlines()
                for reader1 in reader: 
                        for number in reader1.split(" "):
                                arr.append(float(number))
        f.close()
        return arr
    
def f(x,y,z,U_H):
    print(fabs(x+U_H),math.exp(-x*x-y*y-z*z))
    return fabs(x+U_H)* math.exp(-x*x-y*y-z*z)

k = 1
for t in range(0,k):
        arr =[]
        arr1 =[]
        arr2 = []
        name1 ='data0/'+str(0)+str(0)+str(t)+'.txt'
        #name2 ='data0/'+str(0)+str(0)+str(0)+'_.txt'

        arr = read_fiel(name1)
        #arr2 = read_fiel(name2)
        #dx=arr2[0]/arr2[1]
        #v_x = dx
        for i in range(0,len(arr)-1,3):
                d= f(arr[i],arr[i+1],arr[i+2],-26)
                if arr[i] < -24.5: 
                    arr1.append(d)
                    arr2.append(arr[i])
        #        v_x = v_x + dx
        print(len(arr2))
        print(len(arr1))
        plt.plot(arr2,arr1,'ro',arr2,arr1)
        #name = 'CH = ' + str(arr2[0])+ ' N_CH = ' + str(arr2[1]) + ' U_H = ' + str(arr2[2]) + ' L = ' + str(arr2[3]) + 'v_x = '+str(arr2[4])+' delta = '+ str(arr2[5])+' Cp = '+str(arr2[6])            
        #plt.xlabel( name )                        
        plt.savefig('data0/'+str(0)+str(0)+str(t)+'.png')
        


print('end')
