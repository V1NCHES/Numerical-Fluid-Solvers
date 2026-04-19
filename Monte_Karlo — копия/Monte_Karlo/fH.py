from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import datetime
import math
from math import fabs
import numpy as np



def read_fiel(name):
        arr =[]
        i=0
        with open(name,"r") as f:
                reader = f.read().splitlines()
                for reader1 in reader: 
                        for number in reader1.split(" "):
                                #print(number)
                                arr.append(float(number))
        f.close()
        return arr

def read_fielA(name,arr1,arr2):
        arr =[]
        i=0
        with open(name,"r") as f:
                reader = f.read().splitlines()
                for reader1 in reader: 
                        for number in reader1.split(" "):
                                if i == 0:
                                        arr1.append(float(number))
                                        i=1
                                elif i==1:
                                        arr2.append(float(number))
                                        i=2
                                elif i == 4:
                                        i=0
                                else:
                                        i+=1
                                        
        f.close()

        
def Fp(vx, Up,cp):
        return math.exp(-((vx -Up)**2)/cp**2)/(cp*math.sqrt(math.pi))
def FH(vx, UH):
        return math.exp(-((vx - UH)**2))/(math.sqrt(math.pi))

def Fp_0(vx, cp):
        return math.exp(-((vx)**2)/cp**2)/(math.sqrt(math.pi))

def integrate(c,d, vx,f, dvy , n,cp):
        integral = 0
        vy = c
        vz = c
        for i in range(n):
                vy = c
                for j in range(n):
                                integral += f(vx,vy,vz,cp)
                                vy+=dvy
                vz+=dvy	
        return dvy * dvy* integral
def method_of_rectangles(c,d, vx,f,cp):
        n = 100
        sum1 = 0
        sum2 = 0
        dvy = (d - c)/n
        sum2 = integrate(c,d, vx,f, dvy ,n,cp)
        while True:
                sum1 = sum2; sum2 = 0
                n *= 2
                dvy = (d - c)/n
                sum2 = integrate(c,d, vx,f, dvy ,n,cp)
                if not (math.fabs(sum1 - sum2) > 1e-5) and (n < 1e4):
                        break
        if (n > 1e4):
                print("ERROR n > N")
        return sum2


def grafic_fp_0(a,b, dvx, cp , f):
        fp = []
        v = a + dvx/2
        while v < b :
                fp.append(f(v,cp))
                v +=dvx
        return fp


def grafic_fp(a,b, dvx, Up,cp , f):
        fp = []
        v = a + dvx/2
        while v < b :
                fp.append(f(v,Up,cp))
                v +=dvx
        return fp

def grafic_fH(a,b, dvx, cp , f):
        fH = []
        v = a + dvx/2
        while v < b :
                fH.append(f(v,cp))
                v +=dvx
        return fH

def integral(a,b,dx,Arr):
        SUM = 0
        v = a
        for i in range(0,len(Arr),1):
                SUM += math.fabs(v)*arr[i]*dx
                #i+=1

        return SUM
        

def ox(dx,a, b):
        arr = []
        x = a+dx/2
        while x < b  :
                arr.append(x)
                x +=dx
        return arr

def index_name(index,string):
        name = string# +str(int(i/100))+str(int(i/10))+str(0)+'.txt'
        if index < 10:
                name += str(0)+str(0)+str(int(index))
        elif index < 100 and index > 9:
                name += str(0)+str(int(index / 10))
                index = index -  int(index / 10) * 10
                name += str(int(index))
        elif index > 99 and index < 1000:
                name += str(int(index / 100))
                index = index - int(index / 100) * 100
                name += str(int(index/10))
                index = index - int(index / 10) * 10
                name += str(int(index))
        return name
def picters(fp,fH, name1,dx,xt,arr2,t,arr1, name2):
        arrAx =[]
        arrAy = []
        arr =[]
        arrV = []
        #arrN = []
        #arrNx = []
        N_arr =0
        N_arrV = 0
        print(name1)
        arr = read_fiel(name1)

        read_fielA(name2,arrAx,arrAy)
        N_Arr = integrate_R(arrAy,dx)
        N_arr = integrate_R(arr,dx)      
        #arrN.append(N_arr)
        #arrNx.append(xt)
        #xt+=arr2[7]   
        SUM1 = 0
        #arrAyN = []
        for i in range(len(arr)):
                arrAy[i] = arrAy[i]/N_Arr
                arr[i] = arr[i]/N_arr
                SUM1+=dx*arr[i]
        plt.figure(figsize=(12, 12), constrained_layout=True)
        plt.xticks(np.arange(min(arr1)-dx/2, max(arr1)+dx/2, 0.5))
        plt.yticks(np.arange(min(arr)-0.05, max(arr)+0.05, 0.05))
        plt.plot(arr1,arr,'r',marker='o',label=r"$point \ f_H$")
        
        #plt.plot(arrAx,arrAy,'r',marker='o',label=r"$point \ f_HA$")
        
        plt.plot(arr1,fp,'y',marker='o',label=r"$f_p= \frac{1}{c_p\sqrt{\pi}}\cdot exp(-\frac{w_x^2}{c_p^2})$")
        #plt.plot(arr1,fH,'g',marker='o',label=r"$f_H=\frac{1}{\sqrt{\pi}}\cdot exp(-(v_x-U_H)^2)$")
        
        plt.plot(arrAx,arrAy,'b',marker='o',label=r"$A$")
        
        plt.legend(fontsize=14)
        plt.grid (True)
        plt.tick_params(labelsize=12)
        x = round(t*arr2[7],3)        
        name = 'L = ' + str(-arr2[15]) + "; x = " + str(round(x- arr2[15],4))+'; UH =  '+str(round(arr2[16],4))+'; U_p = '+ str(round(arr2[17],4))+ '; Cp = '+str(round(arr2[18],4))  + '; N = ' + str(round(arr2[14],0))         
        plt.xlabel( name, fontsize=20 )
        plt.savefig(name2+'1.png')
        plt.grid ( True )
        #plt.show()
        plt.close()
       # return xt;
        
def integrate_R(arr,dx):
        Sum = 0
        for i in range(len(arr)):
                Sum += arr[i]
        return dx*Sum;

string = 'data/dat5/'
string1 = 'dat1/dat5/'
name2 = string + str(0) + '.txt'
arr2 = []
arr1 =[]
arr3 = []
arr4 = []
arr2 = read_fiel(name2)
dx = arr2[3]
x = arr2[0]
arr1.append(x)



fH = grafic_fH(arr2[0],arr2[1], arr2[3],arr2[16],FH)
print("fH")
fp = grafic_fp(arr2[0] ,arr2[1], arr2[3],arr2[17], arr2[18],Fp)

#fp = grafic_fp_0(arr2[0] ,arr2[1], arr2[3],arr2[18],Fp_0)
print("fp")

arr1 = ox(arr2[3],arr2[0],arr2[1])

#plt.plot(arr1,fp,'y',marker='o',label=r"$f_p= \frac{1}{c_p\sqrt{\pi}}\cdot exp(-\frac{w_x^2}{c_p^2})$")
#plt.plot(arr1,fH,'g',marker='o',label=r"$f_H=\frac{1}{\sqrt{\pi}}\cdot exp(-(v_x-U_H)^2)$")
#plt.show()


#arr4 = ox(arr2[3],-2*(arr2[18]) + arr2[17],2*(arr2[18]) + arr2[17])
xt = -arr2[15]

#data = np.loadtxt('l1' + '.txt', delimiter=' ')

#read_fielA('data/data234/l1.txt',arrAx,arrAy)

arrN = []
arrNx = []

ARR_v_F=[]
SUM_v_F = 0
for t in range(0,int(arr2[6])+1):
        
        arr =[]
        arrV = []
        fp_N = []
        
        N_arr =0
        N_arrV = 0
        name1 = index_name(t,string) + '.txt'
        
        print(name1)
        arr = read_fiel(name1)
        N_arr = integrate_R(arr,dx)
       
        arrN.append(N_arr)
        arrNx.append(xt)
        xt+=arr2[7]      
        SUM1 = 0
        SUM_v_F = integral(-6,5,arr2[3],arr)
        ARR_v_F.append(SUM_v_F)
        print(SUM_v_F)
        
        #print(len(fp))
        for i in range(0,len(fp),1):                
                fp_N.append(fp[i]*N_arr)
                SUM1+=dx*fp_N[i]
        plt.figure(figsize=(12, 12), constrained_layout=True)
        plt.xticks(np.arange(min(arr1)-dx/2, max(arr1)+dx/2, 0.5))
        plt.yticks(np.arange(min(arr)-0.05, max(arr)+0.05, 0.05))

        plt.plot(arr1,arr,'r',marker='o',label=r"$point \ f_H$")
        
        plt.plot(arr1,fp_N,'y',marker='o',label=r"$f_p= \frac{1}{c_p\sqrt{\pi}}\cdot exp(-\frac{w_x^2}{c_p^2})$")

        #plt.plot(arr1,fH,'g',marker='o',label=r"$f_H=\frac{1}{\sqrt{\pi}}\cdot exp(-(v_x-U_H)^2)$")
        plt.legend(fontsize=14)
        plt.grid ( True )
        plt.tick_params(labelsize=12)
        x = round(t*arr2[7],3)
        name = 'L = ' + str(-arr2[15]) + "; x = " + str(round(x- arr2[15],4))+'; UH =  '+str(round(arr2[16],4))+'; U_p = '+ str(round(arr2[17],4))+ '; Cp = '+str(round(arr2[18],4))  + '; N = ' + str(round(arr2[14],0))         
        plt.xlabel( name, fontsize=20 )
        plt.savefig(index_name(t,string)+'.png')
        plt.grid ( True )
        plt.close()


'''
arrNx1 = []
arrN1 = []
        
plt.figure(figsize=(12, 12), constrained_layout=True)
plt.plot(arrNx,arrN,'r',marker='o',label=r"$n$")
#plt.xticks(np.arange(min(arrNx1), max(arrNx1), 0.5))
#plt.yticks(np.arange(min(arrN1), max(arrN1), 0.01))
plt.grid ( True )
plt.savefig(string +'concentration'+'.png')
plt.show()
'''

#print(ARR_v_F)
print('end')
