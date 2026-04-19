import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import math
import numpy as np
from math import fabs

def read_fiel(name):
        arr =[]
        with open(name,"r") as f:
                reader = f.read().splitlines()
                for reader1 in reader: 
                        for number in reader1.split(" "):
                                #print(number)
                                arr.append(float(number))
        f.close()
        return arr

def read_fiel_int(name):
        arr =[]
        with open(name,"r") as f:
                reader = f.read().splitlines()
                for reader1 in reader: 
                        for number in reader1.split(" "):
                                #print(number)
                                arr.append(int(number))
        f.close()
        return arr
def index_name(index,string):
        name = string# +str(int(i/100))+str(int(i/10))+str(0)+'.txt'
        if index < 10:
                name += str(0)+str(0)+str(index)
        elif index < 100 and index > 9:
                name += str(0)+str(int(index / 10))
                index = index -  int(index / 10) * 10
                name += str(index)
        elif index > 100 and index < 1000:
                name += str(index / 100)
                index = index - int(index / 100) * 100
                name += str(index/10)
                index = index -int(index / 10) * 10
                name += str(index)
        return name
def expx(x,U):
        return math.fabs(x)*math.exp(-(x-U)**2)
def expVx(x,y,z,U):
        return math.fabs(x)*math.exp(-(x-U)**2 - y**2 - z**2)/16.705
def expVX(x,U):
        return math.pi *math.fabs(x)*math.exp(-(x-U)**2)/16.705
def expy(y,U):
        return math.exp(-y**2)
def expz(z,U):
        return math.exp(-z**2)

def expL(x,U):
        return math.exp(-x)

def grafic_fp(a,b, dvx,U , f):
        fp = []
        v = a
        while v < b-dvx/2  :
                fp.append(f(v,U))
                v +=dvx
        return fp

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
                if not (math.fabs(sum1 - sum2) > 1e-5) and (n < 1e5):
                        break
        if (n > 1e4):
                print("ERROR n > N")
        return sum2

def grafic_Vx(a,b, dvx, c,d,cp , f):
        fvx = []
        v = a
        while v < b-dvx/2 :
                fvx.append(method_of_rectangles(c,d, v,f,cp))
                v +=dvx
        return fvx
def ox(dx,x, b):
        arr = []
        v = x
        while v < b-dx/2  :
                arr.append(v)
                v +=dx
        return arr
def distogramma(name, namefile , arr, n):
        fig = plt.figure(figsize = (6,4))
        ax = fig.add_subplot()
        ax.hist(arr,n)
        ax.grid()
        plt.xlabel( name )
        plt.savefig(namefile)
        #plt.show()
        plt.close()

def grafic_distogramma(namefile ,colar, arrV ,arrOX ):
        plt.plot(arrOX,arrV,colar)
        plt.grid ( True )
        plt.savefig(namefile)
        plt.close()       
def check_velocity(string):      
        name2 = string + str(0)+str(0)+ str(0)+ '.txt'
        arrx =[]
        arry =[]
        arrz =[]
        arr =[]
        arrOX = []
        arrV = []
        arr = read_fiel(name2)
        print(len(arr))
        for i in range(0,len(arr),3):
                arrx.append(arr[i])
                arry.append(arr[i+1])
                arrz.append(arr[i+2])
        print(len(arrx))
        arr =[]
        arr = read_fiel('data/dat2/'+ str(0) + '.txt')
        #arrvx = grafic_fp(arr[0],arr[1],arr[3],arr[16],expx)

        distogramma('vx.n = ' + str(int(arr[2]))+"; N = " + str(len(arrx)), 'check/'+ str(0)+str(0)+ str(0)+'.png' , arrx, int(arr[2]))
        grafic_distogramma('check/'+ str(0)+str(0)+ str(0)+'_1.png' ,'g', grafic_fp(arr[0],arr[1]-2,arr[3],arr[16],expx) ,ox(arr[3],arr[0], arr[1]-2))
       # grafic_distogramma('check/'+ str(0)+str(0)+ str(0)+'_.png' ,'g', grafic_fp(arr[0],arr[1],arr[3],arr[16],expx) ,ox(arr[3],arr[0], arr[1]))

        distogramma('vy.n = ' + str(int(arr[2]))+"; N = " + str(len(arrx)), 'check/'+ str(0)+str(0)+ str(1)+'.png' , arry, int(arr[2]))
        grafic_distogramma('check/'+ str(0)+str(0)+ str(1)+'_1.png' ,'r', grafic_fp(arr[8],arr[9]-2,arr[3],arr[16],expy) ,ox(arr[3],arr[8], arr[9]-2))
        #grafic_distogramma('check/'+ str(0)+str(0)+ str(1)+'_.png' ,'r', grafic_fp(arr[8],arr[9],arr[3],arr[16],expy) ,ox(arr[3],arr[8], arr[9]))

        distogramma('vz.n = ' + str(int(arr[2]))+"; N = " + str(len(arrx)), 'check/'+ str(0)+str(0)+ str(2)+'.png' , arrz, int(arr[2]))
        grafic_distogramma('check/'+ str(0)+str(0)+ str(2)+'_1.png' ,'b', grafic_fp(arr[11],arr[12]-2,arr[3],arr[16],expz) , ox(arr[3],arr[11], arr[12]-2))
        #grafic_distogramma('check/'+ str(0)+str(0)+ str(2)+'_.png' ,'b', grafic_fp(arr[11],arr[12],arr[3],arr[16],expz) , ox(arr[3],arr[11], arr[12]))

        print(arr[16],arr[0],arr[1])
        
def check_leght_velocity(string):
        arrx =[]
        arry =[]
        arrz =[]
        arr =[]
        name2 = string + str(0)+str(0)+ str(2)+ '.txt'
        arr = read_fiel(name2)
        print(len(arr))
        for i in range(0,len(arr),3):
                arrx.append(arr[i])
                arry.append(arr[i+1])
                arrz.append(arr[i+2])
        print(len(arrx))
        arr =[]
        N = 100
        distogramma('Wx + N = ' + str(N), 'check/'+ str(0)+str(0)+ str(4)+'_.png' , arrx, N)
        distogramma('Wy + N = ' + str(N), 'check/'+ str(0)+str(0)+ str(5)+'_.png' , arry, N)
        distogramma('Wz + N = ' + str(N), 'check/'+ str(0)+str(0)+ str(6)+'_.png' , arrz, N)
       
        name2 = string + str(0)+str(0)+ str(1)+ '.txt'
        arr = read_fiel(name2)
        N = 10000
        distogramma('Leght + N = ' + str(N), 'check/'+ str(0)+str(0)+ str(3)+'_.png' , arr, 10000)
        grafic_distogramma('check/'+ str(0)+str(0)+ str(4)+'_1.png' ,'g', grafic_fp(0,6,0.01,1,expL) ,ox(0.01,0, 6))

def V(string):
        name2 = string + str(0)+str(0)+ str(0)+ '.txt'
        arr =[]
        arr1 = []
        arr = read_fiel_int(name2)
        arr1 = read_fiel('data/dat5/'+ str(0) + '.txt')
        print(len(arr))
        arr2 = []
        print(arr1)
        A = ox(arr1[3],arr1[0], arr1[1])
        print(len(A))
        for i in range(0,len(arr)):
                arr2.append( arr[i] / (arr1[14] * arr1[3]))
        #print(arr2)
        plt.plot(A,arr2,'ro',A,arr2)
        plt.grid ( True )
        #plt.savefig(name2+".png")
        #plt.close()
        
       # grafic_distogramma('check1/'+ str(0)+str(0)+ str(0)+'_1.png' ,'g', grafic_Vx(arr1[0],arr1[1],arr1[3],arr1[8],arr1[9],arr1[16],expVx),A )
        grafic_distogramma('check1/'+ str(0)+str(0)+ str(0)+'_1.png' ,'g', grafic_fp(arr1[0],arr1[1],arr1[3],arr1[16],expVX),A )
        plt.close()
        #plt.savefig(name2+'.png')
        
string = 'check1/'
V(string)




#check_velocity(string)
#check_leght_velocity(string)
#arrXi =[]
#arreps =[]
#name = string +'xi.txt'
#arrXi = read_fiel(name)
#name = string +'eps.txt'
#arreps = read_fiel(name)

#distogramma('Xi', 'check/xi.png' , arrXi, 1000)
#distogramma('eps', 'check/eps.png' , arreps, 1000)
print('end')
