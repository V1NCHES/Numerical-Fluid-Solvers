import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import math
import numpy as np
from math import fabs
def method_of_rectangles(c,d,vx,vy,vz, wx,f,Up, u_):
        n = 100
        sum1 = 0
        sum2 = 0
        dx = (d - c)/n
        sum2 = integrate(c,d,dx,vx,vy,vz, wx, n,Up,u_,f)
        while True:
                sum1 = sum2; sum2 = 0
                n *= 2
                dx = (d - c)/n
                sum2 = integrate(c,d,dx,vx,vy,vz, wx, n,Up, u_,f)
                if not (math.fabs(sum1 - sum2) > 1e-5) and (n < 1e5):
                        break
        if (n > 1e5):
                print("ERROR n > N")
        return sum2

def integrate(c,d,dx,vx,vy,vz, wx, n,Up,u_,f):
        integral = 0
        wy = c
        wz = c
        for i in range(n):
                wy = c
                for j in range(n):
                                integral += dx*dx*f(vx,vy,vz,wx,wy,wz,Up, u_)
                                wy+=dx
                wz+=dx	
        return integral

def method_of_rectanglesW(c,d,x, wx,f):
        n = 100
        sum1 = 0
        sum2 = 0
        dx = (d - c)/n
        sum2 = integrateW(c,d,dx,x, wx, n,f)
        while True:
                sum1 = sum2; sum2 = 0
                n *= 2
                dx = (d - c)/n
                sum2 = integrateW(c,d,dx,x, wx, n,f)
                if not (math.fabs(sum1 - sum2) > 1e-5) and (n < 1e5):
                        break
        if (n > 1e5):
                print("ERROR n > N")
        return sum2

def integrateW(c,d,dx,x, wx, n,f):
        integral = 0
        wy = c
        wz = c
        for i in range(n):
                wy = c
                for j in range(n):
                                integral += f(x,wx,wy,wz)
                                wy+=dx
                wz+=dx	
        return dx* dx* integral
def read_fiel(name,arr):
        with open(name,"r") as f:
                reader = f.read().splitlines()
                for reader1 in reader: 
                        for number in reader1.split(" "):
                                #print(number)
                                arr.append(float(number))
        f.close()

def read_fiel_int(name,arr):
        with open(name,"r") as f:
                reader = f.read().splitlines()
                for reader1 in reader: 
                        for number in reader1.split(" "):
                                #print(number)
                                arr.append(int(number))
        f.close()

def G(vx,vy,vz, wx,wy,wz):
        return math.sqrt((vx-wx)**2+(vy-wy)**2+(vz-wz)**2)
def FUH(UH):
        return   math.exp(-UH**2)*math.pi *( 2 + math.exp(UH**2)*math.sqrt(math.pi)* UH*(1 + math.erf(UH) - math.erfc(UH)) )/2

def expVX(x,U,FU):
        return math.pi *math.fabs(x)*math.exp(-(x-U)**2)/ FU
def Arcos(x,U):
        return math.sin(x)/2;
def expVp(vx,vy,vz, wx,wy,wz,U, u_):
        return G(vx,vy,vz, wx,wy,wz)*math.exp(-((wx-U)**2+wy**2 + wz**2))/(u_)
def expVH(x,U):
        return math.exp(-(x-U)**2)/ math.sqrt(math.pi)
def expL(x,U):
        return U*math.exp(-U*x)

def grafic_fp(a,b, dvx,U , f):
        fp = []
        #FU = FUH(U)
        #print('FU = ',FU)
        v = a
        while v < b-dvx/2  :
                fp.append(f(v,U))
                v +=dvx
        #fp.append(f(v,U))
        return fp


#def grafic_fpData(data, f):
#        fp = []
#        v = data[0]
#        while v < data[1]-data[2]/2  :
#                fp.append(f(v,U))
#                v +=data[2]
#        return fp
def grafic_fW(a,b, dvx,c,d, x , f):
        fp = []
        v = a
        while v < b-dvx/2  :
                fp.append(method_of_rectanglesW(c,d,x,v,f))
                v +=dvx
        return fp


def grafic_VxP(a,b,dvx,c,d,Up,u_,vx,vy,vz, f):
        fvxP = []
        wx = a
        while wx < b-dvx/2 :
                fvxP.append(method_of_rectangles(c,d,vx,vy,vz,wx,f,Up,u_))
                wx +=dvx
        return fvxP
def ox(dx,a, b, arr):
        v = a
        while v < b-dx/2  :
                arr.append(v)
                v +=dx  
def grafic_distogramma(namefile,name ,colar, arrV ,arrOX,arr2, label1,label2,label3):
        plt.figure(figsize=(10, 10), constrained_layout=True)
        plt.plot(arrOX,arr2,'r',marker='o',label=label1)
        plt.plot(arrOX,arr2,'b',label=label2)
        plt.xticks(np.arange(min(arrOX), max(arrOX),1.0))
        plt.plot(arrOX,arrV,colar,label=label3)
        plt.legend(fontsize=12 )
        plt.grid ( True )
        plt.xlabel( name, fontsize=20 )
        plt.tick_params(labelsize=14)
        plt.savefig(namefile)
        plt.close()

def data(index_txt, arr1, arr2, A):
        arr =[]
        read_fiel_int(string + str(0)+str(0)+ str(index_txt)+ '.txt', arr)
        read_fiel(string + str(0)+str(0)+ str(index_txt)+ '_.txt',arr1)
        print(len(arr))
        print(arr1)
        n = arr[len(arr)-1]
        ox(arr1[2],arr1[0], arr1[1], A)
        print(len(A))
        for i in range(0,len(arr)-1):
                arr2.append( arr[i] / (n * arr1[2]))
        return n
def normal_extractionem(string):
        print('normal_extractionem')
        arr1 = []       
        arr2 = []
        A = []
        read_fiel_int(string + '.txt', arr)
        read_fiel_int(string + '_.txt', arr1)
        #print(len(arr1),len(arr2),len(A))
        n = data(0, arr1, arr2, A)
        print(len(arr1),len(arr2),len(A), n)

def V(string):        
        arr1 = []       
        arr2 = []
        A = []
        print(len(arr1),len(arr2),len(A))
        n = data(0, arr1, arr2, A)
        print(len(arr1),len(arr2),len(A), n)
        #F = "$\rho(V_x) = \pi*|v_x|exp(-(v_x-U_H)^2)/16.705$"
        name = 'a = ' + str(arr1[0]) + "; b = " + str(arr1[1])+'; d =  '+str(arr1[2])+'; n = '+str(n) +"; U_H = "+str(arr1[3])+'\nCp = '+str(arr1[4])         
        #name = 'a = ' + str(arr1[0]) + "; b = " + str(arr1[1])+'; d =  '+str(arr1[2])+'; n = '+str(n) +"; U_H = "+str(arr1[6])+'\nCp = '+str(arr1[7])         
        grafic_distogramma('check1/'+ str(0)+str(0)+ str(0)+'.png' , name,'g', grafic_fp(arr1[0],arr1[1],arr1[2],arr1[3],expVX),A,arr2,'point Vx','grafic Vx',r"$\rho(V_x) = \pi|V_x|\cdot\frac{exp(-(V_x-U_H)^2)}{f(U_H)}$" )        
        #grafic_distogramma('check1/'+ str(0)+str(0)+ str(3)+'.png' , name,'g', grafic_fp(arr1[0],arr1[1],arr1[2],arr1[6],expVX),A,arr2,'point Vx','grafic Vx',r"$\rho(V_x) = \frac{\pi*|v_x|\cdot exp(-(v_x-U_H)^2)}{16.705}$" )        

def L(string):
        arr1 = []       
        arr2 = []
        A = []
        print(len(arr1),len(arr2),len(A))
        n = data(11, arr1, arr2, A)
        print(len(arr1),len(arr2),len(A), n)
        name = 'a = ' + str(arr1[0]) + "; b = " + str(arr1[1])+'; d =  '+str(arr1[2])+'; n = '+str(n)  + '; Cp = ' + str(arr1[7]) +'\nVx = ' + str(arr1[3])+'; Vy = ' + str(arr1[4])+'; Vz = '+str(arr1[5])            
        print(arr1[6])        
        grafic_distogramma('check1/'+ str(0)+str(0)+ str(11)+'.png',name ,'g', grafic_fp(arr1[0],arr1[1],arr1[2],arr1[6],expL),A,arr2,'point L','grafic L',r"$\rho(L) = \frac{\nu}{|V_H|}exp(-\frac{\nu}{|V_H|}L)$" )

def Vp(string):
        arr1 = []       
        arr2 = []
        A = []
        #print(len(arr1),len(arr2),len(A))
        n = data(22, arr1, arr2, A)
        #print(len(arr1),len(arr2),len(A), n)
        print(arr1)
        name = 'a = ' + str(arr1[0]) + "; b = " + str(arr1[1])+'; d =  '+str(arr1[2])+'; n = '+str(n)  + '; U_H = ' + str(arr1[6]) +'\nVx = ' + str(arr1[3])+'; Vy = ' + str(arr1[4])+"; Vz = "+str(arr1[5])+'; Cp = '+str(arr1[7])+"; x = "+str(arr1[8])           
        #print(arr1[6])
        filename = 'check1/'+ str(0)+str(0)+ str(22)+'.png'
        Grafic_VxP = grafic_VxP(arr1[0],arr1[1],arr1[2],-5,5,arr1[6],math.pi*math.sqrt(math.pi)*arr1[8],arr1[3],arr1[4],arr1[5], expVp)
        grafic_distogramma(filename,name ,'g',Grafic_VxP,A ,arr2,'point Vp','grafic Vp',r"$\rho(\vec{V}) = \frac{|\vec{U}|}{\pi\sqrt{\pi}\cdot u_{*}}\cdot exp( - (V_x-U_P)^2- V_y^2-V_z^2)$")

def VH(string):        
        arr1 = []       
        arr2 = []
        A = []
        print(len(arr1),len(arr2),len(A))
        n = data(3, arr1, arr2, A)
        print(len(arr1),len(arr2),len(A), n)
        name = 'a = ' + str(arr1[0]) + "; b = " + str(arr1[1])+'; d =  '+str(arr1[2])+'; n = '+str(n)  + '; Cp = ' + str(arr1[7]) +'\nVx = ' + str(arr1[3])+'; Vy = ' + str(arr1[4])+'; Vz = '+str(arr1[5])
        grafic_distogramma('check1/'+ str(0)+str(0)+ str(3)+'.png' , name,'g', grafic_fp(arr1[0],arr1[1],arr1[2],arr1[6],expVH),A,arr2,'point Vx','grafic Vx',r"$\rho(V_x) = \frac{ exp(-(v_x-U_H)^2) }{\sqrt{\pi}}$" )        


def Xi(string):
        arr = []       
        arr2 = []
        A = []
        print(len(arr2),len(A))
        arr =[]
        read_fiel_int(string + str(0)+str(0)+ str(4)+ '.txt', arr)
        n = arr[len(arr)-1]
        ox(0.01,0, math.pi, A)
        print(len(A))
        for i in range(0,len(arr)-2):
                arr2.append( arr[i] / (n * 0.01))

        name = 'a = ' + str(0) + "; b = " + str(math.pi)+'; d =  '+str(0.01)+'; n = '+str(n)
        Grafic_fp = grafic_fp(0,math.pi,0.01,0,Arcos)
        print(len(A),len(Grafic_fp))
        x = np.arange(-1 +0.01, 1 - 0.01, 0.01)
        #arrX=[]
        #for i in range(len(x)):
        #        arrX.append(math.acos(x[i])/math.pi)

                
        plt.figure(figsize=(10, 10), constrained_layout=True)
        plt.plot(A,arr2,'r',marker='o',label='point Xi')
        plt.plot(A,arr2,'b',label='grafic Xi')
        plt.xticks(np.arange(min(A), max(A), 0.5))
        plt.plot(A,Grafic_fp,'g',label=r"$Arccos(Xi)$")
        plt.legend(fontsize=12)
        plt.grid ( True )
        plt.xlabel( name, fontsize=20 )
        plt.tick_params(labelsize=14)
        plt.savefig('check1/'+ str(0)+str(0)+ str(4)+'.png')
        plt.close()                
        #plt.plot(arrX, x ,'k',label="Acos(x)" )
        #plt.show()
string = 'check1/'
#V(string)
#L(string)
#Vp(string)
#VH(string)
#Xi(string)
normal_extractionem('check1/chek/000')
print('end')
