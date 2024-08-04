import os
import sys
import cv2
import math
import glob
import scipy
import cmath
import shutil
import numpy as np
from PIL import Image
import plotly.express as px
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from IPython.display import HTML
import plotly.graph_objects as go
from matplotlib import animation, rc
from scipy.special import genlaguerre
from scipy.special import eval_hermite
from plotly.subplots import make_subplots

class Tool:
    
    # calculate complex conjugate
    def ast(self, x):
        return x.conjugate()

    # calculate argument of complex
    def arg(self, x):
        return cmath.phase(x)

    # calculate exponential
    def exp(self, x):
        return np.exp(x)

    # calculate pi
    def pi(self):
        return np.pi
    
    # 3次元上に表示
    def plot_three_dimention(self, x, y, z, min_x, max_x, min_y, max_y, name="title", color="limegreen"):
        fig = go.Figure(data=[go.Surface(z=z, x=x, 
                                         y=y, name=name)])
        fig.update_traces(contours_z=dict(show=True, usecolormap=True, 
                                          highlightcolor="limegreen", project_z=True), showscale=False)
        fig.update_layout(title="", title_font_size=10, 
                  scene=dict(xaxis=dict(title='x', tickfont=dict(size=12)),
                             yaxis=dict(title='t', tickfont=dict(size=12)),
                             zaxis=dict(title='z', tickfont=dict(size=12)),
                             aspectratio=dict(x=1, y=1, z=1)),
                  scene_camera=dict(
                  eye=dict(x=0.5, y=-1.5, z=1.6)),
                  width=800, height=800, 
                  margin=dict(l=65, r=50, b=65, t=90))
        fig.update_layout(
        scene=dict(
            xaxis=dict(nticks=10, range=[min_x, max_x]),yaxis=dict(nticks=10,
            range=[min_y, max_y]),))
        fig.show()

    def plot_contour(self, x, y, z, min_x, max_x, min_y, max_y):
        fig = plt.figure(figsize = (8, 8))
        ax = fig.add_subplot(111)
        ax.set_title("", size = 15)
        ax.grid()
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)
        ax.set_xlabel("x", size = 14)
        ax.set_ylabel("t", size = 14)
        contour = ax.contour(x, y, abs(z)**2,
                      levels=20, cmap='viridis')
        ax.clabel(contour, inline=True, fontsize=8)
        ax.legend()
        plt.show()

    # ２次元上に表示
    def plot_matrix(self, x, x_range, y_range, min_x, max_x, min_t, max_t):
        plt.figure(figsize=(7, 7))
        plt.xlabel('x')
        plt.ylabel('t')
        plt.xticks(ticks=np.linspace(0, x_range, 9), labels=np.linspace(min_x, max_x, 9))
        plt.yticks(ticks=np.linspace(0, y_range, 5), labels=np.linspace(max_t, min_t, 5))
        plt.gca().invert_yaxis()
        plt.imshow(x)
        plt.colorbar()
        plt.show()   

class Bright_soliton:

    def __init__(self, x, t, a, b, c, k, alpha):
        self.x = x
        self.t = t
        self.a = a
        self.b = b
        self.c = c
        self.k = k
        self.alpha = alpha
        self.tool = Tool()

    def b_sigma(self):
        return self.a*self.c-abs(self.b)**2
    
    def b_theta(self, j):
        return self.k[j]*self.x+1j*self.k[j]**2*self.t
    
    def b_r(self, i, j):
        Sigma = self.b_sigma()
        return (self.a*self.alpha[1][i]*self.tool.ast(self.alpha[1][j]) \
                +self.a*Sigma*self.alpha[2][i]* \
                self.tool.ast(self.alpha[2][j]))/(self.k[i]+self.tool.ast(self.k[j]))**2

    def b_g1(self, j):
        theta1 = self.b_theta(1)
        theta2 = self.b_theta(2)
        return self.alpha[j][1]*self.tool.exp(theta1) \
            +self.alpha[j][2]*self.tool.exp(theta2)
    
    def b_g3(self, j):
        theta1 = self.b_theta(1)
        theta2 = self.b_theta(2)        
        return self.b_beta(j, 1)*self.tool.exp(theta1+self.tool.ast(theta1)+theta2) \
                +self.b_beta(j, 2)*self.tool.exp(theta2+self.tool.ast(theta2)+theta1)
    
    def b_f2(self):
        theta1 = self.b_theta(1)
        theta2 = self.b_theta(2)
        return self.b_r(1, 1)*self.tool.exp(theta1+self.tool.ast(theta1)) \
                +self.b_r(1, 2)*self.tool.exp(theta1+self.tool.ast(theta2)) \
                +self.b_r(2, 1)*self.tool.exp(self.tool.ast(theta1)+theta2) \
                +self.b_r(2, 2)*self.tool.exp(theta2+self.tool.ast(theta2))
    
    def b_f4(self):
        theta1 = self.b_theta(1)
        theta2 = self.b_theta(2)
        mu = ((self.tool.ast(self.k[2])-self.tool.ast(self.k[1]))*(self.b_r(1, 2)*self.b_beta(1, 1) \
                *(self.tool.ast(self.k[1])+self.k[2])-self.b_r(1, 1)*self.b_beta(1, 2) \
                *(self.k[2]*self.tool.ast(self.k[2])))) \
                /(self.alpha[1][1]*(self.k[2]+self.tool.ast(self.k[2])) \
                *(self.tool.ast(self.k[1])+self.k[2]))
        return mu*self.tool.exp(theta1+self.tool.ast(theta1)+theta2+self.tool.ast(theta2))
    
    def b_beta(self, i, j):
        return (self.k[j]-self.k[3-j])*((self.b_r(3-j, j)*self.alpha[i][j]) \
                /(self.k[j]+self.tool.ast(self.k[j]))-(self.b_r(j, j) \
                *self.alpha[i][3-j])/(self.k[3-j]+self.tool.ast(self.k[j])))
    
class Dark_soliton:

    def __init__(self, x, t, a, b, c, k, rho_1, rho_2, omega, theta_0, s_1, s_2, xi_0):
        self.x = x
        self.t = t
        self.a = a
        self.b = b
        self.c = c
        self.k = k
        self.rho_1 = rho_1
        self.rho_2 = rho_2
        self.rho = [0, rho_1, rho_2]
        self.omega = omega
        self.theta_0 = theta_0
        self.s = [0, s_1, s_2]
        self.xi_0 = xi_0
        self.tool = Tool()

    def d_sigma(self):
        return self.a*self.c-abs(self.b)**2
    
    def d_lamuda(self):
        return -2*self.a*self.rho_1**2-2*self.a*self.d_sigma()    
    
    def d_theta(self):
        return self.k*self.x+self.omega*self.t+self.theta_0
    
    def d_xi(self, j):
        return self.s[j]*self.x-(self.s[j]**2+self.d_lamuda())*self.t+self.xi_0[j]
    
    def d_g0(self, j):
        return self.rho[j]*self.tool.exp(1j*self.d_xi(j))

    def d_B(self, j):
        return (-self.k**2+2*1j*self.k*self.s[j]+1j*self.omega)/(self.k**2+2*1j*self.k*self.s[j]+1j*self.omega)

    def d_g1(self, j):
        return self.d_B(j)*self.tool.exp(self.d_theta())
    
    # def d_r(self, i, j):
    #     Sigma = self.d_sigma()
    #     return (self.a*self.alpha[1][i]*self.tool.ast(self.alpha[1][j]) \
    #             +self.a*Sigma*self.alpha[2][i]* \
    #             self.tool.ast(self.alpha[2][j]))/(self.k[i]+self.tool.ast(self.k[j]))**2
    
    # def d_g3(self, j):
    #     theta1 = self.d_theta(1)
    #     theta2 = self.d_theta(2)        
    #     return self.d_beta(j, 1)*self.tool.exp(theta1+self.tool.ast(theta1)+theta2) \
    #             +self.d_beta(j, 2)*self.tool.exp(theta2+self.tool.ast(theta2)+theta1)
    
    def d_f1(self):
        return self.tool.exp(self.d_theta())
    
    # def d_f2(self):
    #     theta1 = self.d_theta(1)
    #     theta2 = self.d_theta(2)
    #     return self.d_r(1, 1)*self.tool.exp(theta1+self.tool.ast(theta1)) \
    #             +self.d_r(1, 2)*self.tool.exp(theta1+self.tool.ast(theta2)) \
    #             +self.d_r(2, 1)*self.tool.exp(self.tool.ast(theta1)+theta2) \
    #             +self.d_r(2, 2)*self.tool.exp(theta2+self.tool.ast(theta2))
    
    # def d_f4(self):
    #     theta1 = self.d_theta(1)
    #     theta2 = self.d_theta(2)
    #     mu = ((self.tool.ast(self.k[2])-self.tool.ast(self.k[1]))*(self.d_r(1, 2)*self.d_beta(1, 1) \
    #             *(self.tool.ast(self.k[1])+self.k[2])-self.d_r(1, 1)*self.d_beta(1, 2) \
    #             *(self.k[2]*self.tool.ast(self.k[2])))) \
    #             /(self.alpha[1][1]*(self.k[2]+self.tool.ast(self.k[2])) \
    #             *(self.tool.ast(self.k[1])+self.k[2]))
    #     return mu*self.tool.exp(theta1+self.tool.ast(theta1)+theta2+self.tool.ast(theta2))
    
    # def d_beta(self, i, j):
    #     return (self.k[j]-self.k[3-j])*((self.d_r(3-j, j)*self.alpha[i][j]) \
    #             /(self.k[j]+self.tool.ast(self.k[j]))-(self.d_r(j, j) \
    #             *self.alpha[i][3-j])/(self.k[3-j]+self.tool.ast(self.k[j])))