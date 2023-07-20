from numpy import *
from math import *
import pandas as pd
import matplotlib.pyplot as plt 
from .SuspecFunction import Chi1xx,Chi1xy,Chi2xx,Chi2xy



############
class Susceptibility_Lorentz_Fit:
    def __init__(self, lorentz_path, paramp, paramm):
        self.Lorentz = pd.read_csv(lorentz_path, sep=' ')
        self.paramp = paramp
        self.paramm = paramm
        self.xrange = self.Lorentz['x']#arange(-30,500)

#複素感受率の虚部を計算しdataframeに格納
    def calc_susceptibility(self): 
        Chi1p = array([(Chi1xx(2*pi*f*10**9,*self.paramp)+ 1j*Chi1xy(2*pi*f*10**9,*self.paramp)).imag/(f*10**9) for f in self.xrange])  #
        Chi1m = array([(Chi1xx(2*pi*f*10**9,*self.paramm)- 1j*Chi1xy(2*pi*f*10**9,*self.paramm)).imag/(f*10**9) for f in self.xrange])  #
        Chi2p = array([(Chi2xx(2*pi*f*10**9,*self.paramp)+ 1j*Chi2xy(2*pi*f*10**9,*self.paramp)).imag/(f*10**9) for f in self.xrange])  #
        Chi2m = array([(Chi2xx(2*pi*f*10**9,*self.paramm)- 1j*Chi2xy(2*pi*f*10**9,*self.paramm)).imag/(f*10**9) for f in self.xrange])  #
        data = {'x':self.xrange,'1p':Chi1p,'1m':Chi1m, '2p':Chi2p, '2m':Chi2m, 'p':Chi1p+Chi2p, 'm':Chi1m+Chi2m} #副格子1の右回り+と左回り-の磁気感受率の虚部
        self.sus_imag = pd.DataFrame(data)

#複素感受率の虚部をプロット
    def plot_susceptibility(self):
        fig1, ax1 = plt.subplots()
        ax1.set_title('susceptibility')
        ax1.plot(self.sus_imag['x'], self.sus_imag['p'],label = 'p')
        ax1.plot(self.sus_imag['x'], self.sus_imag['m'],label = 'm')
        ax1.legend()
        plt.show()

#ローレンツ関数と複素感受率の虚部のピークの大きさをあわせる
    def adjust_peak(self):
        if self.sus_imag['p'].max() >= 0: #副格子1の右回り+の磁気感受率の虚部の最大値が正かどうかで場合分け
            self.peak_ratio1 = self.Lorentz['L1'].max()/abs(self.sus_imag['p'].max())
        else:
            self.peak_ratio1 = (-1)*self.Lorentz['L1'].max()/abs(self.sus_imag['p'].min())
        
        if self.sus_imag['m'].min() >= 0: #左回りについても同様
            self.peak_ratio2 = self.Lorentz['L2'].max()/abs(self.sus_imag['m'].max())
        else:
            self.peak_ratio2 = (-1)*self.Lorentz['L2'].max()/abs(self.sus_imag['m'].min())

        self.sus_imag['p'] = self.sus_imag['p'] * self.peak_ratio1 
        self.sus_imag['m'] = self.sus_imag['m'] * self.peak_ratio2


#複素感受率とローレンツ関数を重ねてプロット
    def plot_sus_Lorentz(self):
        fig2, ax2 = plt.subplots()
        ax2.set_title('each peak')
        ax2.plot(self.sus_imag['x'],self.sus_imag['p'],label='p')
        ax2.plot(self.sus_imag['x'],self.sus_imag['m'],label='m')
        ax2.plot(self.Lorentz['x'],self.Lorentz['L1'],label='L1')
        ax2.plot(self.Lorentz['x'],self.Lorentz['L2'],label='L2')
        ax2.legend()

        fig3, ax3 = plt.subplots()
        ax3.set_title('sum')
        ax3.plot(self.sus_imag['x'],self.sus_imag['p']+self.sus_imag['m'],label='p+m')
        ax3.plot(self.Lorentz['x'],self.Lorentz['Lsum'],label='Lorentz sum')
        ax3.legend()