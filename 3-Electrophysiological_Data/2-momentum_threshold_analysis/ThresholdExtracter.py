


import numpy as np
import pandas as pd
import os
from scipy import interpolate 
from sklearn.utils import resample

from matplotlib.patches import Rectangle

from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import functools
import seaborn as sns
from tqdm import tqdm
from scipy import stats as scpstats
import random
from scipy.signal import savgol_filter
from shapely.geometry import LineString
from matplotlib.ticker import FixedLocator, FixedFormatter,AutoMinorLocator
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm

def intersect(line1,line2):    
    firstline=LineString(line1)
    secondline = LineString(line2)   
    try:
        intersectionA = firstline.intersection(secondline)
        vth=intersectionA.xy[1][0]
    except:
        vth=np.nan
    return(vth)

def df_to_numpy(df):
    df=df.astype(float)
    df=df.T
    df=df.to_numpy()
    return(df)

def FiringRateBis(spiketrain, inf, sup):
    xx = np.arange(inf, sup, 0.001)
    yy = np.zeros_like(xx)
    for spike in spiketrain:
        yy += scpstats.norm.pdf(xx, loc=spike, scale=0.05) #bandwidth de 50 ms
    return(xx,yy)

Vitesses_dico ={}
Vitesses_dico["R7.5_v5"]  =pd.read_csv("Volumes/2_Volume_R7.5_V5.csv",header=None)
Vitesses_dico["R7.5_v10"] =pd.read_csv("Volumes/2_Volume_R7.5_V10.csv",header=None)
Vitesses_dico["R7.5_v20"] =pd.read_csv("Volumes/2_Volume_R7.5_V20.csv",header=None)
Vitesses_dico["R7.5_v40"] =pd.read_csv("Volumes/2_Volume_R7.5_V40.csv",header=None)
Vitesses_dico["R15_v5"] =pd.read_csv("Volumes/2_Volume_R15_V5.csv",header=None)
Vitesses_dico["R15_v10"]=pd.read_csv("Volumes/2_Volume_R15_V10.csv",header=None)
Vitesses_dico["R15_v20"]=pd.read_csv("Volumes/2_Volume_R15_V20.csv",header=None)
Vitesses_dico["R15_v40"]=pd.read_csv("Volumes/2_Volume_R15_V40.csv",header=None)
Vitesses_dico["R25_v5"] =pd.read_csv("Volumes/2_Volume_R25_V5.csv",header=None)
Vitesses_dico["R25_v10"]=pd.read_csv("Volumes/2_Volume_R25_V10.csv",header=None)
Vitesses_dico["R25_v20"]=pd.read_csv("Volumes/2_Volume_R25_v20.csv",header=None)
Vitesses_dico["R25_v40"]=pd.read_csv("Volumes/2_Volume_R25_v40.csv",header=None)

def return_velo(V,r,interp=False):
    NTIMES=101
    string= "R"+ str(r)+"_v"+str(V)
    T=Vitesses_dico[string] #Rechercher clé dans le dico
    X0=np.linspace(0,6/V,NTIMES) #axs[j].plot(np.linspace(0,6/V,NTIMES),T[:,NPOSES//2],label=str(V))
    X=X0-4/V #mise à zéro
    idx=np.min(np.argwhere(X > 0))
    
    if interp==True:
        tCFD,fCFD=X[:idx+1],T.loc[:idx]
        indices = np.logical_not(np.isnan(np.asarray(fCFD))).flatten()
        tCFD = np.asarray(tCFD)[indices]
        fCFD = np.asarray(fCFD)[indices].squeeze()
        temp = interpolate.interp1d(tCFD, fCFD, fill_value="extrapolate") 
        xnew = np.arange(np.nanmin(X[:idx+1]), 0,0.0001)
        ynew = temp(xnew) 
        return(xnew,ynew)
    
    else:
        return(X[:idx+1],T.loc[:idx] )
    
    
clusters=pd.read_csv('ACF_spikes_copy_2.csv')

path_to_folder = 'cricket'

spikes_data_root = f'{path_to_folder}/spikes/'
spikes_files = [spikes_data_root+x for x in os.listdir(spikes_data_root) if x[-3:]=='txt']

valves_data_root = f'{path_to_folder}/stims/'
valves_files = [valves_data_root+x for x in os.listdir(valves_data_root) if x[-3:]=='txt']

codes_data_root = f'{path_to_folder}/codes/'
codes_files = [codes_data_root+x for x in os.listdir(codes_data_root) if x[-4:]=='xlsx']

points_data_root = f'{path_to_folder}/points/'
points_files = [points_data_root+x for x in os.listdir(points_data_root) if x[-3:]=='txt']


def return_valves():
    return(valves_files)


class SplitFile:
    def __init__(self,index):
        self.index = index #a voir
        self.root=os.getcwd()
        os.chdir(self.root)
        print(valves_files[self.index])
        self.stimulus = pd.read_csv(valves_files[self.index], sep='\t')
        self.points = pd.read_csv(points_files[self.index], sep='\t',header=None)
        self.codes = pd.read_excel(codes_files[self.index],header=0)
        self.name=valves_files[self.index][14:-11] #a verifier
        try:
            self.cluster=clusters.loc[clusters['file']==self.name,'ID'].iloc[0]
        except:
            self.cluster=-1
        
    def retrieve_times(self):
        i=0
        control_S=[]
        control_E=[]
        looming_S=[]
        looming_E=[]
        
        #self.codes.drop(self.codes[self.codes['keep'] !=1].index, inplace = True)
        #self.codes=self.codes.reset_index()
        duration=len(self.codes.index)
        for k in range(duration):
            #print(self.codes.loc[k])
            if self.codes['keep'][k]==1:
                if self.codes['stim'][k]=="R":
                    control_S+=[self.points[0][i]]
                    control_E+=[self.points[0][i+1]]
                    i+=2

                elif self.codes['stim'][k]=="L":
                    looming_S+=[self.points[0][i]] 
                    looming_E+=[self.points[0][i+1]]
                    i+=2 #si 0.2 ça marche car ligne verticale pour fin stimulus

                else:
                    print("error")
            else:
                i+=2
                
        self.CS=control_S
        self.CE=control_E
        self.SS=looming_S
        self.SE=looming_E
        
        A=pd.read_csv(spikes_files[self.index], sep='\t')
        spikestimes=A.loc[2:]
        SPK=df_to_numpy(spikestimes)
        self.SPK=SPK
        
        
    def extract_valves(self, plot=False, single=True):
        self.tab=pd.DataFrame(columns=['file','v','R','maxFR','timingFR','spikes','stimulus','repet']) #,'FR_pre'

        os.chdir(self.root)
        A=pd.read_csv(spikes_files[self.index], sep='\t')
        spikestimes=A.loc[2:]
        SPK=df_to_numpy(spikestimes)
        self.SPK=SPK
        self.FRs={}
        
        
        self.LINES=[]
        self.LIST_VTH=[]
        self.LIST_TPEAK=[]
        
        #ici pas dans retrieve times
        self.codes.drop(self.codes[self.codes['keep'] !=1].index, inplace = True)
        self.codes=self.codes.reset_index()
        duration=len(self.codes.index)
        cs=0
        cl=0    
        control_S=self.CS
        control_E=self.CE
        looming_S=self.SS
        looming_E=self.SE
        
        lowlim,upplim=0,0.2 
        NAME=str(codes_files[self.index][14:-10])
        
        subfolder1='indivspik/'
        subfolder='indivspik/'+NAME
        os.chdir(f'{self.root}/{subfolder1}/')
    
        if os.path.exists(f'{self.root}/{subfolder}/'):
            os.chdir(f'{self.root}/{subfolder}/')
        else:
            os.mkdir(f'{self.root}/{subfolder}/')
            os.chdir(f'{self.root}/{subfolder}/')
        
        if plot==True:
            pp = PdfPages('single_plots'+NAME+'.pdf')
        
        
        
        for k in range(duration):
            if k%2==0:
                try:
                    plt.show()
                    pp.savefig(fig)
                except:
                    plt.show()
                    
                if single:
                    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,4))

            if self.codes['stim'][k]=="L":

                duration=looming_E[cs]-looming_S[cs]
                
                SPKcopy=SPK[SPK<looming_E[cs]+upplim]
                SPKcopy=SPKcopy[SPKcopy>looming_S[cs]-lowlim]-(looming_E[cs]) #-(looming_S[cs]-lowlim)
                
                filterstim=self.stimulus.loc[(self.stimulus['Time'] >= looming_S[cs]-lowlim) & (self.stimulus['Time'] <= looming_E[cs]+upplim)]
                
                C=np.asarray((filterstim['Time']-(looming_S[cs]+duration),filterstim["3 IN 4     "]/20))  
                
                
                #B=FiringRateBis(SPKcopy,0, looming_E[cs]-looming_S[cs]+upplim+lowlim)
                B=FiringRateBis(SPKcopy,-lowlim-duration, upplim)
                self.FRs[int(self.codes['v'][k])] = B
                #self.SPKcopy=SPKcopy
                
                if int(self.codes['v'][k]) in [5,10,20,40] and int(self.codes['D'][k]) in [75,150,250]:
                    if int(self.codes['D'][k])==75:
                        tCFD,fCFD=return_velo(int(self.codes['v'][k]),int(self.codes['D'][k])/10)
                    else:
                        tCFD,fCFD=return_velo(int(self.codes['v'][k]),int(self.codes['D'][k]/10))
                    line1=np.column_stack((tCFD[:len(fCFD.dropna())],fCFD.dropna()))
                    maxFR_tim=B[0][B[1].argmax()] if max(B[1])!=0 else 0
                    self.B=B
                    line2=np.column_stack((np.repeat(maxFR_tim,len(fCFD.dropna())),np.linspace(0,1e-7,len(fCFD.dropna()))))

                    vth=intersect(line1,line2)
                    
                    self.LINES+=[line1]
                    self.LIST_VTH+=[vth]
                    self.LIST_TPEAK+=[maxFR_tim]
                    
                if single:
                    ax1.eventplot(SPKcopy,color='grey',alpha=0.05)
                    ax1.plot(filterstim['Time']-(looming_S[cs]+duration),filterstim["3 IN 4     "]/20,linewidth=3,color='black')
                    ax1.set_title('LOOM_v'+str(self.codes['v'][k])+"_D"+str(int(self.codes['D'][k])))  
                
                    try:
                        ax1.plot(B[0],B[1]/max(B[1]),linestyle='dashed',alpha=0.5,color='grey')
                        
                    except:
                        ax1.plot(B[0],B[1],linestyle='dashed',alpha=0.5,color='grey')
                
                new_row = {'file': NAME[:-1],
                           'R':int(self.codes['D'][k]),
                           'v': self.codes['v'][k],
                           'maxFR': max(B[1]),
                           'timingFR':B[0][B[1].argmax()] if max(B[1])!=0 else 0,
                           'spikes':SPKcopy,
                           'stimulus':'LOOM',
                           'repet':0
                          }
                self.tab=pd.concat([self.tab, pd.DataFrame([new_row])], ignore_index=True)
                np.savetxt(NAME+'LOOM_v'+str(self.codes['v'][k])+"_D"+str(int(self.codes['D'][k]))+'_'+str(k)+'.txt',C)
                cs+=1
                
                
            elif self.codes['stim'][k]=="R":
                
                duration=control_E[cl]-control_S[cl]
                SPKcopy=SPK[SPK<control_E[cl]+upplim]
                
                SPKcopy=SPKcopy[SPKcopy>control_S[cl]-lowlim] -(control_E[cl])
                
                filterstim=self.stimulus.loc[(self.stimulus['Time'] >= control_S[cl]-lowlim) & (self.stimulus['Time'] <= control_E[cl]+upplim)]
                
                
                C=np.asarray((filterstim['Time']-(control_S[cl]+duration),filterstim["3 IN 4     "]*2.5/20))  
                
                B=FiringRateBis(SPKcopy,-lowlim-duration, upplim)
                
                
                if single:
                    ax2.eventplot(SPKcopy,color='grey',alpha=0.05)
                    ax2.plot(filterstim['Time']-(control_S[cl]+duration),filterstim["3 IN 4     "]*2.5/20,linewidth=3,color='black')
                    ax2.set_title('REC_v'+str(self.codes['v'][k])+"_D"+str(int(self.codes['D'][k])))
                    try:
                        ax2.plot(B[0],B[1]/max(B[1]),linestyle='dashed',alpha=0.5,color='grey')
                    except:
                        ax2.plot(B[0],B[1],linestyle='dashed',alpha=0.5,color='grey')
                    np.savetxt(NAME+'REC_v'+str(self.codes['v'][k])+"_D"+str(int(self.codes['D'][k]))+'_'+str(k-1)+'.txt',C)
                cl+=1
            else:
                print('error')
        os.chdir(self.root)
        
        if plot==True:
            fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,4))
            plot_neuron('timingFR',S,ax1)
            plot_neuron('maxFR',S,ax2)
            pp.savefig(fig)
            pp.close()
        plt.show()
        
        
        
def plot_neuron(metric,S,ax):
    sns.lineplot(x='v', y=metric, data=S.tab[S.tab['R']==75],ax=ax, marker='o', color='tab:blue')  #lienplot hue="file",alpha=alpha, 
    sns.lineplot(x='v', y=metric, data=S.tab[S.tab['R']==150],ax=ax, marker='o', color='tab:orange')  #lienplot hue="file",alpha=alpha, 
    sns.lineplot(x='v', y=metric, data=S.tab[S.tab['R']==250],ax=ax, marker='o', color='tab:green')  #lienplot hue="file",alpha=alpha, 


def find_threshold_v1(S,ax): #filtrage

    T=pd.DataFrame(columns=['time','m','b','S','pval','R2','adjR','mu','NPTS','TPTS','residuals','keeppoints','var'])
    for delta in np.linspace(0,0.2,500): #200ms
        
        
        Y=[]
        X=[]
        for k in range(len(S.LINES)):
            l1=S.LINES[k]
            l2=np.column_stack((np.repeat(S.LIST_TPEAK[k]-delta,len(l1)),np.linspace(0,0.401e-7,len(l1))))
            vth=intersect(l1,l2)

            X+=[S.LIST_TPEAK[k]-delta]
            Y+=[vth]
        lentotalpoints=len(X)
        indices = np.logical_not(np.logical_or(np.isnan(X), np.isnan(Y)))
        X = np.array(X)[indices]
        Y = np.array(Y)[indices]
        
        try:
            X2 = sm.add_constant(X)
            model = sm.OLS(Y, X2)
            results = model.fit()
            influence = results.get_influence()
            standardized_residuals = abs(influence.resid_studentized_internal)
        except:
            standardized_residuals = np.nan

        try:
            b,m=results.params #attention ordre inversé
            pval=results.pvalues[1]
            adjR=results.rsquared_adj
            R2=results.rsquared
        except: 
            b,m=np.nan,np.nan #attention ordre inversé
            pval=np.nan
            adjR=np.nan
            R2=np.nan
        mu=np.mean(Y)
        SUM=0
        for vth in Y:
            SUM+= (vth - mu)**2
        
        if len(X)>0:
            SUM=np.sqrt(SUM/len(X))      #RMSE
        else:
            SUM=np.nan
            print('error')

        new_row = {'time':delta,'m':m,'b':b,'S':SUM,'pval':pval,'R2':R2,'adjR':adjR,'mu':mu,'NPTS':len(X),'TPTS':lentotalpoints, 'residuals':standardized_residuals, 'keeppoints': np.sum(standardized_residuals<3), 'var':np.sqrt(np.var(Y))}
        T=pd.concat([T, pd.DataFrame([new_row])], ignore_index=True)

    try:
        TBIS=T.loc[T['NPTS']==max(T['NPTS'])]
        K=T.iloc[abs(TBIS['S']).idxmin()]
    except:
        K=T.iloc[abs(T['S']).idxmin()]
        print('pbm')
        """
        K={'time':np.nan,
           'm':np.nan,
           'b':np.nan,
           'S':np.nan,
           'pval':np.nan,
           'R2':np.nan,
           'adjR':np.nan,
           'mu':np.nan,
          'NPTS':len(X),
          'TPTS':lentotalpoints,
          'residuals':np.nan,
          'keeppoints':np.nan,
          'var':np.nan}
        """
    #print('residuals',K['residuals'],K['TPTS'], K['keeppoints'])
    #print('other',K['var'],K['S'])
    LABELS=S.codes.loc[(S.codes['keep']==1) & (S.codes['stim']=='L')]
    LABELS=LABELS[LABELS['v'].isin([5, 10,20,40])]
    LENGTH=np.arange(0,len(S.LINES),1)
    VTH=[]
    for k,j,l in zip(LENGTH,LABELS['v'],LABELS['D']):
        
        if l==250:
            color='tab:green'
        elif l==150:
            color='tab:orange'
        else:
            color='tab:blue'
            
            
        if j==40:
            marker="$40$"
        elif j==20:
            marker="$20$"
        elif j==10:
            marker="$10$"
        else:
            marker="$05$"

        
        #ax.plot(S.LINES[k][:,0],S.LINES[k][:,1],color='grey')
        ax.plot(S.LINES[k][:,0],S.LINES[k][:,1],alpha=0.2,color=color)
        l1=S.LINES[k]
        l2=np.column_stack((np.repeat(S.LIST_TPEAK[k]-K['time'],len(l1)),np.linspace(0,0.401e-7,len(l1))))
        vth=intersect(l1,l2)
        VTH+=[vth]
        ax.scatter(S.LIST_TPEAK[k]-K['time'],vth,label='v'+str(j)+'_D'+str(int(l/10)),color=color,marker=marker)
        ax.scatter(S.LIST_TPEAK[k]-K['time'],vth,color='grey',alpha=0.1)

    VTH=np.array([v for v in VTH if str(v) != 'nan'])
    #lowb,uppb= scpstats.t.interval(confidence=0.95, df=len(S.LINES)-1, loc=K['mu'], scale=scpstats.sem(VTH)) 
    xb = np.array([ resample(VTH).mean() for j in range(10000)])
    lowb, uppb = np.quantile(xb, [0.025, 0.975])
    
    ax.hlines(VTH.mean(),-0.8,0,color='grey',linestyle='dashed')
    
    #print(lowb,uppb)
    
    
    #ax.add_patch(Rectangle((-0.8, lowb), 0.8, uppb-lowb,facecolor = 'grey',fill=True,alpha=0.1,label='Bootstrap 95% CI on mean'))
    #ax.add_patch(Rectangle((-0.8, K['mu']-lowb), 0.8, lowb,facecolor = 'grey',fill=True,alpha=0.1))
    #ax.set_title('Cluster'+ str(S.cluster)+'___'+str(K['NPTS'])+'/'+str(K['TPTS'])+' pts'+'___'+'RMSE:'+str(np.round(K['S']*1e8,5))+'e-8'+'___'+r'$\delta$ ='+str(np.round(K['time']*1e3,1)))
    ax.set_ylim([0,0.6e-7])
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Linear Momentum')
    ax.legend(loc='best')

    #print(r'$\bar{v_{th}}$='+str(np.round(K['mu']*1e7,4)))
    #print(r'$10^4 * \Sigma (v - \bar{v_{th}})^2$=' +str(np.round(K['S']*10000,4)))
    #print(r'$\delta$='+str(np.round(K['time'],4)))
    #print('m='+str(np.round(K['m']*1e7,8)))
    #print('pval='+str(np.round(K['pval'],2)))
    return(T,K,ax)

def plot_find_threshold_v0(T): #juste1
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(4)

    ax1.plot(T['time'],T['pval'],color='grey')
    ax2.plot(T['time'],T['S'],color='grey')
    ax3.plot(T['time'],abs(T['m']),color='grey')


    ax0.plot(T['time'],T['NPTS'])
    ax0.hlines(T['TPTS'],min(T.loc[T['pval']>0,"time"]),max(T.loc[T['pval']>0,"time"]),color='tab:blue',linestyle='dashed')
    
    
    ax1.hlines(0.05,min(T.loc[T['pval']>0,"time"]),max(T.loc[T['pval']>0,"time"]),color='grey',linestyle='dashed') #sinon mauvaise normalisation
    ax1.set_ylim([0,1])
    ax1.set_ylabel('pvalue')
    ax2.set_ylabel('MSE')
    ax3.set_ylabel('pente')
    ax3.set_xlabel('screening neural delay')
    plt.subplots_adjust(wspace=0.3,hspace=0.3)
    plt.show()
    
    
def plot_find_threshold_v1(T0,T1,pp=None):
    fig, axs = plt.subplots(nrows=4,ncols=2,figsize=(16,8))

    #axs[1,0].plot(T0['time'],T0['pval'],color='grey')
    #axs[2,0].plot(T0['time'],T0['S'],color='grey')
    #axs[3,0].plot(T0['time'],abs(T0['m']),color='grey')
    
    axs[0,0].plot(T0['time'],T0['NPTS'])
    axs[1,0].plot(T0['time'],T0['pval'],color='grey') #pour continuité
    axs[1,0].plot(np.where(T0['pval'] > 0.05 , T0['time'], None), np.where(T0['pval']> 0.05, T0['pval'], None), color="red")
    axs[1,0].plot(np.where(T0['pval'] <= 0.05, T0['time'], None), np.where(T0['pval'] <= 0.05, T0['pval'], None), color="grey")
    axs[2,0].plot(T0['time'],T0['S'],color='grey') #pour continuité
    axs[2,0].plot(np.where(T0['pval'] > 0.05 , T0['time'], None), np.where(T0['pval']> 0.05, T0['S'], None), color="red")
    axs[2,0].plot(np.where(T0['pval'] <= 0.05, T0['time'], None), np.where(T0['pval'] <= 0.05, T0['S'], None), color="grey")
    axs[3,0].plot(T0['time'],abs(T0['m']),color='grey') #pour continuité
    axs[3,0].plot(np.where(T0['pval'] > 0.05 , T0['time'], None), np.where(T0['pval']> 0.05, abs(T0['m']), None), color="red")
    axs[3,0].plot(np.where(T0['pval'] <= 0.05, T0['time'], None), np.where(T0['pval'] <= 0.05, abs(T0['m']), None), color="grey")

    
    
    try:
        axs[1,0].hlines(0.05,min(T0.loc[T0['pval']>0,"time"]),max(T0.loc[T0['pval']>0,"time"]),color='grey',linestyle='dashed') #sinon mauvaise normalisation
    except:
        pass
    axs[1,0].set_ylim([0,1])
    #axs[0,0].set_ylim([0,T0['TPTS']])
    axs[0,0].set_ylabel('how many points?')
    axs[1,0].set_ylabel('pvalue')
    axs[2,0].set_ylabel('MSE')
    axs[3,0].set_ylabel('pente')
    axs[3,0].set_xlabel('screening neural delay')
    
    axs[0,1].plot(T1['time'],T1['NPTS'])
    axs[1,1].plot(T1['time'],T1['pval'],color='grey') #pour continuité
    axs[1,1].plot(np.where(T1['pval'] > 0.05 , T1['time'], None), np.where(T1['pval']> 0.05, T1['pval'], None), color="red")
    axs[1,1].plot(np.where(T1['pval'] <= 0.05, T1['time'], None), np.where(T1['pval'] <= 0.05, T1['pval'], None), color="grey")
    axs[2,1].plot(T1['time'],T1['S'],color='grey') #pour continuité
    axs[2,1].plot(np.where(T1['pval'] > 0.05 , T1['time'], None), np.where(T1['pval']> 0.05, T1['S'], None), color="red")
    axs[2,1].plot(np.where(T1['pval'] <= 0.05, T1['time'], None), np.where(T1['pval'] <= 0.05, T1['S'], None), color="grey")
    axs[3,1].plot(T1['time'],abs(T1['m']),color='grey') #pour continuité
    axs[3,1].plot(np.where(T1['pval'] > 0.05 , T1['time'], None), np.where(T1['pval']> 0.05, abs(T1['m']), None), color="red")
    axs[3,1].plot(np.where(T1['pval'] <= 0.05, T1['time'], None), np.where(T1['pval'] <= 0.05, abs(T1['m']), None), color="grey")

    
    try:
        axs[0,1].hlines(T1['TPTS'],min(T1.loc[T1['pval']>0,"time"]),max(T1.loc[T1['pval']>0,"time"]),color='tab:blue',linestyle='dashed')
        axs[1,1].hlines(0.05,min(T1.loc[T1['pval']>0,"time"]),max(T1.loc[T1['pval']>0,"time"]),color='grey',linestyle='dashed') #sinon mauvaise normalisation
    except:
        pass
    axs[1,1].set_ylim([0,1])
    #axs[0,1].set_ylim([0,T1['TPTS']])
    axs[0,1].set_ylabel('how many points?')
    axs[1,1].set_ylabel('pvalue')
    axs[2,1].set_ylabel('MSE')
    axs[3,1].set_ylabel('pente')
    axs[3,1].set_xlabel('screening neural delay')
    
    
    plt.subplots_adjust(wspace=0.3,hspace=0.3)
    plt.show()
    
    try:
        pp.savefig(fig)
    except:
        pass
    
def error_threshold(v,D,eps, vth,plotting=False, ax=None, color1='tab:grey', color2='tab:grey',interp=False):
    tCFD,fCFD=return_velo(v,D,interp)
    if interp==False:
        line1=np.column_stack((tCFD[:len(fCFD.dropna())],fCFD.dropna()))

        line2=np.column_stack((tCFD[:len(fCFD.dropna())],np.repeat(vth,len(fCFD.dropna()))))
        firstline=LineString(line1)
        secondline = LineString(line2)   
        intersectionA = firstline.intersection(secondline)
        time_th=intersectionA.xy[0][0]
        mome_th=intersectionA.xy[1][0]

        line2=np.column_stack((tCFD[:len(fCFD.dropna())],np.repeat(vth-eps,len(fCFD.dropna()))))
        firstline=LineString(line1)
        secondline = LineString(line2)   
        intersectionA = firstline.intersection(secondline)
        time_th_low=intersectionA.xy[0][0]
        mome_th_low=intersectionA.xy[1][0]

        line2=np.column_stack((tCFD[:len(fCFD.dropna())],np.repeat(vth+eps,len(fCFD.dropna()))))
        firstline=LineString(line1)
        secondline = LineString(line2)   
        intersectionA = firstline.intersection(secondline)
        time_th_sup=intersectionA.xy[0][0]
        mome_th_sup=intersectionA.xy[1][0]
    else:
        line1=np.column_stack((tCFD,fCFD))

        line2=np.column_stack((tCFD,np.repeat(vth,len(fCFD))))
        firstline=LineString(line1)
        secondline = LineString(line2)   
        intersectionA = firstline.intersection(secondline)
        time_th=intersectionA.xy[0][0]
        mome_th=intersectionA.xy[1][0]

        line2=np.column_stack((tCFD,np.repeat(vth-eps,len(fCFD))))
        firstline=LineString(line1)
        secondline = LineString(line2)   
        intersectionA = firstline.intersection(secondline)
        time_th_low=intersectionA.xy[0][0]
        mome_th_low=intersectionA.xy[1][0]

        line2=np.column_stack((tCFD,np.repeat(vth+eps,len(fCFD))))
        firstline=LineString(line1)
        secondline = LineString(line2)   
        intersectionA = firstline.intersection(secondline)
        time_th_sup=intersectionA.xy[0][0]
        mome_th_sup=intersectionA.xy[1][0]
        
    
    if plotting:
        ax.plot(tCFD,fCFD,color='grey')

        ax.hlines(vth, min(tCFD),time_th,linestyle='dashed',color=color1)
        ax.hlines(vth+eps, min(tCFD),time_th_sup,linestyle='dotted',color=color1)
        ax.hlines(vth-eps, min(tCFD),time_th_low,linestyle='dotted',color=color1)

        ax.vlines(time_th_low, 0,mome_th_low,linestyle='dashed',color=color2)
        ax.vlines(time_th_sup, 0,mome_th_sup,linestyle='dashed',color=color2)
        
        ax.set_xlabel('Time')
        ax.set_ylabel('Momentum')
    
    return(time_th_low,time_th_sup)

