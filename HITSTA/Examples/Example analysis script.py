# -*- coding: utf-8 -*-
#%%
"""
Created on Sun Oct 13 21:24:49 2024

@author: can
"""

import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy import stats
import seaborn as sns
from scipy.optimize import curve_fit
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
import ruptures as rpt
from sklearn.linear_model import RANSACRegressor, LinearRegression
from scipy.special import gammainc, gamma

def bandedge_func(x,a_const,b_const,l_const):

    return a_const*(0.5-0.5*np.tanh( (x-b_const)/l_const ))

def reflectance_func(x,a_const,b_const,l_const,r_const,a2_const,b2_const,l2_const,r2_const,c_const):
    return bandedge_func(x,a_const,b_const,l_const,r_const)+bandedge_func(x,a2_const,b2_const,l2_const,r2_const)+c_const

def import_HITSTA(filenames):
    """
    Expects "filenames" to be a single str or else a list of str (paths to data files)
    
    Returns a dictionary with keys "ID1", "ID2", "ID3", etc. and values which are also dictionaries with keys "Reflectance", "Times", ... etc.
    
    """
    #If a single filename is received, make it a one-element list for compatibility with multi-file imports.
    if type(filenames) == str:
        filenames = [filenames]
    
    global exp_func, ex_func, expstetch_func, linear_func, plfit_func, gaussian_func #Make fitting functions available outside this function
    
    
    def exp_func(x,a_const,t_const,c_const):
        return a_const*np.exp(-x/t_const)+c_const
    
    def ex_func(x,a_const,t_const,c_const):
        return -a_const*np.exp(-x/t_const)+a_const+c_const
    
    def expstretch_func(x,a_const,t_const,c_const,r_const):
        return a_const*np.exp(-(x/t_const)**r_const)+c_const
    
    def stretched_exp_definite_integral(t0, t1, tau, beta, A=1.0):
        x0 = (np.array(t0) / tau)**beta
        x1 = (np.array(t1) / tau)**beta
        pre = tau / beta
        gamma_term = gamma(1 / beta)
        return A * pre * gamma_term * (gammainc(1 / beta, x1) - gammainc(1 / beta, x0))
    
    def linear_func(x,t_const,b_const):
        return x/t_const+b_const
    
    def plfit_func(x,ashort_const,tshort_const,along_const,tlong_const):
        return ashort_const*np.exp(-x/tshort_const)+along_const*(x**0.5)*np.exp(-x/tlong_const)
    
    def smooth(y,N):
        ysmooth = []
        for row in y:
            ysmooth.append(np.convolve(row, np.ones(N)/N, mode='same'))
        return np.array(ysmooth)
    
    def gaussian_func(wavelengths, A, mean, std, constant):
        return A*np.exp(-(wavelengths-mean)**2/std**2)+constant
    
    def parse_section(section):
        section = section.strip("\t\n")
        section = section.replace(",",".")
        ID = section.split("\n")[0]
        IDnum = int(ID[2:])
        csvStringIO = StringIO(section)
        dataframe = pd.read_csv(csvStringIO,delimiter="\t",header=1)
        dataframe = dataframe.loc[:, ~dataframe.columns.str.contains('^Unnamed')]
        return ID, dataframe
    
    all_sections = []
    for filename in filenames:
        with open(filename, errors = "ignore") as f:
            contents = f.read()
        sections = contents.split("#")
        all_sections.append(sections)
    
    sections_merged = [all_sections[0][0]] #Take header from first file
    for i in range(1,len(all_sections[0])): #Skip header
        for j in range(0,len(filenames)):
            print("Merging section",i,j)
            if j==0:
                IDstring, df = parse_section(all_sections[0][i])
            else:
                _, df_add = parse_section(all_sections[j][i])
                df_add["Time (h)"] = df_add["Time (h)"]+df["Time (h)"].iloc[-1]
                df = pd.concat([df,df_add])
        sections_merged.append(IDstring + "\n" + df.to_csv(sep="\t"))
    
    header = sections[0]
    exp = {}
    #Read in the raw data
    for section in sections_merged[1:50]:
        IDstring, df = parse_section(section)
        wavelengths = df.columns[4:].to_numpy(dtype=float)
        wavelengths = np.linspace(min(wavelengths),max(wavelengths),len(wavelengths))
        exp[IDstring] = { 
            "Rounds" : df[df["Type (D/W/L)"]=="D"].to_numpy()[:,0]-1,
            "Times" : df[df["Type (D/W/L)"]=="D"]["Time (h)"].to_numpy(),       
            "Wavelengths" : wavelengths,
            "Dark Raw" : df[df["Type (D/W/L)"]=="D"].iloc[:,4:].to_numpy(dtype=float)/np.tile( df[df["Type (D/W/L)"]=="D"].iloc[:,3].to_numpy(dtype=float),(1600,1)).transpose(),
            "Reflectance Raw" : df[df["Type (D/W/L)"]=="W"].iloc[:,4:].to_numpy(dtype=float)/np.tile(df[df["Type (D/W/L)"]=="W"].iloc[:,3].to_numpy(dtype=float),(1600,1)).transpose(),
            "Laser Raw" : df[df["Type (D/W/L)"]=="L"].iloc[:,4:].to_numpy(dtype=float)/np.tile(df[df["Type (D/W/L)"]=="L"].iloc[:,3].to_numpy(dtype=float),(1600,1)).transpose(),
            }
        print(IDstring,len(exp[IDstring]["Reflectance Raw"]),len(exp[IDstring]["Dark Raw"]))
        
    #Time-dependent quantities
    for IDstring in exp.keys():
        ## --General-- 
        rounds = len(exp[IDstring]["Rounds"])
        inds = (wavelengths>650)&(wavelengths<850)
        wavelengths_cut = wavelengths[inds]
        
        ## --PL--
        exp[IDstring]["PL"] = exp[IDstring]["Laser Raw"]-exp[IDstring]["Dark Raw"]
        inds_PLsubtr = wavelengths>920
        exp[IDstring]["PL Subtr"] = exp[IDstring]["PL"]-np.tile(np.expand_dims(np.mean(exp[IDstring]["PL"][:,inds_PLsubtr],axis=1),1),(1,1600))
        exp[IDstring]["PL Peak Intensity"] = np.max(exp[IDstring]["PL Subtr"], axis=1 )
        exp[IDstring]["PL Peak Wavelength"] = exp[IDstring]["Wavelengths"][np.argmax(exp[IDstring]["PL Subtr"], axis=1 )]
        exp[IDstring]["PL Fitted"] = []
        exp[IDstring]["PL Fit Parameters"] = []
        exp[IDstring]["PL Fitted"] = []
        exp[IDstring]["PL Fit Parameters"] = []
    
        for idt, PL_measurement in enumerate(exp[IDstring]["PL Subtr"]):
            try:
                if idt==0:
                    #gaussian_func(wavelengths, A, mean, std, constant):
                    popt, pcov= curve_fit(gaussian_func,wavelengths_cut, PL_measurement[inds] , maxfev=20000, p0=[2,750,100,0],bounds=([0,500,10,-0.005],[200,900,150,0.005]))
                else:
                    popt, pcov= curve_fit(gaussian_func,wavelengths_cut, PL_measurement[inds] , maxfev=20000, p0=popt,bounds=([0,500,10,-0.005],[200,900,150,0.005]))
                exp[IDstring]["PL Fitted"].append(gaussian_func(wavelengths,*popt))
                exp[IDstring]["PL Fit Parameters"].append(popt)
            except:
                print([IDstring," failed PL fit"])
        exp[IDstring]["PL Fitted"] = np.array(exp[IDstring]["PL Fitted"])
        exp[IDstring]["PL Fit Parameters"] = np.array(exp[IDstring]["PL Fit Parameters"])
        if  exp[IDstring]["PL Fit Parameters"][0,0]>1:
            exp[IDstring]["Bandgap (initial PL)"] = 1240/exp[IDstring]["PL Fit Parameters"][0,1]
        else:
            exp[IDstring]["Bandgap (initial PL)"] = np.nan

        #PL self-similarity
        PL = exp[IDstring]["PL Subtr"]
        PL0 = np.tile(exp[IDstring]["PL Subtr"][0,:],(len(exp[IDstring]["Times"]),1))
        exp[IDstring]["PL Self-Similarity"] = np.sum(PL*PL0,axis=1)/( np.sqrt(np.sum(PL0*PL0,axis=1))*np.sqrt(np.sum(PL*PL,axis=1)))
        
        ## --Reflectance--
        exp[IDstring]["Reflectance Raw Subtr"] = exp[IDstring]["Reflectance Raw"]-exp[IDstring]["Dark Raw"]
        Refl_ID2_0 = np.tile(exp["ID2"]["Reflectance Raw"][0],(len(exp[IDstring]["Rounds"]),1))
        Dark_ID2_0 = np.tile(exp["ID2"]["Dark Raw"][0],(len(exp[IDstring]["Rounds"]),1))
        Refl_denominator = np.abs(exp["ID2"]["Reflectance Raw"]-exp["ID2"]["Dark Raw"] + 0.01) 
        exp[IDstring]["Reflectance"] = exp[IDstring]["Reflectance Raw Subtr"]/Refl_denominator
        Reflectance = exp[IDstring]["Reflectance"][0]
        Rmid = 0.5*(np.max(Reflectance[inds])+np.min(Reflectance[inds]))
        index_bandedge = np.argmin( np.abs(Reflectance[inds]-Rmid))
        wavelength_bandedge = wavelengths_cut[index_bandedge]
        R0_bandedge = Reflectance[inds][index_bandedge]
        exp[IDstring]["Bandedge"] = (wavelength_bandedge,R0_bandedge)
        exp[IDstring]["Rmid"] = Rmid
        bandedge_interval = 100
        R_slopes = []
        for Reflectance in exp[IDstring]["Reflectance"]:
            x=wavelengths_cut[max( [(index_bandedge-bandedge_interval),0] ):min( [(index_bandedge+bandedge_interval),len(wavelengths_cut)] )]
            y=Reflectance[inds][max( [(index_bandedge-bandedge_interval),0] ):min( [(index_bandedge+bandedge_interval),len(wavelengths_cut)] )]
            y = np.array(y,dtype=float)
            regression = stats.linregress(x,y)
            if np.isnan(regression.slope):
                print("nan value")
            R_slopes.append(regression.slope) 
        exp[IDstring]["R_slopes (raw)"] = np.array(R_slopes)
        exp[IDstring]["R_slopes (norm.)"] = np.array(R_slopes)/R_slopes[0]
        exp[IDstring]["Eg (estimate)"] = 1240/(wavelength_bandedge-Rmid/R_slopes[0])
        exp[IDstring]["Eg (estimate)"] = 1240/(wavelength_bandedge-50)
        
        #Step near 530
        ind_min = np.argmax(exp[IDstring]["Wavelengths"]>520)
        ind_max = np.argmax(exp[IDstring]["Wavelengths"]>560)
        exp[IDstring]["Short-Wavelength Step"] = exp[IDstring]["Reflectance"][:,ind_max]-exp[IDstring]["Reflectance"][:,ind_min]
        
        
    #Metrics
    early = 3
    for IDstring in exp.keys():  
        #Short-wavelength step fit
        X = exp[IDstring]["Times"]
        Y = exp[IDstring]["Short-Wavelength Step"]
        # X = X[~np.isnan(Y)&~np.isinf(Y)]
        # Y = Y[~np.isnan(Y)&~np.isinf(Y)]
        print(len(X)==len(Y))
        try:
            popt, pcov = curve_fit(ex_func,X,Y,p0=[0.1,5,0.05],bounds=([0,0.05,0],[0.5,1000,0.5]))
            exp[IDstring]["SWS Fit Covariance"] = pcov
            exp[IDstring]["SWS Fit Parameters"] = popt
            exp[IDstring]["SWS Fit"] = (X,ex_func(X,*popt))
        except Exception as e:            
            exp[IDstring]["SWS Fit Covariance"] = np.nan
            exp[IDstring]["SWS Fit Parameters"] = [np.nan]*2
            exp[IDstring]["SWS Fit"] = [np.nan]*2
            print("SWS fit failed on",IDstring,"with exception:",e)
        SWS_score = 100*(1-Y/0.6)
        if SWS_score[0] < 80:
            exp[IDstring]["T80_SWS_score"] = 0
        elif all(SWS_score>80):
            popt, pcov = curve_fit(linear_func,X,SWS_score,p0=[-100,100],bounds=([-1000,90],[-50,100]))
            exp[IDstring]["T80_SWS_score"] = (80-popt[1])*popt[0]            
        else:
            exp[IDstring]["T80_SWS_score"] = X[np.argmax(SWS_score<80)]
        
        #PL lifetime fits
        skip = 0
        ind_min =  max([np.argmax(exp[IDstring]["PL Fit Parameters"][skip:,0])+skip,skip])
        # ind_max = np.argmax( np.diff( exp[IDstring]["PL Fit Parameters"][:,1] )<-20 )
        ind_max = len(exp[IDstring]["Times"])
        # if ind_max == 0:
        #     ind_max = len(exp[IDstring]["Times"])
        X = exp[IDstring]["Times"][ind_min:ind_max]-exp[IDstring]["Times"][ind_min]
        Y = exp[IDstring]["PL Fit Parameters"][ind_min:ind_max,0]
        X_test = exp[IDstring]["Times"].reshape(-1,1)
        Y_test = np.log(exp[IDstring]["PL Fit Parameters"][:,0]).reshape(-1,1)
        try:            
            model = "l2" 
            algo  = rpt.Pelt(model=model).fit(Y_test) # PELT is O(N)
            sigma2 = np.var(y[-len(y)//4:])
            pen    = 3 * sigma2 * np.log(len(y))
            bkps = algo.predict(pen=pen)
            ind_min = bkps[0]
                    
            X = exp[IDstring]["Times"][ind_min:ind_max]-exp[IDstring]["Times"][ind_min]
            Y = exp[IDstring]["PL Fit Parameters"][ind_min:ind_max,0]
            popt, pcov = curve_fit(expstretch_func,X,Y,p0=[10,10,0.5,0.75],bounds=([0,0,0.49,0.5],[100,1000,0.51,1]))
            exp[IDstring]["PL Lifetime Fit Parameters"] = popt
            exp[IDstring]["PL Lifetime Fit"] = (X+exp[IDstring]["Times"][ind_min],expstretch_func(X,*popt))
            #def stretched_exp_definite_integral(t0, t1, tau, beta, A=1.0):
            tmin = exp[IDstring]["Times"][ind_min]
            tmax = 1000
            exp[IDstring]["PL Avg. 1000h"] = 1/1000*(
                                                np.trapz(exp[IDstring]["PL Fit Parameters"][:ind_min,0],x=exp[IDstring]["Times"][:ind_min])
                                              + stretched_exp_definite_integral(tmin,tmax, popt[1], popt[3], A=popt[0])
                                              )

        except Exception as e:            
            exp[IDstring]["PL Lifetime Fit Parameters"] = [np.nan]*4
            exp[IDstring]["PL Lifetime Fit"] = [np.nan]*2
            print("PL lifetime fit failed on",IDstring,"with exception:",e)  
            
            
            
        # bandedge_func(x,a_const,b_const,l_const,r_const)    
        ind_min = 0
        X = exp[IDstring]["Times"][ind_min:ind_max]-exp[IDstring]["Times"][ind_min]
        Y = exp[IDstring]["R_slopes (norm.)"][ind_min:ind_max]
        try:
            popt, pcov = curve_fit(bandedge_func,X,Y,p0=[1,10,20],bounds=([0.9,0,1],[1.2,150,100]))
            # popt, pcov = curve_fit(expstretch_func,X,Y,p0=[8,0.3,0.5,5],bounds=([0,001,0,2],[20,5,20,100]))
            exp[IDstring]["R_slopes Fit Parameters"] = popt
            exp[IDstring]["R_slopes Fit"] = (X+exp[IDstring]["Times"][ind_min],bandedge_func(X,*popt))
            exp[IDstring]["R_slopes Fit T80"] = popt[1]+popt[2]*np.arctanh(-0.3/0.5)
        except Exception as e:            
            exp[IDstring]["R_slopes Fit Parameters"] = [np.nan]*4
            exp[IDstring]["R_slopes Fit"] = [np.nan]*2
            exp[IDstring]["R_slopes Fit T80"] = np.nan
            print("R_slopes Fit failed on",IDstring,"with exception:",e)  

             
        
        #Deg rate based on BES
        try:
            X = exp[IDstring]["Times"]
            Y = exp[IDstring]["R_slopes (norm.)"]
            inds = (X>1)
            popt, pcov = curve_fit(linear_func,X[inds], Y[inds], maxfev=20000, p0=[1,0.9],bounds=([0.01,0],[1e3,1]))

            exp[IDstring]["BES Fit"] = (X,linear_func(X,*popt))
            exp[IDstring]["BES Fit Parameters"] = popt
        except Exception as e:
            print("Failed deg. rate fitting on ",IDstring,"with exception: ",e)
            exp[IDstring]["BES Fit"] = np.nan
            exp[IDstring]["BES Fit Parameters"] = np.nan
        exp[IDstring]["BES Final"] = exp[IDstring]["R_slopes (norm.)"][-1]
        exp[IDstring]["SWS Final"] = exp[IDstring]["Short-Wavelength Step"][-1]
        

    return exp

exp = import_HITSTA([r"HITSTA/EXAMPLES/Example data .txt"])

#%% Transflectance spectrum (single sample) over time

colormap = np.linspace(0,1,30)
colors = plt.cm.coolwarm(colormap)  
fig, ax = plt.subplots(dpi=600)

ID = "ID4"
for rnd in exp[ID]["Rounds"]:
   plt.plot(exp[ID]["Wavelengths"],exp[ID]["Reflectance"][rnd], color = colors[rnd],label=rnd, linewidth=1.5)
   plt.plot(exp[ID]["Bandedge"][0],exp[ID]["Bandedge"][1],'o')


plt.xlim([600,900])
plt.ylim([0,1])
plt.legend(fontsize=4, loc="center right")
plt.ylabel ("Transflectance", fontsize=10)
plt.xlabel ("Wavelength (nm)", fontsize = 10)


#%% PL Spectrum (single sample) over time

fig, ax = plt.subplots(dpi=600)

colormap = np.linspace(0,1,10)
colors = plt.cm.coolwarm(colormap) 
 
labels= ["As dep","30 min aging-HITSTA"]

ID = "ID21"
testX = []
testY = []
for rnd in exp[ID]["Rounds"]:
    plt.plot(exp[ID]["Wavelengths"],exp[ID]["PL Subtr"][rnd], color =colors[rnd],label=str(rnd), linewidth=1.5)

plt.xlim([600,900])
plt.legend(loc="center right", fontsize =4)
plt.ylabel ("PL intensity (counts)", fontsize = 10)
plt.xlabel ("Wavelength (nm)", fontsize = 10)
    
#%% PL lifetime fit

fig, ax = plt.subplots(dpi=600)
colormap = np.linspace(0,1,6)
colors = plt.cm.tab20b(colormap)  

excluded_ids = ["ID1", "ID2"]
included_ids =["ID13","ID16","ID18","ID20","ID21"]
labels = ["MACl","PbCl2","FACl","CsCl","Control"]


for idx, ID in enumerate(included_ids):
    if ID in excluded_ids:
        continue
    if exp[ID]["PL Fit Parameters"][0,0]>0.0:
        plt.plot(exp[ID]["Times"], exp[ID]["PL Fit Parameters"][:,0], 'o-', linewidth=2,
                 label=labels[idx], markersize=4, color=colors[idx])

plt.ylabel ("PL Peak Intensity", fontsize = 10)
plt.xlabel ("Time (h)", fontsize = 10)
plt.legend(loc="upper right", fontsize=8)


#%% Anneal time - PL intensity (initial)

df = pd.read_excel(r'PATH/TO/RUNSHEET.xlsx')
fig, ax1 = plt.subplots(dpi=600)
X=[]
Y=[]
Z=[]
for ID in exp.keys():
    if ID not in ["ID1", "ID2","ID21","ID22","ID23","ID24"]:
        Y.append(exp[ID]["PL Peak Intensity"][0])
        row = df[df["HITSTA"]==int(ID[2:])]
        X.append(row["Anneal time [m]"].iloc[0])
        Z.append(ID)       

plt.scatter(X,Y)
plt.xlabel("Anneal time [m]",fontsize = 10)
plt.ylabel("PL Peak Intensity",fontsize = 10)

for i, point in enumerate(zip(X,Y)):
    plt.annotate(Z[i],point,fontsize=7)
    

#%% Catplot

# Extract relevant data from exp
pl_data = []
for ID in exp:
    if ID not in ["ID1", "ID2"]:
        hitsta_id = int(ID[2:])
        pl_peak = exp[ID]["PL Peak Intensity"][0]
        pl_data.append((hitsta_id, pl_peak))

# Create DataFrame from PL data
pl_df = pd.DataFrame(pl_data, columns=["HITSTA", "PL Peak Intensity"])

# Merge with experiment dataframe
merged_df = df.merge(pl_df, on="HITSTA")

# Focus only on parameters of interest and melt
melted = merged_df.melt(
    id_vars=["HITSTA", "PL Peak Intensity"],
    value_vars=["FAI [M]", "DMAI [M]", "PbI2 [M]", "CsI [M]"],
    var_name="Parameter",
    value_name="Value"
)

# Plot using catplot
sns.set(style="whitegrid")
g = sns.catplot(
    data=melted,
    x="Value",
    y="PL Peak Intensity",
    col="Parameter",
    kind="strip",  # or "box", "violin", "swarm"
    col_wrap=2,
    height=4,
    aspect=1.2,
    jitter=True
)
g.set_titles("{col_name}")
plt.tight_layout()
plt.show()
