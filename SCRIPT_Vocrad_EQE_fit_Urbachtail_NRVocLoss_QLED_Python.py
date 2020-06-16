###################################################################################################
## Evaluation of EQE measurement data + Urbach tail to determine the radiative open-circuit voltage, non-radiative voltage loss, and external luminescence quantum efficiency
# by: Lisa Krückemeier & Dane W. deQuilettes  
# (supplemented version of 'SCRIPT_Vocrad_EQE_fit_Urbachtail.m' by Lisa Krückemeier) 

# translated to Python 2.x and 3.x by Christian Wolff
###################################################################################################
##
###################################################################################################
filename = "EQE_Liu_ACSEnergyLett_19_recipeB.dat"

header_lines = 0 #change this if you have header line(s)

measured_VOC = 1.26 #in [V]
measured_ELQY_PLQY = 9.8395e-2 # absolute number, not %
measured_PCE = 20.8
temperature = 300 #  in [°K]

## some optionals:
###########################################################################
## !!!! VARIABLE: ADJUST URBACH ENERGY IF NECESSARY
E_Urbach=0.015;                        # [eV] Urbach energy (usually between 13.5 to 16 meV for metal-halide perovskites)
###########################################################################
# Decide if you want to choose the EQE value for Urbach Tail attachment manually or automatically
manually='no';                          # choose 'yes' or 'no' if you want to chose the transition point manually 'yes' if a manual correction is necessary or automatically 'no'. 'no' attaches Urbach tail midway between data on logscale.
EQE_level=0.01;                         # EQE value at which the Urbach Tail should be attached (only if manually 'yes').
###########################################################################
# Decide if you want to compare your Voc to a survey of literature values
comparsion_literature='yes';            # choose 'yes' or 'no'





###########################################################################

#below there is nothing you need to modify; just SHIFT+ENTER or RUN

from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
import numpy as np
import pandas as pd
import csv


plt.close('all')
coords = []


if header_lines == 0: # in case you have a header
    try :
        df = pd.read_csv(filename, header=None,sep='\t')     
    except IndexError:
        df = pd.read_csv(filename, header=None) 

else:
    try:
        df = pd.read_csv(filename, header=int(header_lines-1),sep='\t') 
        
    except IndexError:
        df = pd.read_csv(filename, header=int(header_lines-1)) 
        


x = df[df.columns[0]];y = df[df.columns[1]];x = np.array(x);y = np.array(y);


if any(x > 10): #check if energy (eV) or wavelength (nm)
    x =1240/x
    
if any(y > 10): #check if EQE is given in (%), if so it's translated to abs. numbers
    y = y/100

if x[1]-x[2]>0: # bring both arrays into correct order (i.e. w.r.t eV increasing) if one started with e.g. wavelength in increasing order e.g. 300nm, 305nm,...
    x.sort()
    y = np.flip(y)
    

xold = x
x = np.linspace(min(x),max(x),500,endpoint = True)
y= np.interp(x,xold,y)
    
coords = [];

def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

#plt.waitforbuttonpress()
if manually == 'no':    
    fig = plt.figure(1);plt.plot(x,y);plt.yscale("log");plt.ylabel('EQE');plt.xlabel('energy [eV]')#;cid = fig.canvas.mpl_connect('button_press_event', onclick)
    
    while True:
        pts = []
        while len(pts) < 2:
            tellme('Select two points for fit')
            pts = np.asarray(plt.ginput(2))
            #if len(pts) < 2:
             #   tellme('Too few points, starting over')
              #  time.sleep(3)  # Wait a second
    
        ph = plt.plot(pts[:, 0], pts[:, 1], '*r', lw=2)
        
        tellme('Happy? press (any) Key for yes,\n mouse click to start over')
        coords = pts
        if plt.waitforbuttonpress():
            plt.close('all')
            break
    



    


####here you push the last part of the script
q=1.602176462e-19;                      #% [As], elementary charge
h_Js=6.62606876e-34;                    #% [Js], Planck's constant
h_eVs=4.135667662e-15;                  #% [eVs], Planck's constant
k = 1.38064852e-23;                     #% [(m^2)kg(s^-2)(K^-1)], Boltzmann constant
T = temperature
VT=(k*T)/q;                             #% [V], 25.8mV thermal voltage at 300K
c=299792458;                            #% [m/s], speed of light c_0
vor=((h_Js*c)/q)/(1e-9);                #% prefactor for converting between energy and wavelength in (eVnm)




x_interp = np.linspace(min(x),max(x),1000,endpoint = True)
y_interp = np.interp(x_interp,x,y)
#y_interp = savgol_filter(y_interp, 51,4)

if manually == 'yes':

    x_interp = x_interp[y_interp>=EQE_level]
    y_interp = y_interp[y_interp>=EQE_level]
    x_extrap = np.linspace(-1,0,500,endpoint = False)+ min(x_interp)
    
    y_extrap = y_interp[0]*np.exp((x_extrap-min(x_interp))/E_Urbach)
else:
    i_x1 = np.abs(x-coords[0][0]).argmin();i_x2 = np.abs(x-coords[1][0]).argmin();
    x1 = x[min(i_x1,i_x2)];x2 = x[max(i_x1,i_x2)];y1 = y[min(i_x1,i_x2)];y2 = y[max(i_x1,i_x2)];

    coords = coords[-2:];
    yfit = y[min(i_x1,i_x2):max(i_x1,i_x2)]
    xfit = x[min(i_x1,i_x2):max(i_x1,i_x2)]-x1
    yfit2 = np.log(yfit)
    log_slope = np.mean(np.diff(yfit2)/np.diff(xfit))
    E_Urbach = 1/log_slope
    x_extrap = np.linspace(-1,0,500,endpoint = False)+ x2
    y_extrap = y2*np.exp((x_extrap-x2)*log_slope)
    
    x_interp = np.linspace(x2,max(x),1000,endpoint = True)
    y_interp = np.interp(x_interp,x[max(i_x1,i_x2):],y[max(i_x1,i_x2):])
    

    
x_new=np.append(x_extrap, x_interp)
y_new=np.append(y_extrap, y_interp)

fig2 = plt.figure();
fig2.canvas.manager.window.move(0,0);plt.plot(x_new,y_new, '-ob', label='inter-/extrapolated $Q_e$');plt.plot(x,y, '-*r', label='original $Q_e$');plt.yscale("log");plt.ylabel('$Q_e$');plt.xlabel('energy [eV]'); plt.ylim((1e-6, 2))
legend = plt.legend(loc='lower right', fontsize='7');

phi_BB = (2*3.14159265*q**3*(x_new)**2)/(h_Js**3*c**2*(np.exp(x_new/VT)-1))
EL = phi_BB*y_new
j0rad = []
#for i in range(len(EL)):
    
j0rad = integrate.cumtrapz(EL,x_new)
j0rad*= q



filename2 = "AM15G.dat"
df3 = pd.read_csv(filename2, header=None) 
 
#print(df3)
energy_AM15 = df3[df3.columns[1]];
energy_AM15 = np.array(energy_AM15);
spectrum_AM15 = df3[df3.columns[2]];
spectrum_AM15 = np.array(spectrum_AM15);

spectrum_AM15G_interp = np.interp(x_new, energy_AM15,spectrum_AM15)
JSC_calc= integrate.cumtrapz(y_new*spectrum_AM15G_interp,x_new)
JSC_int = max(JSC_calc*q*1e4)

VOC_rad_calc= VT *np.log(JSC_int/max(j0rad))
   
dEQE_interp = np.diff(y_new)/np.diff(np.flip(-x_new))
E_G = x_new[dEQE_interp.argmax()]
#print('band gap is max(d/dE (EQE) = %s eV' %(E_G))


df2 = pd.read_excel("LiteratureSurvey.xlsx") 
VOC_ref = df2[df2.columns[2]];
VOC_ref = np.array(VOC_ref);
VOC_nr = df2[df2.columns[3]];
VOC_nr = np.array(VOC_nr);
ERE_ref = df2[df2.columns[5]];
ERE_ref = np.array(ERE_ref);
PCE_ref = df2[df2.columns[8]];
PCE_ref = np.array(PCE_ref);
ERE_calc= np.exp((measured_VOC-VOC_rad_calc)/VT)*100
fig3 = plt.figure()
fig3.canvas.manager.window.move(45,45);ax = fig3.add_subplot(131);ax2 =fig3.add_subplot(132); ax.plot(x_new,y_new, label='$Q_e$'),ax.plot(x_new[1:],dEQE_interp/max(dEQE_interp),'-', label='$d/dE (EQE)$');ax.set_ylabel('$Q_e$');ax.set_xlabel('energy [eV]');ax.set_ylim((0,1));
ax2.get_yaxis().set_visible(False)
ax2.plot(x_new,y_new,label='$Q_e$');ax2.set_ylim((0,1)); ax2.plot(x_new,EL/max(EL), label='norm. $Q_{e, lum.}$') ;ax2.set_ylabel('$Q_e$');ax2.set_xlabel('energy [eV]');
ax3 =fig3.add_subplot(133);ax3.plot(x_new,y_new,label='$Q_e$'), ax3.plot(x_new,EL/max(EL), label='norm. $Q_{e, lum.}$');ax3.set_xlabel('energy [eV]');ax3.set_yscale('log');ax3.set_ylim((1e-4,2));ax3.yaxis.tick_right()

legend = ax.legend(loc='lower center', fontsize='7');legend2 = ax2.legend(loc='lower center', fontsize='7');legend3 = ax3.legend(loc='lower right',fontsize='7')


if comparsion_literature =='yes':
    fig5 = plt.figure()
    fig5.canvas.manager.window.move(90,90)
    ax7 = fig5.add_subplot(111)
    
    
    ax7.semilogx(ERE_ref,PCE_ref, '*', label='literature'), ax7.semilogx(measured_ELQY_PLQY*100,measured_PCE, 'or', label='this work given $Q_e$');
    ax7.semilogx(ERE_calc,measured_PCE, 'ok', label='this work calcuated $Q_e$');
    ax7.set_ylabel('PCE[%]')
    ax7.set_xlabel('ERE[%]')
    ax7.set_xlim(1e-9,200)
    legend = ax7.legend(loc='lower right', fontsize='7')
    # Set scond x-axis
    ax8 = ax7.twiny()
    ax8.set_xscale("log")
    # Decide the ticklabel position in the new x-axis,
    # then convert them to the position in the old x-axis
    newlabel = list([0.6, 0.5, 0.4, 0.3, 0.2, 0.1]) # labels of the xticklabels: the position in the new x-axis
    new_pos_value = np.linspace(0.6,0.1,6);
    newpos  = 100*np.exp(-1*new_pos_value/VT)   # position of the xticklabels in the old x-axis
    #print(newpos, newlabel)
    ax8.set_xticks(newpos)
    ax8.set_xticklabels(newlabel)
    
    ax8.set_xlabel('$\Delta V_{OC}^{nr}$ [V]')
    ax8.set_xlim(ax7.get_xlim());

width1 = 0.35          
fig6 = plt.figure()
ax9 = fig6.add_subplot(111)
fig6.canvas.manager.window.move(135,135)
ax9.bar([1], E_G, width1, label='$E_G$=%.3f eV'%(E_G))
ax9.bar([1], VOC_rad_calc, width1, label='$V_{OC, rad.}$=%.3f V'%(VOC_rad_calc))
ax9.bar([1], measured_VOC, width1, label='$V_{OC, meas.}$=%.3f V'%(measured_VOC))
ax9.arrow(1-width1/2+0.07, (VOC_rad_calc)*0.999, 0, (measured_VOC-VOC_rad_calc)*0.5, width = 0.0025,facecolor ='k') 
ax9.text(1-width1/2+0.085, measured_VOC*1.01 , '$\Delta V_{OC, nr.}$=%.3f V'%(VOC_rad_calc-measured_VOC), fontsize=7)
ax9.set_xlim((0.5,1.5))
legend = ax9.legend(loc='upper right', fontsize='7')
ax9.set_xticks([])
ax9.set_ylabel('photovoltage [V]')

plt.show()

array1 = JSC_int, max(j0rad), VOC_rad_calc, measured_VOC,measured_ELQY_PLQY*100,ERE_calc,E_G,E_Urbach

with open('overview_results.txt', 'w') as f:
    csv.writer(f, delimiter=',').writerow(['J_SC [A/m2]', 'J0_rad [A/m2]', 'VOC_rad [V]', 'measured VOC [V]', 'ERE_measured [%]', 'ERE_calc[%]', 'E_G [eV]', 'E_Urbach [eV]'] )
    csv.writer(f, delimiter=',').writerow(array1 )


array2 = [np.transpose(x_new[:-1]),np.transpose(y_new[:-1]),np.transpose(j0rad),np.transpose(JSC_calc)]

with open('detailed_results.txt', 'w') as f:
    csv.writer(f, delimiter=',').writerow(['energy [eV]', 'EQE_PV [unitless]', 'J0_rad [A/m2]', 'J_SC [A/m2]'] )
    for i in range(len(j0rad)):
        csv.writer(f, delimiter=',').writerow([x_new[i],y_new[i],j0rad[i],JSC_calc[i]*q*1e4])
    





