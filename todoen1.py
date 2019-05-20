import numpy as np
from scipy.optimize import curve_fit
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as plt
import math

def Gauss(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def Gain(el,g,tau):
	return g * np.exp(-tau/np.sin(el))

def temp_model(el,trec,tau):
    return trec + 253 * (1. - np.exp(-tau/np.sin(el)))

#----temperaturas de brillo---------------------------------
frec_central_8Hz= 8450000000                      #frec central Hz
frec_central_32Hz = 31951000000                   #frec central Hz
k = 1.3806488 * 10**(-16)                         #cte de Boltzmann [erg/K]

c = 3 * 10**10                                    #velocidad de la luz cm/s
HPBW_8Hz = 58.4 * c / (frec_central_8Hz* 3500)
HPBW_32Hz = 58.4 * c / (frec_central_32Hz * 3500)
print HPBW_8Hz*60, HPBW_32Hz*60
omega_8Hz = 1.133 * HPBW_8Hz**2 * ((np.pi / (180.))**2) #* 1.133        #Tools of radioastronomy-angulo solido del haz de antena [str], considero HPBW=1arcmin
omega_32Hz = 1.133 * HPBW_32Hz**2 * ((np.pi / (180.))**2)

S_0521365_8GHz_140 = 4.016                        #flujo de la fuente 0521-365 [Jy] DOY=140
S_0521365_8GHz_162 = 3.807                        #flujo de la fuente 0521-365 [Jy] DOY=162
S_0521365_8GHz = S_0521365_8GHz_140 + ((S_0521365_8GHz_162 - S_0521365_8GHz_140)/(162. - 140.)) * (157. - 140.)
S_0521365_32GHz_141 = 5.117
S_0521365_32GHz_273 = 4.253
S_0521365_32GHz = S_0521365_32GHz_141 + ((S_0521365_32GHz_273 - S_0521365_32GHz_141)/(273. - 141.)) * (163. - 141.)
print 'S_8GHz = ',S_0521365_8GHz
print 'S_32GHz = ',S_0521365_32GHz
T_0521365_8GHz = S_0521365_8GHz * 10**(-23) * c**2 /(2 * frec_central_8Hz**2 * k * omega_8Hz)
T_0521365_32GHz = S_0521365_32GHz * 10**(-23) * c**2 /(2 * frec_central_32Hz**2 * k * omega_32Hz)
print 'Tb_8GHz = ',T_0521365_8GHz
print 'Tb_32GHz = ',T_0521365_32GHz


#----lectura de la tabla------------------------------------
b0 = "pot_157RR3.txt"
powfile = open(b0,'r')
scan = powfile.readlines()
total_rows = len(scan)
powfile.close()

diodo_ON = np.zeros((13,total_rows))
diodo_OFF = np.zeros((13,total_rows))

for q in range(total_rows):   
	diodo = np.fromstring(scan[q], dtype=float, sep=' ')
	diodo_ON[:,q] = diodo[:13]
	diodo_OFF[:,q] = diodo[13:26]

#----como varia la potencia en cada obs(t)------------------
for q in range(total_rows):
	plt.plot(range(1,14),diodo_ON[:,q],'.',color='green')
	plt.plot(range(1,14),diodo_OFF[:,q],'.',color='orange')
plt.ylabel('potencia')
plt.xlabel('ciclos')
plt.plot(range(1,14),diodo_ON[:,0],'.',color='green',label='diodo ON')
plt.plot(range(1,14),diodo_OFF[:,0],'.',color='orange',label='diodo OFF')
plt.legend(loc='center right')
plt.savefig('/home/regina/Facultad/tesis/escritos/latex/imagenes/pot_diodoONOFF.png')
plt.show()

Tcal = 23.
print 2 * Tcal * diodo_OFF[0,0] / (diodo_ON[0,0] - diodo_OFF[0,0]) 

diodo_ON_med = np.sqrt(np.mean(diodo_ON,0))
diodo_ON_std = np.std(diodo_ON,0)

diodo_OFF_med = np.sqrt(np.mean(diodo_OFF,0))
diodo_OFF_std = np.std(diodo_OFF,0)

for q in range(total_rows):
	if (np.log(diodo_ON_med[q]/diodo_OFF_med[q])<=0):
		print 'log < 0 --->',q

#---obs antes del mapeo-------------------------------------
#----comparacion diodo on/off siendo on/off source y Tsis---
#-----------------------------------------------------------
comp = np.zeros(10)
Tsis = np.zeros(10)
Tcal = 23.    #[k]de lo que dice el documento

for q in range(10):
	comp[q] = (diodo_ON_med[2*q] - diodo_OFF_med[2*q]) / (diodo_ON_med[2*q+1] - diodo_OFF_med[2*q+1])
	Tsis[q] = np.sqrt(2.) * Tcal * diodo_OFF_med[2*q+1] / (diodo_ON_med[2*q+1] - diodo_OFF_med[2*q+1])
	Ta = Tsis[q] * (diodo_OFF_med[2*q] - diodo_OFF_med[2*q+1]) / (diodo_OFF_med[2*q+1])
print '------------------'
print '---COMPARACION----'
print '------------------'
print 'cociente =', np.mean(comp), '+/-', np.std(comp)
print 'Tsis =', np.mean(Tsis), '+/-', np.std(Tsis)
print 'Ta =', np.mean(Ta), '+/-', np.std(Ta)

#----mapeo-------------------------------------------------
print '------------------'
print '----MAPEO---------'
print '------------------'

Tsis_ref = 2 * Tcal * diodo_OFF_med[19] / (diodo_ON_med[19] - diodo_OFF_med[19]) #obs off de referencia
print 'Tsis de referencia =', Tsis_ref
alfa = np.arange(-1,1.5,0.5)

mapeo = Tsis_ref * (diodo_OFF_med[20:45] - diodo_OFF_med[19]) / (diodo_OFF_med[19]) #temperatura de antena
map_RR_Ta = np.reshape(mapeo, (5,5), order='F')
map_RR_Ta = np.array((map_RR_Ta[::-1,0],map_RR_Ta[:,1],map_RR_Ta[::-1,2],map_RR_Ta[:,3],map_RR_Ta[::-1,4]))
map_RR_Ta = np.transpose(map_RR_Ta)

Tsis = 2 * Tcal * diodo_OFF_med[20:45] / (diodo_ON_med[20:45] - diodo_OFF_med[20:45])
map_Tsis = np.reshape(Tsis, (5,5), order='F')
map_Tsis = np.array((map_Tsis[::-1,0],map_Tsis[:,1],map_Tsis[::-1,2],map_Tsis[:,3],map_Tsis[::-1,4]))
map_Tsis = np.transpose(map_Tsis)

plt.suptitle('Mapeo PKS0521-365 a 8.45GHz \n Ta [K]')
plt.xlabel(u'\u03B1 Cos(\u03B4) [arcmin]')
plt.ylabel(u'\u03B4 [arcmin]')
plt.contourf(alfa,alfa,map_RR_Ta,cmap=plt.cm.bone)
plt.axvline(x=0,linestyle='--',color='w')
plt.plot(alfa,np.zeros(5),linestyle='--',color='w')
plt.plot(alfa,np.zeros(5),'o',color='k')
plt.plot(np.zeros(5),alfa,'o',color='k')
plt.text(-1.0,0.05,'(-1,0)',color='w')
plt.text(-0.5,0.05,'(-0.5,0)',color='w')
plt.text(0.05,0.05,'(0,0)',color='k')
plt.text(0.45,0.05,'(0.5,0)',color='w')
plt.text(0.85,0.05,'(1,0)',color='w')
plt.text(0.05,-1.0,'(0,-1)',color='w')
plt.text(0.05,-0.5,'(0.-0.5)',color='w')
plt.text(0.05,0.9,'(0,1)',color='w')
plt.text(0.05,0.5,'(0,0.5)',color='w')
plt.colorbar()
plt.savefig('/home/regina/Facultad/tesis/escritos/correcciones/imagenes/map_X_dOFF.png')
plt.show()

plt.suptitle('Mapeo PKS0521-365 a 8.45GHz \n Tsis [K]')
plt.xlabel(u'\u03B1 Cos(\u03B4) [arcmin]')
plt.ylabel(u'\u03B4 [arcmin]')
plt.contourf(alfa,alfa,map_Tsis)
plt.colorbar()
plt.savefig('/home/regina/Facultad/tesis/escritos/latex/imagenes/map_X_dTsis.png')
plt.show()

temps = open('/home/regina/Facultad/tesis/pruebas/nuevos_crudos/Obs06Jun18Radiometer_clean.csv')
tabla_tsis = np.genfromtxt(temps, skip_header=1, dtype=None, usecols=(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))

tsis_radiom = np.mean(tabla_tsis[22:47][:],1)
map_tsis_radiom = np.reshape(tsis_radiom, (5,5), order='F')
map_tsis_radiom = np.array((map_tsis_radiom[::-1,0],map_tsis_radiom[:,1],map_tsis_radiom[::-1,2],map_tsis_radiom[:,3],map_tsis_radiom[::-1,4]))
map_tsis_radiom = np.transpose(map_tsis_radiom)

plt.suptitle('Mapeo PKS0521-365 a 8.45GHz \n Tsis [K] (radiomentro)')
plt.xlabel(u'\u03B1 Cos(\u03B4) [arcmin]')
plt.ylabel(u'\u03B4 [arcmin]')
plt.contourf(alfa,alfa,map_tsis_radiom)
plt.colorbar()
plt.savefig('/home/regina/Facultad/tesis/escritos/latex/imagenes/map_X_Tsis.png')
plt.show()


#-------apuntado--ajuste-gaussiano-------------------------
cruz = open('coor_cruz.txt')
coord_a, coord_d = np.loadtxt(cruz, delimiter=' ', usecols=(0,1), unpack='True')
cruz.close()
coord_a = coord_a[0:5] * 15 * math.cos(math.radians(coord_d[0]))* 60.
coord_a = coord_a[::-1]
coord_d = coord_d[5:10]* 60.

mean_a = np.sum(coord_a * map_RR_Ta[2,:]) / np.sum(map_RR_Ta[2,:])      #weighted arithmetic mean
sigma_a = np.sqrt(np.sum(map_RR_Ta[2,:]*(coord_a - mean_a)**2) / np.sum(map_RR_Ta[2,:]))

popt_a, pcov_a = curve_fit(Gauss,coord_a,map_RR_Ta[2,:],p0=(max(map_RR_Ta[2,:]),mean_a,sigma_a))

mean_d = np.sum(coord_d * map_RR_Ta[:,2]) / np.sum(map_RR_Ta[:,2])
sigma_d = np.sqrt(np.sum(map_RR_Ta[:,2]*(coord_d - mean_d)**2) / len(coord_d))

popt_d, pcov_d = curve_fit(Gauss,coord_d,map_RR_Ta[:,2],p0=(max(map_RR_Ta[:,2]),mean_d,sigma_d))
xa = np.arange(min(coord_a),max(coord_a),0.00001)
xd = np.arange(min(coord_d),max(coord_d),0.00001)

print '--------------------'
print '-----APUNTADO-------'
print '--------------------'
print 'alfa',popt_a[1],pcov_a[1,1]
print 'sigma_a',popt_a[2],pcov_a[2,2]
print 'corr alfa [arcmin]= ', (popt_a[1]-coord_a[2]) 
print 'delta',popt_d[1], pcov_d[1,1]
print 'sigma_d',popt_d[2], pcov_d[2,2]
print 'corr delta [arcmin]= ', (popt_d[1]-coord_d[2])

plt.suptitle('Apuntado')
plt.subplot(2,1,1)
plt.plot(coord_a, map_RR_Ta[2,:],'o', coord_a, Gauss(coord_a,*popt_a),'+', xa, Gauss(xa,*popt_a), 'r-')
plt.ylabel(r'$T_A$ [K]')
plt.xlabel(u'\u03B1 Cos(\u03B4) [arcmin]')
plt.ylim(-0.2,2.5)
plt.xlim(coord_a[0]-0.5,coord_a[4]+0.5)
plt.text(coord_a[0]-0.1,0.3,'(-1,0)')
plt.text(coord_a[1]-0.2,1,'(-0.5,0)')
plt.text(coord_a[2]-0.1,2.1,'(0,0)')
plt.text(coord_a[3]+0.1,1.2,'(0.5,0)')
plt.text(coord_a[4],0.4,'(1,0)')
plt.axvline(x=coord_a[2],linestyle='--')
plt.subplot(2,1,2)
plt.plot(coord_d, map_RR_Ta[:,2],'o', coord_d, Gauss(coord_d,*popt_d),'+', xd, Gauss(xd,*popt_d), 'r-')
plt.ylabel(r'$T_A$ [K]')
plt.xlabel(u'\u03B4 [arcmin]')
plt.ylim(-0.2,2.5)
plt.xlim(coord_d[0]-0.5,coord_d[4]+0.5)
plt.text(coord_d[0]-0.1,0.3,'(0,-1)')
plt.text(coord_d[1]-0.2,1,'(0,-0.5)')
plt.text(coord_d[2]-0.1,2.1,'(0,0)')
plt.text(coord_d[3]+0.1,1,'(0,0.5)')
plt.text(coord_d[4],0.3,'(0,1)')
plt.axvline(x=coord_d[2],linestyle='--')
plt.savefig('/home/regina/Facultad/tesis/escritos/correcciones/imagenes/apuntado_d157.png')
plt.show()

#---ganancia-----------------------------------------------

V_sOFF_dOFF = diodo_OFF_med[46::2]   #off source, diodo off
V_sOFF_dON = diodo_ON_med[46::2]     #off source, diodo on
V_sON_dOFF = diodo_OFF_med[45::2]    #on source, diodo off
V_sON_dON = diodo_ON_med[45::2]

print len(V_sOFF_dOFF)

comp_gan = (V_sON_dON - V_sON_dOFF) / (V_sOFF_dON - V_sOFF_dOFF)

Tsis_OFF = 2. * Tcal * V_sOFF_dOFF / (V_sOFF_dON - V_sOFF_dOFF)
Tsis_ON = 2. * Tcal * V_sON_dOFF / (V_sON_dON - V_sON_dOFF)
Ta = Tsis_OFF * (V_sON_dOFF - V_sOFF_dOFF) / V_sOFF_dOFF
#Ta = Tsis_ref * (V_sON_dOFF - diodo_OFF_med[19]) / (diodo_OFF_med[19])

altfile = 'elev_rad.txt'
alt_rad_ON, alt_rad_OFF = np.genfromtxt(altfile, dtype=None, usecols=(0,1),unpack=True)

alt_ON = alt_rad_ON * 180. / np.pi
alt_OFF = alt_rad_OFF * 180. / np.pi

ganancia, cov = curve_fit(Gain,alt_rad_ON,Ta,p0=(1.,-0.01))

temp_mod, cov2 = curve_fit(temp_model,alt_rad_ON,Tsis_ON,p0=(20.,0.01))
#temp_mod, cov2 = curve_fit(temp_model,alt_rad[40::],Tsis_ON[40::],p0=(40.,200.,0.01))

airmass = 1.0/np.sin(alt_rad_ON)
#coef_lin, err_lin = P.polyfit(airmass[40::],Tsis_ON[40::],1,full=True)
coef_lin, err_lin = P.polyfit(airmass,Tsis_ON,1,full=True)
poli_lin = P.Polynomial(coef_lin)

print '--------------------'
print '-----GANANCIA-------'
print '--------------------'
print 'mod atm', ganancia, cov
print 'temp modelo', temp_mod, cov2
print 'pol lin', coef_lin, err_lin

plt.ylabel('Tsis ONsource [K]')
plt.xlabel(u'El[\u00B0]')
plt.plot(alt_ON,Tsis_ON,'.')
plt.plot(alt_ON,temp_model(alt_rad_ON,*temp_mod))
plt.text(43,42.4,r'$T_{sis} = T_{r} + T_{atm}^o (1-e^{-\tau_o\,A})$'+'\n'+r'$T_{r} = 41.952 \pm 0.002$'+'\n'+r'$T_{atm} =253$'+'\n'+r'$\tau=0.0024 \pm 4e-8$')
plt.savefig('/home/regina/Facultad/tesis/escritos/correcciones/imagenes/Tsis_atm_d157.png')
#plt.plot(alt,Tsis_OFF,'.')
plt.show()

plt.xlabel('Masa de aire')
plt.ylabel('Tsis ONsource')
plt.plot(airmass,poli_lin(airmass))
plt.plot(airmass,Tsis_ON,'.')
plt.text(0.013,42.4,'y = ax + b \n'+r'$a=59.38\pm0.08$'+'\n'+r'$b=41\pm1$')
plt.savefig('/home/regina/Facultad/tesis/escritos/latex/imagenes/Tsis_atmlin_d157.png')
plt.show()

plt.ylabel(r'$T_A$ [K]')
plt.xlabel(u'El[\u00B0]')
plt.plot(alt_ON,Ta,'.',color='r',label=r'$T_A^{inst}$')
#plt.plot(alt,Gain(alt_rad_ON,*ganancia))
#plt.plot(alt_ON, Ta*np.exp(4.23243219e-04/np.sin(alt_rad_ON)),'.',label=r'$T_A^c \tau=0.0004$',color='b')
plt.plot(alt_ON, Ta*np.exp(0.0024/np.sin(alt_rad_ON)),'+',label=r'$T_A^c \tau=0.0024$',color='k')
plt.legend(loc='upper right') 
plt.text(55,3.2,r'$T_A^c = T_A^{inst}\, e^{\tau_o\,A}$')
plt.savefig('/home/regina/Facultad/tesis/escritos/correcciones/imagenes/Ta_cor_d157.png')
plt.show()

plt.ylabel('Tsis OFFsource [K]')
plt.xlabel(u'El[\u00B0]')
plt.plot(alt_OFF,Tsis_OFF,'.')
plt.savefig('/home/regina/Facultad/tesis/escritos/correcciones/imagenes/TsisOFF_atm_d157.png')
#plt.plot(alt,Tsis_OFF,'.')
plt.show()

plt.ylabel('comp')
plt.xlabel('tiempo')
plt.plot(alt_ON,comp_gan,'.')
plt.show()


plt.plot(alt_ON,V_sON_dOFF-V_sOFF_dOFF,'.',color='m')
plt.plot(alt_ON,V_sON_dON-V_sOFF_dON,'.',color='c')

plt.show()



V_sOFF_dOFF = diodo_OFF[:,46::2]   #off source, diodo off
V_sOFF_dON = diodo_ON[:,46::2]     #off source, diodo on
Tsis_OFF = 2. * Tcal * V_sOFF_dOFF / (V_sOFF_dON - V_sOFF_dOFF)
for q in range(13):
    print Tsis_OFF[q,:]
    plt.plot(alt_OFF,Tsis_OFF[q,:],'.',color='g')
plt.ylim(38,39)
plt.ylabel('Tsis OFF')
plt.xlabel('elevacion')
plt.savefig('para_manuel.png')
plt.show()
