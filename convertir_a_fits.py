import numpy as np
from scipy import signal
import struct
from astropy.io import fits
import datetime
import matplotlib.pyplot as plt
from Tkinter import *
import Tkinter, Tkconstants, tkFileDialog
import os, errno
import sys

#agregue un comentario

def choose_filenameobs():      #funcion para buscar nombre de 'archivo de sesion de obs'
    global appfilename
    global number_pf
    appfilename = tkFileDialog.askopenfilename(initialdir = "/home/regina/Facultad/tesis/pruebas/nuevos_crudos/",\
        title = "Select file",filetypes = (("obs files","*.obs"),("all files","*.*")))  
    label_namefile.configure(text=appfilename,foreground='green')
    nam = open(appfilename)
    renglon = nam.readlines()
    number_pf = int((len(renglon)-17)/15)+1
    fin.set(number_pf)
    nam.close()

def choose_filenamecsv():       #funcion para buscar nombre de 'archivo tabla radiometro'
    global tsis_filename
    tsis_filename = tkFileDialog.askopenfilename(initialdir = "/home/regina/Facultad/tesis/pruebas/nuevos_crudos/",\
        title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))  
    label_namefile2.configure(text=tsis_filename,foreground='green')

def config():           #funcion de configuracion: resolucion, intervalo de observacion, agrupados o no
    global number_pf
    global ini_file
    global fin_file
    global agrup
    global n_canal            #n_canal[0] para banda X(8GHz) y n_canal[1] para banda Ka(32GHz)
    beta = var.get()
    if (beta == 1):
        n_canal = (64,16)
    elif (beta == 2):
        n_canal = (128,32)
    else:
        n_canal = (512,128)
    ini_file = inicio.get()
    fin_file = fin.get()
    agrup = agrupar_state.get()
    if (ini_file<1 or ini_file>number_pf):
        label_check_ini.configure(text='Fuera de rango',foreground='red')
    elif (fin_file<1 or fin_file>number_pf):
        label_check_fin.configure(text='Fuera de rango',foreground='red')
    elif (ini_file>fin_file):
        label_check_ini.configure(text='No es un intervalo',foreground='red')
    else:
        app.destroy()
    
#-------------------------------------------------------------------------------
#------------INICIO VENTANA-----------------------------------------------------
#-------------------------------------------------------------------------------

app = Tk()
app.title(u'Conversi\u00F3n a SDFITS')

number_pf = 1

#-------frame 1-Archivo de sesion de obs----------------------------------------
vp1 = Frame(app)
vp1.grid(column=0, row=0, padx=(100,100), pady=(50,30))
vp1.columnconfigure(0, weight=4)
vp1.rowconfigure(1, weight=4)

label_namefile = Label(vp1, text=u'Seleccione archivo de observaci\u00F3n')
label_namefile.grid(row=0,column=0)
boton_obsfile = Button(vp1,text='Buscar',command=choose_filenameobs)
boton_obsfile.grid(row=0,column=1)

label_namefile2 = Label(vp1, text=u'Seleccione archivo de Tsis')
label_namefile2.grid(row=1,column=0)
boton_csvfile = Button(vp1,text='Buscar',command=choose_filenamecsv)
boton_csvfile.grid(row=1,column=1)

#-------frame 2-Intervalo de obs------------------------------------------------
vp2 = Frame(app)
vp2.grid(column=0, row=1, padx=(100,100), pady=(30,30))

inicio = IntVar()
inicio.set(1)
fin = IntVar()
fin.set(number_pf)
agrupar_state = BooleanVar()
agrupar_state.set(False) #set check state

label_sel = Label(vp2,text='Seleccione observaciones a convertir:')
label_sel.grid(row=0,column=0)

label_inicio = Label(vp2,text='Desde ')
label_inicio.grid(row=1,column=0)
entrada_inicio = Entry(vp2, textvariable=inicio,width=4)
entrada_inicio.grid(row=1,column=1)
label_check_ini = Label(vp2,text='')
label_check_ini.grid(row=1,column=2)

label_fin = Label(vp2,text='Hasta ')
label_fin.grid(row=2,column=0)
entrada_fin = Entry(vp2, textvariable=fin,width=4)
entrada_fin.grid(row=2,column=1)
label_check_fin = Label(vp2,text='')
label_check_fin.grid(row=2,column=2)

agrupar = Checkbutton(vp2, text=u'Agrupar en un \u00FAnico FITS', var=agrupar_state)
agrupar.grid(row=3,column=2)

#-------frame 3-resolucion---------------------------------------------------------
vp3 = Frame(app)
vp3.grid(column=0, row=2, padx=(100,100), pady=(30,30))
var = IntVar()
var.set(3)

label1 = Label(vp3, text=u'Seleccione la resoluci\u00F3n:')
label1.grid(row=0, column=0)

R1 = Radiobutton(vp3, text=u'\u0394v \u2264 1 km/s', variable=var, value=1)#, command=sel_res)
R1.grid(row=1, column=0)

R2 = Radiobutton(vp3, text=u'\u0394v \u2264 0.5 km/s', variable=var, value=2)#, command=sel_res)
R2.grid(row=2, column=0)

R3 = Radiobutton(vp3, text=u'\u0394v \u2264 0.1 km/s', variable=var, value=3)#, command=sel_res)
R3.grid(row=3, column=0)

#-------frame 4-Botones ok y salir-------------------------------------------------------
vp4 = Frame(app)
vp4.grid(column=0, row=3, padx=(100,100), pady=(30,30))

exit = Button(vp4, text='Cancelar', command=quit)
exit.grid(row=0, column=0)

done = Button(vp4, text='Confirmar', command=config)
done.grid(row=0, column=1)

#-------------------------------------------------------------------------------------
app.mainloop()

#---------------------------------------------------------------------------------
#----------------------FIN VENTANA------------------------------------------------
#---------------------------------------------------------------------------------

#-----crear directorio fits-------------------------------------------------------
path = appfilename[:-37] + "fits"
try:  
    os.mkdir(path)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
#-----abrir archivo de sesion de observacion-------------------------------------------

a0 = appfilename
nam = open(a0)
renglon = nam.readlines()
nam.close()

if (renglon[2][6:10] != 'TTCP'):
    print 'comprobar si la observacion se realizo con TTCP, chequear ancho de banda'

#-----abrir archivo auxiliar (plan de observacion)-------------------------------------
#---contiene nombre de objeto(0-8 char), AR(9-24 char), Dec(25-41 char)----------------

#plan_obs = open('/home/reginasar/Documentos/Tesis/programas/plan_observacion/obs_plan157.txt', 'r')
#plan_obs = open('/home/regina/Facultad/tesis/pruebas/plan_observacion/obs_plan_157.txt','r')
plan_obs = open('/home/regina/Facultad/tesis/pruebas/plan_observacion/obs_plan_163.txt','r')
objetos = plan_obs.readlines()
plan_obs.close()

total_exp = len(objetos)
objeto = np.empty((total_exp), dtype=str)
RA_obj = np.empty((total_exp), dtype=str)
DEC_obj = np.empty((total_exp), dtype=str)
for p in range(total_exp):
    objeto[p] = objetos[p][0:8]
#    print objetos[p][9:11],objetos[p][12:14],objetos[p][15:24]
#    print objetos[p][25:28],objetos[p][29:31],objetos[p][32:41]
    RA_obj[p] = float(objetos[p][9:11])*54000 + float(objetos[p][12:14])*900 + float(objetos[p][15:24])*15
    DEC_obj[p] = float(objetos[p][25:28])*3600 + float(objetos[p][29:31])*60 + float(objetos[p][32:41])

#-----abrir archivo radiometro (obs de un minuto)---------------------------------------

#temps = open('/home/reginasar/Documentos/Tesis/Datos/Obs06Jun18Radiometer_clean.csv')
#temps = open('/home/regina/Facultad/tesis/pruebas/nuevos_crudos/Obs06Jun18Radiometer_clean.csv','r')
temps = open(tsis_filename)
tabla_tsis = np.genfromtxt(temps, skip_header=1, dtype=None, usecols=(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19))

#-----fija parametros header, diodo-----------------------------------------------------

header_size = 176                          #bytes por header de PRF
t_diodo_on = 2                             #tiempo con diodo encendido por ciclo
t_diodo_off = 2                            #tiempo con diodo apagado por ciclo
t_ciclo = t_diodo_on + t_diodo_off         #tiempo de cada ciclo

#-----agrupar en unico fits-------------------------------------------------------------
#---si agrupa, genera el header principal del archivo fits------------------------------

if (agrup==True):
    p = ini_file - 1
    a0 = appfilename[:-37] + renglon[17 + p * 15][2:39]
    f0 = open(a0,'rb')
    header = f0.read(header_size)
    label, reclen, recver, stid, spid, ssize, srate, valflag, agflag, \
        rf_to_if_dw, if_to_ch_dw, year, doy, secod, pcosec, chacc_phase, \
        cppc0, cppc1, cppc2, cppc3, emp1, emp2, endlab \
        = struct.unpack_from("<4sI4HI2H2d2HI6d36s40si",header)
    f0.close()
    hdr = fits.Header()                             #genera un header 
    hdr['DATE'] = renglon[14 + p * 15][30:47]       #Fecha de observacion YYYY-DOYTHH:MM:SS
    hdr['ORIGIN'] = 'DSA3'                          #agrega nombre de la antena de origen al header
    hdr['STATID'] = str(stid)                       #agrega identificador de la estacion
    primary_hdu = fits.PrimaryHDU(header=hdr)       #se define este header como el principal
    hdul = fits.HDUList([primary_hdu])              #arma la HDUlist con el header principal

for p in range(ini_file-1,fin_file):

#------GENERAR ARCHIVO FITS----------------------------------------
	#------header principal----------------------    
    if (agrup==False):
        #a0 = '/home/paulabe/Descargas/ttcp-157/2018157/TTCP1/OL/' + renglon[17 + p * 15][2:39]
        a0 = appfilename[:-37] + renglon[17 + p * 15][2:39]
        f0 = open(a0,'rb')
        header = f0.read(header_size)
        label, reclen, recver, stid, spid, ssize, srate, valflag, agflag, \
            rf_to_if_dw, if_to_ch_dw, year, doy, secod, pcosec, chacc_phase, \
            cppc0, cppc1, cppc2, cppc3, emp1, emp2, endlab \
            = struct.unpack_from("<4sI4HI2H2d2HI6d36s40si",header)
        f0.close()
        hdr = fits.Header()                             #genera un header 
        hdr['DATE'] = renglon[14 + p * 15][30:47]       #Fecha de observacion YYYY-DOYTHH:MM:SS
        hdr['ORIGIN'] = 'DSA3'                          #agrega nombre de la antena de origen al header
        hdr['STATID'] = str(stid)                       #agrega identificador de la estacion
        primary_hdu = fits.PrimaryHDU(header=hdr)       #se define este header como el principal
        hdul = fits.HDUList([primary_hdu])              #arma la HDUlist con el header principal

    #------lectura de datos y reduccion----------------------------

    if (rf_to_if_dw < 1e10):                        #fijar N_canal segun la banda
        N_canal = n_canal[0]
    else:
        N_canal = n_canal[1]

    start_time = renglon[14 + p * 15][39:47]
    stop_time = renglon[14 + p * 15][58:66]
    sec_start = int(start_time[0:2]) * 3600 + int(start_time[3:5]) * 60 + int(start_time[6:8])
    sec_stop = int(stop_time[0:2]) * 3600 + int(stop_time[3:5]) * 60 + int(stop_time[6:8])
    t_obs = sec_stop - sec_start - 1
    Pxx = np.zeros((t_obs,2,N_canal*4))
    time_sys = np.zeros(t_obs)
    diode = np.full(t_obs, False, dtype=bool)
    bandera = np.zeros(t_obs)
    tsis_particular = np.zeros(t_obs)
    for l in range(0,14):
        tsis_particular[l*4] = tabla_tsis[p][l]

    #---polarizacion-RR------------    
    for q in range(4):
        #a = '/home/paulabe/Descargas/ttcp-157/2018157/TTCP1/OL/' + renglon[17 + p * 15 + q][2:39]
        a = appfilename[:-37] + renglon[17 + p * 15 + q][2:39]
        f = open(a,'rb')
        print a
        
        for h in range(t_obs):
            header = f.read(header_size)
            label, reclen, recver, stid, spid, ssize, srate, valflag, agflag, \
                rf_to_if_dw, if_to_ch_dw, year, doy, secod, pcosec, chacc_phase, \
                cppc0, cppc1, cppc2, cppc3, emp1, emp2, endlab \
                = struct.unpack_from("<4sI4HI2H2d2HI6d36s40si",header)
            datablock_size = reclen - header_size
            data = f.read(datablock_size)
            datos = struct.unpack_from("<2000000h",data)
            datos = np.reshape(datos, (2,srate), order='F')
            complejo = (datos[0,:] + 0.5) + (datos[1,:] + 0.5) * 1j        
#            frec, Pxx_inter = signal.periodogram(complejo, srate, nfft=N_canal)    #calcula espectro
            frec, Pxx_inter = signal.welch(complejo, srate, nperseg=N_canal, noverlap=0, window=signal.get_window('boxcar',N_canal), nfft=N_canal)
            Pxx[h,0,q*N_canal:(q+1)*N_canal] = np.fft.fftshift(Pxx_inter)    #guarda espectro en un pedacito de array 
            Pxx_inter=np.zeros(N_canal)                                      #limpia el array

        f.close()
    #---polarizacion-LL-------------
    for q in range(4,8):
        #a = '/home/paulabe/Descargas/ttcp-157/2018157/TTCP1/OL/' + renglon[17 + p * 15 + q][2:39]
        a = appfilename[:-37] + renglon[17 + p * 15 + q][2:39]
        f = open(a,'rb')
        print a

        for h in range(t_obs):
            header = f.read(header_size)
            label, reclen, recver, stid, spid, ssize, srate, valflag, agflag, \
                rf_to_if_dw, if_to_ch_dw, year, doy, secod, pcosec, chacc_phase, \
                cppc0, cppc1, cppc2, cppc3, emp1, emp2, endlab \
                = struct.unpack_from("<4sI4HI2H2d2HI6d36s40si",header)
            datablock_size = reclen - header_size
            data = f.read(datablock_size)
            datos = struct.unpack_from("<2000000h",data)
            datos = np.reshape(datos, (2,1000000), order='F')
            complejo = (datos[0,:] + 0.5) + (datos[1,:] + 0.5) * 1j		
#            frec, Pxx_inter = signal.periodogram(complejo, srate, nfft=N_canal)
            frec, Pxx_inter = signal.welch(complejo, srate, nperseg=N_canal, noverlap=0, window=signal.get_window('boxcar',N_canal), nfft=N_canal)

 #           frec = np.fft.fftshift(frec)
 #           plt.plot(frec,Pxx_inter)
 #           plt.show()
            Pxx[h,1,(q-4)*N_canal:(q-3)*N_canal] = np.fft.fftshift(Pxx_inter)
            Pxx_inter = np.zeros(N_canal)

            #---info para header--- 
            bandera[h] += valflag                                #guarda flags 
            time_sys[h] = secod + pcosec * 10**(-12)             #guarda segundos del dia (lo repite 4 veces innecesariamente)
            if (h % t_ciclo < t_diodo_on):                       #guarda estado del calibrador (lo repite 4 veces innecesariamente)
                diode[h] = True
            else:
                diode[h] = False
            bandera[h] += valflag 

        f.close()

#        plt.semilogy(frec, Pxx_medio[:,q])
#        plt.show()
    #---fin reduccion de datos--------------------

    #------columnas SDFITS------------------------
    col0 = fits.Column(name='EXPOSURE', format='1D', array=np.ones(t_obs))         #CORE column. exposicion
    col1 = fits.Column(name='TIMESYS', format='1D', array=time_sys)                #SHARED column. genera columna de tabla fits, segundo del dia
    col2 = fits.Column(name='CAL', format='1L', array=diode)                       #genera columna de tabla fits, estado de calibrador
    col3 = fits.Column(name='TSYS', format='1D', array=tsis_particular)            #CORE column. temperatura de sistema
    col4 = fits.Column(name='DATA', format=str(N_canal*8)+'D', array=Pxx)          #CORE column. data
    col5 = fits.Column(name='FLAGGED', format='1E', array=bandera)                 #genera columna de tabla fits, flags

    t = fits.BinTableHDU.from_columns([col0, col1, col2, col3, col4, col5])              #genera tabla binaria fits dadas las columnas  
    
    #------header extension----------------------
    t.header['TDIM5'] = '(' + str(N_canal*4) + ',2)'    #DATA dimension definition (mandatory)
    t.header['CTYPE1'] = 'FREQ'                         #DATA 1st dimension name
    t.header['CTYPE2'] = 'STOKES'                       #DATA 2nd dimension name
    t.header['EXTNAME'] = 'SINGLE DISH'                 #MANDATORY agrega nombre de la extension al header/cabecera de la extension
    t.header['NMATRIX'] = '1'                           #MANDATORY agrega cantidad de juegos de datos al header de la extension
    t.header['DATE-OBS'] = renglon[14 + p * 15][30:47]  #CORE virtual column
    t.header['TELESCOP'] = 'DSA 3'                      #CORE virtual column
    t.header['SCAN'] = p                                #SHARED virtual column
    t.header['SITELONG'] = -69.398                      #SHARED virtual column
    t.header['SITELAT'] = -35.776                       #SHARED virtual column
    t.header['SITEELEV'] = 1550                         #SHARED virtual column
    t.header['BACKEND'] = 'TTCP'                        #SHARED virtual column
    t.header['BANDWID'] = srate*4                       #CORE virtual column
    t.header['RESTFREQ'] = rf_to_if_dw                   #SHARED virual column
    t.header['FREQRES'] = float(srate) / float(N_canal) #SHARED virual column
    t.header['OBJECT'] = tabla_tsis[p][15]              #CORE virtual column tabla_tsis[p][15] tambien sirve aca
    t.header['OBJ-RA'] = RA_obj[p]                      #virtual column
    t.header['OBJ-DEC'] = DEC_obj[p]                    #virtual column
    t.header['EQUINOX'] = 'J2000.0'                     #virtual column
#    t.header['AZIMUTH'] =                              #SHARED virtual column
#    t.header['ELEVATIO'] = tabla_tsis[p][17]            #SHARED virtual column
#    t.header['TCAL'] =                                 #SHARED virtual column

    hdul.append(t)                                      #agrega la tabla al HDUlist

    #------guarda archivo fits-------------------
    #salid = '/home/reginasar/Documentos/Tesis/Datos/jun06_18' + renglon[17 + p * 15][6:10] + \
	#	renglon[17 + p * 15][23:35] + '_' + str(N_canal*4).zfill(4) + '.fits'     #nombre del archivo con la tabla fits
    if (agrup==False):
        salid = appfilename[:-37] + 'fits/' + renglon[17 + p * 15][6:10] + \
            renglon[17 + p * 15][23:35] + '_' + str(N_canal*4).zfill(4) + '.fits'     #nombre del archivo con la tabla fits
        hdul.writeto(salid, clobber=True)     #guarda la HDUlist en un archivo fits, sobreescribe si ya existe

if (agrup==True):
    p = ini_file - 1
    salid = appfilename[:-37] + 'fits/' + renglon[17 + p * 15][6:10] + \
            renglon[17 + p * 15][23:35] + '_' + str(N_canal*4).zfill(4) + 'multi' + '.fits'     #nombre del archivo con las tablas fits
    hdul.writeto(salid, clobber=True)     #guarda la HDUlist en un archivo fits, sobreescribe si ya existe



