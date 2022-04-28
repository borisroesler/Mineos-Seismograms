#!/usr/bin/env python
import numpy as np
import pandas as pd
import obspy
import os
import glob

model = 'AK135'
EQ_file = 'EQs'
station_file = 'StationsII.csv'
upper_period = 1/500
lower_period = 1/100


#%% Functions

def makeMatrix(MT):
    MT = np.array([[MT[0], MT[3], MT[4]], [MT[3], MT[1], MT[5]], [MT[4], MT[5], MT[2]]])
    return MT


def unmakeMatrix(MT):
    MT = np.array([MT[0,0], MT[1,1], MT[2,2], MT[0,1], MT[0,2], MT[1,2]])
    return MT


def Moment(M):
    
    trace = np.trace(M)
    M_iso = np.diag(np.array([trace,trace,trace]))
    
    M_devi = M - M_iso
    eigenw, _ = np.linalg.eig(M_devi)
    
    M0 = np.sqrt(np.sum((eigenw**2)/2))
    
    return (M0)


def Magnitude(M0):
    return 2/3*(np.log10(M0)-16.1)


def MTtoAngles(MT):
    
    def Strike_and_Dip(v):
        dip = 0.5*np.pi - np.arccos(np.sqrt(v[1]**2 + v[2]**2))
        azimuth = 0.5*np.pi + np.arctan2(v[1],v[2]) 
        
        if v[0] > 0:
            azimuth = azimuth + np.pi
        strike = azimuth + 0.5*np.pi
        
        if strike >= 2*np.pi:
            strike = strike - 2*np.pi
        
        strike = np.degrees(strike)
        dip = np.degrees(dip)
        
        return strike, dip
    
    evals,evecs = np.linalg.eigh(MT)
    i1 = np.argmax(evals)
    i3 = np.argmin(evals)
    
    # eigenvectors
    p = evecs[:,i3]
    t = evecs[:,i1]
    #n = evecs[:,i2]
    
    # fix convention for eigenvectors
    if p[0] < 0: p = -p
    if t[0] < 0: t = -t
    
    pole1 = (t+p)/np.sqrt(2)
    pole2 = (t-p)/np.sqrt(2)
    if p[0] > t[0]: pole1 = -pole1
    
    strike1, dip1 = Strike_and_Dip(pole1)
    strike2, dip2 = Strike_and_Dip(pole2)
    
    rake1 = np.degrees(np.arccos(np.dot(pole2,np.array([0.,-np.cos(np.radians(strike1)),np.sin(np.radians(strike1))]))))
    if pole2[0] < 0: rake1 = -rake1
    
    rake2 = np.degrees(np.arccos(np.dot(pole1,np.array([0.,-np.cos(np.radians(strike2)),np.sin(np.radians(strike2))]))))
    if pole1[0] < 0: rake2 = -rake2
    
    return strike1, dip1, rake1, strike2, dip2, rake2


def makeCMTFiles(filename):
        
    def writeFile(f,**kwargs):
    
        f.write(name + ' ')
        f.write(str(obspy.UTCDateTime(event['date']).year) + ' ')
        f.write(str(obspy.UTCDateTime(event['date']).julday) + ' ')
        date = obspy.UTCDateTime(event['date']+'T'+event['time'])
        f.write(str(date.hour) + ' ' + str(date.minute) + ' ' + '{:.2f}'.format(date.second) + ' ')
        
        if kwargs:
            f.write('{:.2f}'.format(kwargs['latitude']) + ' ' + '{:.2f}'.format(kwargs['longitude']) + ' ')
            f.write(str(kwargs['depth']) + ' ')
        else:
            f.write('{:.2f}'.format(event['latitude']) + ' ' + '{:.2f}'.format(event['longitude']) + ' ')
            f.write(str(event['depth']) + ' ')
        f.write('1.0 ')
        f.write(str(event['half duration']) + ' ')
        
        exp = event['exponent']
            
        MT = np.array([event['Mrr'],event['Mtt'],event['Mpp'],event['Mrt'],event['Mrp'],event['Mtp']])
        f.write('{:.4g}'.format(Moment(makeMatrix(MT))) + ' ')
        MT = MT/10**int(exp)
        f.write(str(round(MT[0],3)) + ' ' + str(round(MT[1],3)) + ' ' + str(round(MT[2],3)) + ' ' + str(round(MT[3],3)) + ' ' + str(round(MT[4],3)) + ' ' + str(round(MT[5],3)) + ' ')
        f.write('1.0e' + str(exp) + ' ')
        
        strike1, dip1, rake1, strike2, dip2, rake2 = MTtoAngles(makeMatrix(MT))
        f.write(str(int(round(strike1))) + ' ')
        f.write(str(int(round(dip1))) + ' ')
        f.write(str(int(round(rake1))) + ' ')
        
        f.write(str(int(round(strike2))) + ' ')
        f.write(str(int(round(dip2))) + ' ')
        f.write(str(int(round(rake2))) + '\n')
    
        f.close()
        
        return
    
    cat = pd.read_csv(filename, delimiter=',')
    
    for i in range(len(cat)):
        event = cat.iloc[i]
        name = event['event name']
        
        # Folder
        try:
            os.mkdir('events/' + name)
        except:
            pass
        
        # Seismograms' CMT files
        file = open('events/' + name + '/CMTSOLUTION', 'w')
        writeFile(file)
    
    return


def Stations(filename):
    
    stations = pd.read_csv(filename, delimiter=',')
    
    f1 = open('aux/stations.site', 'w')
    f2 = open('aux/stations.sitechan', 'w')
    
    for i in range(len(stations)):
        station = stations.iloc[i]
        f1.write('{:<6}'.format(station['Station Code']) + ' ')
        f1.write(' 1980001       -1'+ ' ')
        f1.write('{:9.4f}'.format(station['Latitude']) + ' ')
        f1.write('{:9.4f}'.format(station['Longitude']) + ' ')
        f1.write('{:9.4f}'.format(1e-3*station['Elevation']) + ' ')
        f1.write('                                                   -    -         0.0000    0.0000' + ' ')
        f1.write('01/01/01 00:00:00\n')
        
        f2.write('{:<6}'.format(station['Station Code']) + ' ')
        f2.write('LHE       1980001       -1       -1 n       0.0000   90.0   90.0 -                                                  01/01/01 00:00:00\n')
        f2.write('{:<6}'.format(station['Station Code']) + ' ')
        f2.write('LHN       1980001       -1       -1 n       0.0000    0.0   90.0 -                                                  01/01/01 00:00:00\n')
        f2.write('{:<6}'.format(station['Station Code']) + ' ')
        f2.write('LHZ       1980001       -1       -1 n       0.0000    0.0    0.0 -                                                  01/01/01 00:00:00\n')
    
    f1. close()
    f2. close()
    
    return


def Eigenfunctions(model, resname, ascii=False):
    modes = ['S', 'T']
    fdb_list = open('aux/db_list','w')
    path = os.getcwd()
    os.system('rm -r aux/' + model + '*')
    os.system('rm -r aux/*' + resname + '*')
    for idx, mode in enumerate(modes):
        f = open('parameter_file','w')
        f.write(model + '\n')
        f.write(resname + '_' + mode + '\n')
        f.write('e' + resname + '_' + mode + '\n')
        f.write('1.0e-10 1\n')
        f.write(str(3-idx) + '\n')
        f.write('0 8000 0.0 20.0 0 200\n')
        f.close()
        os.system('bin/minos_bran < parameter_file')

        f = open('parameter_file','w')
        f.write(str(3-idx) + '\n')
        f.write(model + '\n')
        f.write('6371000\n')
        f.write(resname + '_' + mode + '\n')
        f.write('e' + resname + '_' + mode + '\n')
        f.write(model + '_' + mode + '\n')
        f.close()
        
        os.system('bin/eigcon < parameter_file')
        os.system('rm parameter_file')
        os.system('mv *' + resname + '_' + mode + ' aux')
        os.system('mv ' + model + '_' + mode + '* aux')
        
        fdb_list.write(path + '/aux/' + model + '_' + mode + '\n')
    fdb_list.close()
    return


def runGreen(cmtfile):

    parafile = open('parameter_file','w')
    parafile.write('aux/stations\n')
    parafile.write('aux/db_list\n')
    parafile.write(cmtfile + '\n')
    parafile.write('0. 260.\n')

    parafile.write('8000\n')
    parafile.write('green\n')
    parafile.close()
    
    os.system('bin/green < parameter_file')
    if os.path.exists('Syndat.wfdic'):
        os.system('rm -r Syndat.wfdic')
    parafile = open('parameter_file','w')
    parafile.write('CMTSOLUTION\n')
    parafile.write('0\n')
    parafile.write('green\n')
    parafile.write('Syndat\n')
    parafile.write('0\n') # Acceleration
    parafile.close()
    
    os.system('bin/syndat < parameter_file')
    os.system('cp aux/stations.site Syndat.site')
    os.system('cp aux/stations.sitechan Syndat.sitechan')
    os.system('bin/creat_origin ' + ' CMTSOLUTION Syndat')
    os.system('bin/cucss2sac Syndat Syns')
    synall = glob.glob('Syns/*.SAC')
    for syncur in synall:
        #os.rename(syncur, syncur.replace(' ','0'))
        fname = os.path.split(syncur)[-1].rsplit('.', 3)[-3:]
        os.rename(syncur, 'Syns/' + fname[0] + '.' + fname[1] + '.' + fname[2])

    return


#%% Main Program
makeCMTFiles(EQ_file)
Stations(station_file)

os.system('cp Models/' + model + '.txt ' + model + '.txt')
Eigenfunctions(model + '.txt', model)

CMTs = glob.glob('events/*/CMTSOLUTION')
for CMT in CMTs:
    os.system('mv ' + CMT  + ' .') # cp
    runGreen('CMTSOLUTION')
    dire = CMT.replace('/CMTSOLUTION', '')
    try:
        os.mkdir(dire)
    except:
        pass
    os.system('cp Syns/*.SAC ' + dire)
    os.system('rm -r Syn*')
    os.system('rm -r green*')

os.system('rm CMTSOLUTION')
os.system('rm -r aux/*' + model +'*')
os.system('rm aux/db_list')
os.system('rm parameter_file')
os.system('rm ' + model + '.txt')


#%% Processing
cat = pd.read_csv(EQ_file, delimiter=',')
events = cat['event name'].values

data = pd.read_csv(station_file, delimiter=',')
names = data['Station Code'].values
networks = data['Network Code'].values

for idx,event in enumerate(events):
    
    st = obspy.read('events/' + event + '/*.SAC')
    
    st.taper(max_percentage=0.05, type='cosine', max_length=None, side='both')
    st.trim(starttime=st[0].stats.starttime-30*60, endtime=st[0].stats.endtime, pad=True, fill_value=0.)
    st.filter('bandpass', freqmin=upper_period, freqmax=lower_period, corners=4, zerophase=True)
    
    for tr in st:
        tr.stats.network = networks[np.where(names==tr.stats.station)[0][0]] # correct missing network names in file names
        tr.data = 1e-6*tr.data # mm/s**2
        
        tr.write('events/' + event + '/' + tr.id + '.SAC', format = 'SAC')
        os.remove('events/' + event + '/' + tr.stats.station + '.' + tr.stats.channel + '.SAC')

