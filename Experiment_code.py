import time
import ADS1256
import RPi.GPIO as GPIO
import pandas as pd
import csv 
import matplotlib.pyplot as plt
import numpy as np
import asyncio
import queue
import board
import busio
import adafruit_mcp4725
from scipy.signal import chirp, find_peaks, detrend, welch, resample, decimate
import tkinter
from tkinter import Tk, filedialog, simpledialog
import os.path
from multiprocessing import Process, Lock

#===== modal identification

import CESSIPy_modRenan as SSI 
from MRPy import MRPy #Library with modal analysis functions
import auxFunctions_postProcessEMMARM as auxEMMARM #Library with auxiliary functions
from tkinter import Tk, filedialog

###         For modal identification
#=======================================================================================
#Global variables
typeOfProcessing = ['singleFile','batchProcessing']
typeOfSystem = ['uEMMARM','National','old_uEMMARM','RPi']
fileTypes = [('.emm files', '*.emm'), ('text files', '*.txt'), ('CSV files', '*.csv')]

#File naming option
resultFilePreffix = 'RaspberryPi_EMM-ARM' #Preffix to be added to the beginning of the result file
headerResultFiles = '1ST BATCH PAPER \n'

time_break=600


#=======================================================================================
# Selection options
selectedProcessing = typeOfProcessing[0] #Options: ['singleFile','batchProcessing']
selectedSystem = typeOfSystem[3] #Options: ['uEMMARM','National','old_uEMMARM','RPi]
desiredChannel = 1 #Number of variables in the file. For example, National test files may have four channels

#Calibration factor of the sensor used to convert test file reading to g's
calibrationFactor = 1 #This is for National
#calibrationFactor = 1/2400 #This is for commercial ADXL335

# Sampling frequency
#samplingFrequencyOriginal=1200 #In Hz. Not used by uEMMARM, which is based on a fixed duration of measurement session, automatically obtained from test files
# Length of window for PSD/FFT
#nps = samplingFrequencyOriginal +1

#Filter configuration
#See documentation about the function
filterConfiguration = { 1: {'filter': 'detrend', 'type': 'linear'},
						3: {'filter': 'butterworth', 'order': 8, 'type': 'bandpass', 'frequencies': [10,250]}}

#Selection of modal identification methods
modalIdentificationMethod_SingleAnalysis = {'peak-picking':False, 'BFD':False, 'EFDD':False, 'SSI-COV':True}
modalIdentificationMethod_BatchAnalysis = {'peak-picking':True, 'BFD':False, 'EFDD':False, 'SSI-COV':True}

#Peak-picking method
#Interval for averaging around the peak for the average peak-picking method
intervalForAveragingHz = 0.1

#EFDD method
EFFDfint = [0.5, 1.5] #This specifies the interval around the peak frequency used to compute the autocorrelation function. Specified in % of the peak frequency
EFFDtint = np.array([0.05,1]) #This specifies initial and final time interval used to fit the theoretical autocorrelation function. Specified in seconds

#SSI Configuration
numModesToBeConsidered=3 # This variable specifies how many modes will be considered in SSI identification. For standard EMM-ARM, we need only 1 (allegedly if the experiment was correctly performed, it will refer to the first frequency)
i    = 60 # number of time lags. From Peeters: The number of block rows i of Href should be such that r i >= nmax, which is the maximum model order.
startingOrderNumber = 2 # The order number by which the SSI method will start
endOrderNumber = 60 # The order number by which the SSI method will end
incrementBetweenOrder = 1 # The increment between each iteration
refs = [0] # reference sensors (all three are ref.)
stabcrit = {'freq':0.03, 'damping': 0.15, 'mac':0.02}
# This is defined based on what is expectable from tests
valid_range={'freq': [10, 250], 'damping':[0.001, 0.10]}
# tolerance for SSI methods
tol = np.array(([0.001,10, 250],     # [allowed_variation, lower_bound, higher_bound]: eigenfrequencies
				[0.15,0.01,.055],   # [allowed_variation, lower_bound, higher_bound]: damping ratios
				[0.02,0,1]))        # [allowed_variation, lower_bound, higher_bound]: MAC

#Base configuration for plot
plotConfiguration={'typeForPSD': 'None', 'frequencyBand': [0,250], 'ylimForPSD': [1e-10,None], 'typeForEFDD': 'only_ANPSD', 'ylimForEFDD': [None,1], 'typeForBFD':True, 'typeForEFDD': 'Autocorrelation-SVD', 'fontSize': 10, 'fontName':'Times New Roman', 'figSize': (5,2), 'figSizeBFD': (5,5), 'figSizeEFDD': (5,5), 'typeForStabilizationDiagram': 'StabilizationPSD', 'figSizeStabilization': (7,7), 'dpi': 150}

selectedFileType = fileTypes[2]
selectedFileExtension = selectedFileType[1][-4:]

# Validation of numberChannels variable
if selectedSystem == "uEMMARM":
	desiredChannel = 1
elif selectedSystem == "old_uEMMARM":
	desiredChannel = 1

### =====================

i2c = busio.I2C(board.SCL, board.SDA)
dac = adafruit_mcp4725.MCP4725(i2c)
root = Tk()
root.withdraw()
root.attributes('-topmost', True)

led1=19
led2=6
led3=13
led4=16
list_led = [19,6,13,16]

relay1=26
relay2=20
relay3=21

acc1=2
acc2=1
acc3=2
acc4=3
adc_excitation=4 

GPIO.setmode(GPIO.BCM)
GPIO.setup(relay1, GPIO.OUT)
GPIO.setup(relay2, GPIO.OUT)
GPIO.setup(relay3, GPIO.OUT)
GPIO.setup(led1, GPIO.OUT)
GPIO.setup(led2, GPIO.OUT)
GPIO.setup(led3, GPIO.OUT)
GPIO.setup(led4, GPIO.OUT)

#GPIO.output(relay1,0)
#GPIO.output(relay2,0)
#GPIO.output(relay3,0)
GPIO.output(led1,1)
GPIO.output(led2,1)
GPIO.output(led3,1)
GPIO.output(led4,1)
GPIO.output(relay3,1)                                                            

count = 75
count_file = 75

Acceleration_list = []
Volatge_list = []


h = np.linspace(0,45,8000*45)                                                                                
w = chirp(h,f0=15, f1=250 , t1=45, method='linear')
x = (w+1)/2

def leds_off():
	for i in list_led:
		GPIO.output(i,0)

def led_sleep_routine():
	while True:
		for i in list_led:
			GPIO.output(i,1)
			time.sleep(0.1)
		for i in list_led:
			GPIO.output(i,0)
			time.sleep(0.1)

def led_sleep_routine_break(timeAtBreak):
	while (time.time()-timeAtBreak)<time_break:
		for i in list_led:
			GPIO.output(i,1)
			time.sleep(0.1)
		for i in list_led:
			GPIO.output(i,0)
			time.sleep(0.1)

def Write_data(count, folderPath_specimen_1, folderPath_specimen_2, count_file):
	date = time.strftime("%Y-%m-%d-%H-%M-%S")
	if count % 2 == 0:
		name_of_file = "mould_1_"
		path_file = folderPath_specimen_1
	else:
		name_of_file = "mould_3_"
		path_file = folderPath_specimen_2
	name = date + "_" + name_of_file + str(count_file)
	completename=os.path.join(path_file,name+".csv")
	with open(completename, 'a') as csvFile:
		fieldnames=['Acceleration','Voltage']
		writer = csv.DictWriter(csvFile, fieldnames=fieldnames)
		for i, j in zip(Acceleration_list, Volatge_list):
			writer.writerow({'Acceleration': i, 'Voltage': j})
	csvFile.close() 
	#data = pd.read_csv(completename, low_memory=False)
	#print("Total rows: {0}".format(len(data)))

def acceleration(start_time,duration_measurement, count):
	ADC = ADS1256.ADS1256()
	ADC.ADS1256_init()
	while (time.time() - start_time) < 0.5 :
		Output = ADC.ADS1256_GetChannalValue(acc1)*5.0/0x7fffff
		Voltage = 2*(ADC.ADS1256_GetChannalValue(adc_excitation)*5.0/0x7fffff)-5
		time.sleep(0.1)
	while (time.time()-start_time)<duration_measurement+1:
		if count % 2 ==0:
			Output = ADC.ADS1256_GetChannalValue(acc1)*5.0/0x7fffff
			Voltage = 2*(ADC.ADS1256_GetChannalValue(adc_excitation)*5.0/0x7fffff)-5
			Acceleration =round((Output-1.8675)/0.8675, 5)
			Acceleration_list.append(Acceleration)
			Volatge_list.append(Voltage)
		else:
			Output = ADC.ADS1256_GetChannalValue(acc2)*5.0/0x7fffff
			Voltage = 2*(ADC.ADS1256_GetChannalValue(adc_excitation)*5.0/0x7fffff)-5
			Acceleration =round((Output-1.8675)/0.8675, 5)
			Acceleration_list.append(Acceleration)
			Volatge_list.append(Voltage)         

def excitation(count):
	if count % 2 == 0:
		GPIO.output(relay1,1)
		GPIO.output(relay2,0)
		GPIO.output(relay3,1)
		GPIO.output(led1,1)
		#print("relay 1 on")
	else:
		GPIO.output(relay1,0)
		GPIO.output(relay2,1)
		GPIO.output(relay3,0)
		GPIO.output(led2,1)
	time.sleep(0.5)
	while True:
		#print("debut excitation")
		start_time2 = time.time()
		for i in x:
			dac.normalized_value = i
		duration = time.time()-start_time2
		#print("duration excitation: ",duration, " s")
		time.sleep(0.5)

def modalID(accelerations, measurement_frequency,age, count):
	#Apply filter
	accelerationFiltered, samplingFrequencyFiltered  = auxEMMARM.filtering(accelerations, measurement_frequency, filterConfiguration)
	yk = MRPy(accelerationFiltered,samplingFrequencyFiltered)
	#Calculate PSD and plot it
	nps=measurement_frequency+1
	PSD = SSI.SDM(yk, nperseg=nps, plot=plotConfiguration, window='hann', nfft=2*nps) #Default is 50% overlap, not configurable in CESSI.py

	FPP, ignoreThis, PSDAveragedPeakIndex, ignoreThis=auxEMMARM.averagedPeakPickingMethod(PSD, intervalForAveragingHz, verbose="on")
	PSD.pki  = np.array([PSDAveragedPeakIndex], dtype=int)
	PSD.MGi  = np.array([0], dtype=int)
	PSD.svi  = np.array([0], dtype=int)
	PSD.fint = np.array([EFFDfint[0]*FPP, EFFDfint[1]*FPP])
	PSD.tint = EFFDtint
	#Plot the ANPSD. This will populate PSD variable with ANPSD and pki variables
	PSD = SSI.ANPSD_from_SDM(PSD,plot=plotConfiguration, mode="batch")

	yk = SSI.rearrange_data(yk,refs) 
	FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL = SSI.SSI_COV_iterator(yk,i,startingOrderNumber, endOrderNumber, incrementBetweenOrder, plot=False)
	stableModes = SSI.stabilization_diagram(FSSI_MODEL,ZSSI_MODEL,VSSI_MODEL, tol=tol, plot=plotConfiguration, PSD=PSD, verbose=False)
	FSSI, ZSSI, VSSI, numStablePoles = SSI.stable_modes(FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL, stableModes, tol=0.01, spo=10, verbose=True)
	###5.4.2) ORGANIZE THE RESULTS
	###DESCRIPTION: the previous function was implemented for the general purpose of CESSI.py library. This will organize the results for EMM-ARM purpose
	FSSI_CLUSTER=np.zeros(numModesToBeConsidered)
	ZSSI_CLUSTER=np.zeros(numModesToBeConsidered)
	numStablePoles_CLUSTER=np.zeros(numModesToBeConsidered)
	eigenfrequenciesIndices = np.flip(np.argsort(numStablePoles))
	if FSSI.size == 0:
		print('No sufficiently large frequency clusters could be identified') 
	else:
		for r, l in enumerate(np.take_along_axis(FSSI, eigenfrequenciesIndices, 0)[0:3]): FSSI_CLUSTER[r]=l
	if ZSSI.size == 0:
		print('No sufficiently large damping clusters could be identified') 
	else:
		for r, l in enumerate(np.take_along_axis(ZSSI, eigenfrequenciesIndices, 0)[0:3]): ZSSI_CLUSTER[r]=100*l
	if numStablePoles.size == 0:
		print('No sufficiently large clusters could be identified') 
	else:
		for r, l in enumerate(np.take_along_axis(numStablePoles, eigenfrequenciesIndices, 0)[0:3]): numStablePoles_CLUSTER[r]=int(l)
	print("Resonant frequency :", round(FSSI_CLUSTER[0],1), "Hz")
	print("Damping ratio :", round(ZSSI_CLUSTER[0],1), "%")

	#===================================================================================
	###                                     E-modulus
	#===================================================================================
	if n_specimen == 2 and count % 2  != 0: 
		frequencyEmptyTube = 144.5 #(Hz) #second specimen
		# Material properties
		externalTubeDiameter = 50.3/1000 #(m)
		internalTubeDiameter = 44.55/1000
		massOfEmptyTube = 332.74/1000 #(kg)
		totalLengthOfTube = 551/1000 #(m)
		freeLength = 500/1000
		massAtMidSpan = (7.51+5.8+9.15+26.4)/1000
		massOfFilledTube = 2055.5/1000 #(Kg)
	else:
		frequencyEmptyTube = 154.2 #(Hz) #1er specimen
		# Material properties
		externalTubeDiameter = 50.2/1000 #(m)
		internalTubeDiameter = 44.2/1000
		massOfEmptyTube = 358.71/1000 #(kg)
		totalLengthOfTube = 549/1000 #(m)
		freeLength = 500/1000
		massAtMidSpan = (7.51+5.8+9.15+26.4)/1000
		massOfFilledTube = 2033.9/1000 #(Kg)
	#print("frequency tube", frequencyEmptyTube)
	#Compute geometrical and mechanical properties of the mould tube
	tubeMomentOfInertia = np.pi*((externalTubeDiameter**4)-(internalTubeDiameter**4))/64
	materialMomentOfInertia = np.pi*((internalTubeDiameter**4))/64
	tubeEmptyLinearMass = massOfEmptyTube/totalLengthOfTube
	tubeEmptyMassFreeLength = freeLength*tubeEmptyLinearMass
	#tubeModulusInitialGuess = ((freeLength)**3)*(massAtMidSpan+0.24*tubeEmptyMassFreeLength)*((frequencyEmptyTube*2*np.pi)**2)/(3*tubeMomentOfInertia) #Modulus given in Pa
	#Estimate the E-modulus of the mould tube from the transcental equation
	tubeModulusInitialGuess = (((np.pi*frequencyEmptyTube)**2)*(freeLength**3)*(massAtMidSpan+0.49*tubeEmptyMassFreeLength))/(12*tubeMomentOfInertia)
	tubeFlexuralStiffness = auxEMMARM.solveCantileverTranscendentalEquation(tubeModulusInitialGuess*tubeMomentOfInertia, frequencyEmptyTube, tubeEmptyLinearMass, freeLength, massAtMidSpan)
	tubeModulus = tubeFlexuralStiffness/tubeMomentOfInertia

	#Extract only the vibration frequencies
	Frequency = FSSI_CLUSTER[0]
	#Properties of filled tube
	tubeFullLinearMass = massOfFilledTube/totalLengthOfTube
	tubeFullMassFreeLength = freeLength*tubeFullLinearMass

	compositeFlexuralStiffnessInitialGuess = (((np.pi*Frequency)**2)*(freeLength**3)*(massAtMidSpan+0.49*tubeFullMassFreeLength))/(12)
	compositeFlexuralStiffness = auxEMMARM.solveCantileverTranscendentalEquation(compositeFlexuralStiffnessInitialGuess, Frequency, tubeFullLinearMass, freeLength, massAtMidSpan)
	modulusElasticity = (compositeFlexuralStiffness-tubeModulus*tubeMomentOfInertia)/materialMomentOfInertia
	print("Tube E-Modulus=", round(float(tubeModulus)/1000000000,1), "GPa")
   # print("flexural strength initial guess =", round(compositeFlexuralStiffnessInitialGuess,1))
   # print("flexural strength =", round(float(compositeFlexuralStiffness),1))
	print("E modulus =", round(float(modulusElasticity)/1000000000,1), "GPa")
	print("================= END OF RESULTS FOR THIS SPECIMEN ==================== ")
	print("=======================================================================")
	if n_specimen == 2 and count % 2  != 0: 
		with open("Modal_ID_3.csv", 'a') as csvFile:
			fieldnames=['age','Frequency','E_Modulus']
			writer = csv.DictWriter(csvFile, fieldnames=fieldnames)
			writer.writerow({'age':age ,'Frequency': Frequency, 'E_Modulus': float(modulusElasticity)})
		csvFile.close() #second specimen
	else:
		with open("Modal_ID_1.csv", 'a') as csvFile:
			fieldnames=['age','Frequency','E_Modulus']
			writer = csv.DictWriter(csvFile, fieldnames=fieldnames)
			writer.writerow({'age':age ,'Frequency': Frequency, 'E_Modulus': round(float(modulusElasticity)/1000000000,3)})
		csvFile.close() #(Hz) #1er specimen
	
try:
	if __name__ == '__main__':
		l = Process(target=led_sleep_routine)
		l.start()
		folderPath_specimen_1 = filedialog.askdirectory(title="Select folder to save files of first specimen")
		time.sleep(0.5)
		folderPath_specimen_2 = filedialog.askdirectory(title="Select folder to save files of second specimen")
		time.sleep(0.5)
		n_specimen = simpledialog.askinteger(title="Number of specimen", prompt="How many specimen under test?")
		time.sleep(0.5)
		duration_measurement = simpledialog.askinteger(title="Measurement duration", prompt="Measurement duration (s)")
		l.terminate()
		testBeginning_time = time.time()
		while(1):
			leds_off()
			start_time = time.time()
			ExcitationProcess = Process(target = excitation, args=(count,))
			ExcitationProcess.start()
			MeasurementProcess = Process(target=acceleration(start_time,duration_measurement,count))
			MeasurementProcess.start()
			MeasurementProcess.join()
			ExcitationProcess.terminate()
			time_end = time.time()
			breakProcess = Process(target=led_sleep_routine_break, args=(time_end,))
			breakProcess.start() 
			age=(time.time()-testBeginning_time -0.5)/(3600*24) #days
			measurement_frequency=int(len(Acceleration_list)/ (time_end -1 -start_time))
			print('frequency measurement = ',measurement_frequency, " Hz")
			Write_data(count, folderPath_specimen_1, folderPath_specimen_2, count_file)
			print("age = ", age, "days")
			ModalIDProcess = Process(target=modalID(Acceleration_list,measurement_frequency,age, count))
			ModalIDProcess.start()
			Acceleration_list.clear()
			Volatge_list.clear()
			GPIO.output(relay1,0)
			GPIO.output(relay2,0)
			breakProcess.join()
			if n_specimen == 2 and count % 2  != 0:
				count_file +=1
			if n_specimen==2:
				count += 1

except KeyboardInterrupt:
	dac.raw_value = 0
	GPIO.cleanup()
	print("exit")
