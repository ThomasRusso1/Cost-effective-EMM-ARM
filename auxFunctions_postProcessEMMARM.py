# Importing libraries
import numpy   as     np

#Module for selecting files
from tkinter import filedialog

#Module to deal with files from the system
import os
import pandas as pd

#Modal library modules
import CESSIPy_modRenan as SSI 
from MRPy import MRPy #Library with modal analysis functions
from scipy.signal import detrend, welch, resample, decimate, find_peaks, butter, sosfilt

#Module to solve the transcendental function to estimate E-modulus
from scipy.optimize import fsolve

#Module for plotting
import matplotlib.pyplot as plt

def readSingleFile(pathForFile, selectedSystem, desiredChannel=1):
    """
    Read a single file selected in the dialog

    Parameters
    -------       
    pathForFile : str
        Complete path of the file
    selectedSystem: str
        String that defines the type of system selected, so the file can be correctly read
    desiredChannel: byte, optional
        Desired channel in the file. For example, in case of National system, test files can have up to 4 channels for example, so we have to select one.
        Default case is 1, since some test files/systems have only one channel per file

    Returns
    -------    
    acceleration : nparray
        nparray containing the acceleration read.
    """ 
    #Select the type of system
    if selectedSystem == "National":
        acceleration = pd.read_table(pathForFile.name, names=["accel_"+str(i+1) for i in range(0,2)]) #range(0,2) because national will always provide 2 valued accel files
    elif selectedSystem == "old_uEMMARM":
        acceleration = pd.read_table(pathForFile.name, names=["accel_0"])
    elif selectedSystem == "uEMMARM":
        # Create a dtype with the binary data format and the desired column names
        dt = np.dtype([("accel_1", 'i2')])
        data = np.fromfile(pathForFile.name, dtype=dt)
        acceleration = pd.DataFrame(data)
    elif selectedSystem == "RPi":
        acceleration = pd.read_csv(pathForFile.name)
    else:
        raise Exception('ERROR: Selected system is not implemented in this version')

    acceleration=acceleration.to_numpy().T[desiredChannel-1]
    return acceleration #Returns a np.array

def readBatchFile(folderPath, files, selectedSystem, desiredChannel=1):
    """
    Read a single file from the directory fo batch reading

    Parameters
    -------       
    folderPath : str
        Complete path of the file
    files: str
        Name of the file to be read
    selectedSystem: str
        String that defines the type of system selected, so the file can be correctly read
    desiredChannel: byte, optional
        Desired channel in the file. For example, in case of National system, test files can have up to 4 channels for example, so we have to select one.
        Default case is 1, since some test files/systems have only one channel per file

    Returns
    -------    
    acceleration : nparray
        nparray containing the acceleration read.
    """ 
    #Select the type of system
    if selectedSystem == "National":
        acceleration = pd.read_table(folderPath+"/"+files, names=["accel_"+str(i+1) for i in range(0,2)]) #range(0,2) because national will always provide 2 valued accel files
    elif selectedSystem == "old_uEMMARM":
        acceleration = pd.read_table(folderPath+"/"+files, names=["accel_0"])
    elif selectedSystem == "uEMMARM":
        # Create a dtype with the binary data format and the desired column names
        dt = np.dtype([("accel_1", 'i2')])
        data = np.fromfile(folderPath+"/"+files, dtype=dt)
        acceleration = pd.DataFrame(data)
    elif selectedSystem == "RPi":
        acceleration = pd.read_csv(folderPath+"/"+files)
    else:
        raise Exception('ERROR: Selected system is not implemented in this version')
        
    acceleration=acceleration.to_numpy().T[desiredChannel-1]
    return acceleration #Returns a np.array

def getAgeAtMeasurementBatchFile(folderPath, files, firstMeasurementFile, selectedSystem):
    """
    Retrieve the age, in seconds, of the material at the instant of measurement. It considers only direct age (no corrections such as maturity correction)

    Parameters
    -------       
    folderPath : str
        Complete path of the file
    files: str
        Name of the file to be read
    firstMeasurementFile: str
        Name of the first measurement file
    selectedSystem: str
        String that defines the type of system selected, so the file can be correctly read

    Returns
    -------    
    ageAtMeasurement : int
        Age in seconds at the instant of start of the current measurement session under consideration.
    """ 
    from datetime import datetime
    from datetime import timedelta, date

    #Select the type of system
    if selectedSystem == "National":
        #The isntant of measurement for the National system is stored in the file name.
        currentSeconds = 10*int(files[-6]) + int(files[-5])
        currentMinutes = 10*int(files[-9]) + int(files[-8])
        currentHours = 10*int(files[-12]) + int(files[-11])
        currentDay = 10*int(files[-15]) + int(files[-14])
        currentMonth = 10*int(files[-18]) + int(files[-17])
        currentYear = 2000+10*int(files[-21]) + int(files[-20])
        currentTime = datetime(year=currentYear, month=currentMonth, day=currentDay,
                                hour=currentHours, minute=currentMinutes, second=currentSeconds, microsecond=0, tzinfo=None, fold=0)

        #The isntant of measurement for the National system is stored in the file name.
        initialSeconds = 10*int(firstMeasurementFile[-6]) + int(firstMeasurementFile[-5])
        initialMinutes = 10*int(firstMeasurementFile[-9]) + int(firstMeasurementFile[-8])
        initialHours = 10*int(firstMeasurementFile[-12]) + int(firstMeasurementFile[-11])
        initialDay = 10*int(firstMeasurementFile[-15]) + int(firstMeasurementFile[-14])
        initialMonth = 10*int(firstMeasurementFile[-18]) + int(firstMeasurementFile[-17])
        initialYear = 2000+10*int(firstMeasurementFile[-21]) + int(firstMeasurementFile[-20])
        initialTime = datetime(year=initialYear, month=initialMonth, day=initialDay,
                                hour=initialHours, minute=initialMinutes, second=initialSeconds, microsecond=0, tzinfo=None, fold=0)
    elif selectedSystem == "old_uEMMARM":
        a=1
        # TODO: Implement extracting date time from files obtained from the old EMM-ARM system
        # acceleration = pd.read_table(folderPath+"/"+files, names=["accel_0"]) #Gives 25.848842 Hz in National
    elif selectedSystem == "uEMMARM":
        # Create a dtype with the binary data format and the desired column names
        # dt = np.dtype([("accel_1", 'i2')])
        # data = np.fromfile(folderPath+"/"+files, dtype=dt)
        # acceleration = pd.DataFrame(data)
        timeData = pd.read_table(folderPath+"/"+files[0:-3]+"txt")
        #The isntant of measurement for the National system is stored in the file name.
        currentSeconds = int(timeData._values[5][0])
        currentMinutes = int(timeData._values[4][0])
        currentHours = int(timeData._values[3][0])
        currentDay = int(timeData._values[2][0])
        currentMonth = int(timeData._values[1][0])
        currentYear = 2000+int(timeData._values[0][0])
        currentTime = datetime(year=currentYear, month=currentMonth, day=currentDay,
                                hour=currentHours, minute=currentMinutes, second=currentSeconds, microsecond=0, tzinfo=None, fold=0)

        #The isntant of measurement for the National system is stored in the file name.
        timeData = pd.read_table(folderPath+"/"+firstMeasurementFile[0:-3]+"txt")
        initialSeconds = int(timeData._values[5][0])
        initialMinutes = int(timeData._values[4][0])
        initialHours = int(timeData._values[3][0])
        initialDay = int(timeData._values[2][0])
        initialMonth = int(timeData._values[1][0])
        initialYear = 2000+int(timeData._values[0][0])
        initialTime = datetime(year=initialYear, month=initialMonth, day=initialDay,
                                hour=initialHours, minute=initialMinutes, second=initialSeconds, microsecond=0, tzinfo=None, fold=0)

        # Compute the difference in time
        # Account for delay in the beggining of the test
    
    # Compute the difference in time
    # Account for delay in the beggining of the test
    ageAtMeasurement = currentTime - initialTime
    return ageAtMeasurement.total_seconds() #Returns a np.array

def convertToG(accelerationDigital,calibrationFactor):
    """
    Convert from digits (raw data from data-acquisition system) to acceleration in g's, by multiplication by a calibration factor

    Parameters
    -------       
    accelerationDigital : nparray
        Numpy array with the series of acceleration values.
    calibrationFactor: float
        Calibration factor that converts from digits (digital units) to g force (g=~9.81 m/s2)

    Returns
    -------    
    accelerationDigital*calibrationFactor : nparray
        Numpy array with accelerations in g's (g=~9.81 m/s2)
    """ 
    return accelerationDigital*calibrationFactor

def getSamplingFrequency_uEMMARM(folderPath, files, numberOfSamplingPoints):
    """
    Read a the .txt file that accompanies the .emm files in a uEMMARM system test to obtain sampling frequency.
    The uEMMARM system works with a fixed duration of measurement session, which is registered at the companion .txt file.
    The sampling frequency is also preconfigured in the system, but real sampling frequency may suffer minor variations within each test due to system instability.
    Thus, the true sampling frequency of each session is computed as the number of sampling points divided by the duration of the measurement session

    Parameters
    -------       
    folderPath : str
        Complete path of the file
    files: str
        Name of the file to be read
    numberOfSamplingPoints: int
        Number of sampling points in the measurement session. Obtained from the previous read .emm file

    Returns
    -------    
    samplingFrequency : float
        Sampling frequency computed as explained in the description of this function
    """ 
    #Read .txt file associated to the current .emm file
    with open(folderPath+"/"+files[:-4]+".txt") as f:
        lines = f.readlines()
    #Extract the duration:
    sessionDuration = int(lines[13]) #Use 10 for uEMMARM v0.2, and 13 for uEMMARM v0.1
    samplingFrequency = numberOfSamplingPoints/sessionDuration
    return samplingFrequency #Returns a np.array


def filtering(acceleration, samplingFrequency, filterConfiguration):
    """
    Applies the filtering specified in filterConfiguration parameter, which shall be a dictionary. If multiple filters are specified, they will be applied in order.

    Parameters
    -------       
    acceleration: nparray
        Numpy array with the series of acceleration values.
    samplingFrequency: scalar
        Scalar specifying the sampling frequency of the acceleration time series.
    filterConfiguration: nested dictionary
        Nested dictionary that will list, in order of application, the filters to be applied to the acceleration series.
        The filters are simplified versions of scipy.signal functions.
        Currently, the following filter are supported, which require the following configuration.
        1: {'filter': 'detrend', 'type': 'linear' or 'constant'}. For further info, see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.detrend.html#scipy.signal.detrend
        2: {'filter': 'decimation', 'decimationFactor': positive integer - e.g. 10}. For further info, see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.decimate.html#scipy.signal.decimate
        3: {'filter': 'butterworth', 'order': positive integer - e.g. 8, 'type': type of filter - 'highpass' or 'lowpass' or 'bandpass' or 'bandstop', 'frequencies': frequencies of the filter - if highpass or lowpass it is a scalar - if bandpass or bandstop it is a list to specificy loww and high frequencies of the band, 'samplingFrequency': the sampling frequency of the signal in Hz}}. For further info, see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html#scipy.signal.butter

    Returns
    -------    
    acceleration : nparray
        nparray containing the filtered acceleration series.
    """ 
        #Apply as many filters as it was passed in filterConfiguration
    for filter in filterConfiguration:
        currentFilter = filterConfiguration[filter]
        if currentFilter['filter']=='detrend':
            #Apply detrending
            acceleration = detrend(acceleration, type = currentFilter['type'])
        elif currentFilter['filter']=='decimation':
            #Apply the decimation filter
            acceleration = decimate(acceleration,currentFilter['decimationFactor'])
            samplingFrequency = samplingFrequency/currentFilter['decimationFactor']
        elif currentFilter['filter']=='butterworth':
            if (currentFilter['type']=='lowpass') or (currentFilter['type']=='highpass'):
                currentFilter['frequencies']=currentFilter['frequencies'][0]
            designedButterworthFilter = butter(currentFilter['order'],currentFilter['frequencies'],currentFilter['type'], fs=samplingFrequency, output='sos')
            #designedButterworthFilter = butter(8, 15, 'hp', fs=1000, output='sos')
            acceleration = sosfilt(designedButterworthFilter, acceleration)
        else:
            raise Exception("The filter specified is not currently supported. See documentation for supported configuration")

    return acceleration, samplingFrequency

def plotAccelerationTimeSeries(accelerationData, plot={'fontSize': 15, 'fontName':'Times New Roman', 'figSize': (5,2), 'dpi': 150}):
    """
    Function to make a standardized plot of acceleration time series. 
    Multiple acceleration time series may be plotted in a single graph, if desired.
    The data to be plotted is specified by accelerationData, which main contain a nested list with data and metadata from each acceleration time series. 

    Parameters
    -------       
    accelerationData : nested list
        Nested list containing, in each row, the data and metadata from an acceleration time series.
        The format of accelerationData is:
            accelerationData = [[accelerationTimeSeries1[1:n],samplingFrequency1,label1],
                                [accelerationTimeSeries2[1:n],samplingFrequency2,label2],
                                ...
                                [accelerationTimeSeriesN[1:n],samplingFrequencyN,labelN],]
        In which:
            accelerationTimeSeriesN[1:n] is a 1-column nparray containing the acceleration values.
            samplingFrequencyN is a scalar specifying the sampling frequency of acceleratiomTimeSeriesN, so time data can be reconstructed.
            labelN is the label to be used in the plot to identify accelerationTimeSeriesN.
    plot : dictionary, optional #Editted EMM-ARM
        It has the following format:
            plot={'fontSize': 15, 'fontName':'Times New Roman', 'figSize': (5,2), 'dpi': 150}
        In which:
            'fontSize' is a scalar and specifies the base font size of the plot
            'fontName' is a str and specifies the font type of the plot
            'figSize' is a tuple (width, height) and specifies the size of the figure
            'dpi' is a scalar and specifies the DPI of the figure
    Returns
    -------    
    Nothing.

    """ 

    plt.figure(figsize=plot['figSize'], dpi=plot['dpi'])

    for accelerationSeries in accelerationData:
        time=np.arange(0,len(accelerationSeries[0])/accelerationSeries[1],1/accelerationSeries[1])
        plt.plot(time,accelerationSeries[0], label=accelerationSeries[2])

    plt.xlabel("Time (s)", size=plot['fontSize'], fontname=plot['fontName'])
    plt.ylabel("Acceleration (g)", size=plot['fontSize'], fontname=plot['fontName'])
    plt.legend()
    plt.show()

def averagedPeakPickingMethod(PSD, intervalForAveragingInHz, verbose=False):
    #TODO: Implement allowing identification of more than 1 peak
    """
    This method adapts a "crude" version of the peak-picking method for frequency identification by considering a pondering averaged with the PSD intensities around the PSD peak.

    Estimate the first peak in the PSD amplitude and the associated frequency, called yMaxFrequency. The index associated to this frequency is yMaxPeakIndex.
    Around such peak, an average is taken on the intervalForAveragingInHz (ex.: [refFrequency-intervalForAveragingInHz,refFrequency+intervalForAveragingInHz000]), by taking the values of the PSD as pondering factors. The results of the averaging is called averagedFrequency
    
    Not necessairy averagedFrequency exists in the PSD, as it is the product of an average with the frequency values of the PSD. So, the closest frequency of the PSD to the averagedFrequency is also found. This frequency is called PSDAveragedFrequency and the associated index is PSDAveragedPeakIndex. They are not a fruit of the peak-picking method, but useful to start other frequency identification methods, such as EFDD, that require an initial input/estimate of the frequency in the PSD series.

    Parameters
    -------       
    PSD : auxclass_like
        Auxclass object that contains the attributes f and pki.
    intervalForAveragingHz: float
        Defines the value, in Hz, at each side of the peak, used to average and find the natural frequency
    verbose: bool, optional.
        Defines if verbose mode is on, so to print the results of the identification metho

    Returns
    -------    
    averagedFrequency: scalar
        The frequency identified with the method called "Averaged Peak-Picking"
    PSDAveragedFrequency: scalar
        The frequency in the PSD closest to averagedFrequency.
    PSDAveragedPeakIndex: scalar
        The index associated to the PSDAveragedFrequency in the PSD series.
    yMaxPeakIndex: scalar
        The index associated to yMaxFrequency. This is the first peak identified in the PSD.

    """ 
    # Find the peak
    yPeaksIndex, _ = find_peaks(abs(PSD)[0][0])
    yMaxPeakIndex = yPeaksIndex[np.argmax(abs(PSD)[0][0][yPeaksIndex])]

    # Find the index from the maximum peak, and lower and higher boundary for averaging
    indexLowerAvgingBoundary = (np.abs(PSD.f - (PSD.f[yMaxPeakIndex]-intervalForAveragingInHz))).argmin()
    indexHigherAvgingBoundary = (np.abs(PSD.f - (PSD.f[yMaxPeakIndex]+intervalForAveragingInHz))).argmin()

    #Find averaged eigenfrequency around the peak
    rangeOfInterest=range(indexLowerAvgingBoundary, indexHigherAvgingBoundary+1, 1)

    #Prepare vectors to caculate average with vectorized functions
    vectorPSD=np.array(np.abs(PSD)[0][0][rangeOfInterest])
    vectorPSD=vectorPSD.reshape(len(vectorPSD),1) #Reshape 0D vector to a 2D vector
    vectorFrequency=np.array(np.abs(PSD.f[rangeOfInterest]))
    vectorFrequency=vectorFrequency.reshape(len(vectorFrequency),1) #Reshape 0D vector to a 2D vector

    #Use vector product to calculate averaged frequency around the peak
    averagedFrequency=np.dot(vectorPSD.transpose(),vectorFrequency)/np.sum(vectorPSD)

    #Find the closest index to averagedFrequency
    PSDAveragedPeakIndex = (np.abs(PSD.f - averagedFrequency)).argmin()
    PSDAveragedFrequency = PSD.f[PSDAveragedPeakIndex]

    if verbose == True:
        print("=================================================================================")
        print("RESULTS FROM AVERAGED PEAK-PICKING METHOD")
        print("Peak selected as *reference* peak for the averaged peak-picking method:")
        print("{:.3f} Hz".format(PSD.f[yMaxPeakIndex]))
        print("Considering an interval around *reference* peak of {:.3f} Hz.".format(intervalForAveragingInHz))
        print("Averaged peak-picking estimated frequency:")
        print("{:.3f} Hz".format(averagedFrequency[0][0]))
        print("END OF RESULTS FROM AVERAGED PEAK-PICKING METHOD")
        print("=================================================================================")

    #return PSDaveragedPeakIndex, averagedFrequency, yMaxPeakIndex
    return averagedFrequency, PSDAveragedFrequency, PSDAveragedPeakIndex, yMaxPeakIndex

def solveCantileverTranscendentalEquation(initialGuess, vibrationFrequency, linearMass, freeLength, tipMass):
    """
    This method numerically solves the transcendental equation of a cantilever beam under free vibration with a concentrated mass at its free tip, outputing the flexural stiffness (EI) of the beam

    Parameters
    -------       
    initialGuess: float
        E-modulus initial guess, so to start the numerical process.
    vibrationFrequency: float
        Vibration frequency of the beam, in Hz
    linearMass: float
        Mass of the tube filled with the material, in kg/m
    freeLength: float
        Free length of the cantilever beam, in meters
    tipMass: float
        Total mass at the free tip of the cantilever beam, in kg

    Returns
    -------    
    flexuralStiffness: float
        The flexural stiffness (EI) of the beam. If it is a composite beam, it is the composite flexural stiffness (considering the two materials as a perfect composite section)
    """ 

    #Compute the natural frequency in rad/s and store important variables in easy-to-read codes
    ω = 2*np.pi*vibrationFrequency
    L = (freeLength/2)
    mL = linearMass
    mT = tipMass
  #  a = (((mL*(ω**2))/EI)**(1/4))
    k=(8.72*(10**12)) #(N/m)
    h=9999
    #Define the transcendental function structure
  #  f = lambda EI: ((((ω**2)*mL/EI)**(1/4))**3)*(np.cosh((((ω**2)*mL/EI)**(1/4))*L)*np.cos((((ω**2)*mL/EI)**(1/4))*L)+1)+(ω*ω*mT/EI)*(np.cos((((ω**2)*mL/EI)**(1/4))*L)*np.sinh((((ω**2)*mL/EI)**(1/4))*L)-np.cosh((((ω**2)*mL/EI)**(1/4))*L)*np.sin((((ω**2)*mL/EI)**(1/4))*L))
   # f = lambda EI: (1/(2*k))*(EI*(a**3)*mT*(1+(np.cosh(a*L)**2)-(np.sinh(a*L)**2)+(2*np.cos(a*L)*np.cosh(a*L)))+2*EI*(a**6)*((np.sin(a*L)*np.cosh(a*L))+(np.cos(a*L)*np.sinh(a*L)))+4*EI*(a**3)*k*np.cos(a*L)*np.cosh(a*L)+2*k*(ω**2)*mT*((np.cos(a*L)*np.sinh(a*L))+(np.cosh(a*L)*np.sin(a*L))))
    
    f = lambda EI: (1/(2*k))*(EI*((((mL*(ω**2))/EI)**(1/4))**3)*mT*(1+(np.cosh((((mL*(ω**2))/EI)**(1/4))*L)**2)-(np.sinh((((mL*(ω**2))/EI)**(1/4))*L)**2)+(2*np.cos((((mL*(ω**2))/EI)**(1/4))*L)*np.cosh((((mL*(ω**2))/EI)**(1/4))*L)))+2*EI*((((mL*(ω**2))/EI)**(1/4))**6)*((np.sin((((mL*(ω**2))/EI)**(1/4))*L)*np.cosh((((mL*(ω**2))/EI)**(1/4))*L))+(np.cos((((mL*(ω**2))/EI)**(1/4))*L)*np.sinh((((mL*(ω**2))/EI)**(1/4))*L)))+4*EI*((((mL*(ω**2))/EI)**(1/4))**3)*k*np.cos((((mL*(ω**2))/EI)**(1/4))*L)*np.cosh((((mL*(ω**2))/EI)**(1/4))*L)+2*k*(ω**2)*mT*((np.cos((((mL*(ω**2))/EI)**(1/4))*L)*np.sinh((((mL*(ω**2))/EI)**(1/4))*L))+(np.cosh((((mL*(ω**2))/EI)**(1/4))*L)*np.sin((((mL*(ω**2))/EI)**(1/4))*L))))

  #  f = lambda EI: (a**3)*(1+(np.cosh(a*L)*np.cos(a*L)))+(((ω**2)*mT)/EI)*((np.cos(a*L)*np.sinh(a*L))-(np.cosh(a*L)*np.sin(a*L)))+((EI*(a**4))/k)*((np.cos(a*L)*np.sinh(a*L))-(np.cosh(a*L)*np.sin(a*L)))-2*((a*mT*ω*ω)/k)*np.sin(a*L)*np.sinh(a*L)


    #Solve the transcendental equation
    flexuralStiffness = fsolve(f, initialGuess)

    return flexuralStiffness