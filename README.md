# Cost-effective-EMM-ARM
The Elastic Modulus Measurement through Ambient Response Method (EMM-ARM) is designed to continuously monitor the elastic modulus of hardening construction materials such as concrete, cement paste, mortars, stabilized soils, and epoxy resin.
By using a Raspberry Pi for the computing device and cost-effective electronic components, the cost of the system was one-twentieth of the traditional one, without compromising the measurement reliability.

Regarding the codes, several files were are necessary:
- ADS1256.py contains the functions to register parameters and manage measurements with the ADC.
- auxFunctions_postProcessEMMARM contains functions for the modal identification of the resonant frequencies.
- CESSIPy_modRenan contains functions associated with calculations of the resonant frequencies.
- config_ADC contains the configuration of the ADC.
- Experiment_code represents the main code, calling all the other presented codes.
- MRPy also contains functions associated with calculations of the resonant frequencies.

A 3D printed box was designed in order to host all the components of the device (STEP files available)
<img width="323" alt="image" src="https://user-images.githubusercontent.com/97248962/233833350-5f6a63ad-1edb-40da-8da8-1ee0c3c2cf22.png">

Associated publication: https://www.mdpi.com/2075-5309/13/5/1117
