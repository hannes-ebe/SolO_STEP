########################################################
###Testing of average pitch angle for different pixel###
########################################################

import numpy as np
import matplotlib.pyplot as plt
from virtual_detector import *

flow_vector = np.array([[-0.8412,  0.4396,  0.3149], [-0.8743,  0.457,   0.1635], [-0.8862,  0.4632, -0.    ], [-0.8743,  0.457,  -0.1635], [-0.8412,  0.4396, -0.315 ], [-0.7775,  0.5444,  0.3149], [-0.8082,  0.5658,  0.1635], [-0.8191,  0.5736,  0.    ], [-0.8082,  0.5659, -0.1634], [-0.7775,  0.5444, -0.3149], [-0.7008,  0.6401,  0.3149], [-0.7284,  0.6653,  0.1634], [-0.7384,  0.6744, -0.    ], [-0.7285,  0.6653, -0.1635], [-0.7008,  0.6401, -0.315 ]])

def pw(self,flow,B):
    '''Übergebe den particle flow-Vektor als Geschwindigkeit und den Magnetfeldvektor (am besten in SRF) und berechne die Pitchwinkel über das Skalarprodukt.'''
    len_flow = np.sqrt(flow[0]**2 + flow[1]**2 + flow[2]**2)
    len_B = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    argument = (flow[0]*B[0] + flow[1]*B[1] + flow[2]*B[2])/len_flow/len_B
    result = np.arccos(argument)
    return result

### Skript von Lars nutzen und dann mal vergleichen mit alter Berechnung

det = VDSTEP()

# plot_cosmu_cov(det)
# plot_pic_frac(det)
plot_FoV(det)
# plot_illumination(det)
# plot_geomf(det)