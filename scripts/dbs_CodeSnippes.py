import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator 
import seaborn as sns
import re
import pyqtgraph as pg
from scipy import integrate
import pandas as pd
import numpy as np
import os
theme_bw = "/Users/au652733/Python/theme_bw.mplstyle"
plt.style.use(theme_bw)

# global parameter
dcolor = dict({'pH': '#1CC49A', 'sig pH': '#4B5258', 'para2': '#196E94', 'sig para2': '#314945',
               'para3': '#DCA744', 'sig para3': '#89621A', 'total para': '#A86349', 'background': '#383e42', 
               'font': 'white'})

ls = dict({'target': '-.', 'simulation': ':'})
save_type = ['jpg']
fs, fs_font = 10, 12
sns.set_context('paper', font_scale=1.)
pg.setConfigOption('background', dcolor['background'])
pg.setConfigOption('foreground', 'white')

# known acud/base systems
systems = dict({'TDS': dict({'base': ['HS-', 'HS', 'hs-', 'hs'], 'acid': ['H2S', 'h2s']}),
                'TAN': dict({'base': ['NH3', 'nh3'], 'acid': ['NH4+', 'NH4', 'nh4+', 'nh4']})})

# convert everything to seconds
dic_tunits2sec = dict({'s': 1, 'ms': 10**-3, 'µs': 10**-6, 'min': 60, 'h': 60*60})

# global parameter for electrochemical sensor calibration
n = 1
R = 8.314                    # [J/mol*K]
F = 96485                    # [C/mol]
s_nernst = 2.303
Tconvert = 273.15
temp_std = 25

# fixed parameter depending on literature research
ph_res = 1e-3                # resolution of the pH sensor
ph_E0 = 0.22                 # standard potential of the reference electrode of the pH meter; V
tsteps = 1e-3                # time steps for pH and measuring sensor (theory); seconds

# electrochemical sensor · 2
conc_calib = 1               # reaction quotient for calibrating the measuring sensor
sens_E = 0.09                # potential of the measuring sensor at certain concentration (conc_calib); V
sens_res = 1e-9              # resolution of the measuring sensor

# others
_integral_counter = 0        # helper scalar counting the integration
ls_lines = list()            # list to collect integration boundaries in pyqtgraph
dres, df_res = dict(), None  # simulation results

# get the local path for relative directories
loc_path = os.getcwd() 


# ---------------------------------------------------------------------------------------------------------------
def plot_phSensor(df_res, fgs=(5, 3.5)):
    fig, ax_ph = plt.subplots(figsize=fgs)
    ax_ph.set_xlabel('Time, sec', fontsize=fs)
    ax_ph.set_ylabel('pH value', fontsize=fs)
    ax_ph.set_title('pH sensor · target vs response signal', fontsize=fs)

    # pH sensor response curve
    ax_ph.plot(df_res['target pH'], color=dcolor['pH'], label='pH target')
    ax_ph.plot(df_res['pH calc'], color='k', ls='-.', label='pH calc')
    ax_ph.legend(fontsize=fs*0.75)

    # axis layout
    ax_ph.tick_params(labelsize=fs)

    for axis in ['top','bottom','left','right']:
        ax_ph.spines[axis].set_linewidth(1)
    fig.tight_layout()
    return fig, ax_ph


def plot_AcidBasePair(para2_grp, para3_grp, df_res, fgs=(5, 3.5)):
    # get the current unit
    unit = df_res.filter(like='target').filter(like='Acid').columns[0].split('_')[1].split(' ')[0]

    fig, ax_pair = plt.subplots(figsize=fgs)
    ax_pair1 = ax_pair.twinx()
    ax_pair.set_xlabel('Time, sec', fontsize=fs)
    ax_pair.set_ylabel(para2_grp + ', ' + unit, fontsize=fs, color=dcolor['para2'])
    ax_pair1.set_ylabel(para3_grp + ', ' + unit, fontsize=fs, color=dcolor['para3'])
    ax_pair.set_title('Acid-Base-Pair · target vs response signal', fontsize=fs)

    # acid/base pair sensor response curve

    lns1 = ax_pair.plot(df_res.filter(like='target', axis=1).filter(like='Acid'), lw=1., color=dcolor['para2'], 
                        label='acid target')
    lns2 = ax_pair.plot(df_res['Acid calc'], color='k', lw=1., ls='-.', label='acid calc')
    lns3 = ax_pair1.plot(df_res.filter(like='target', axis=1).filter(like='Base'), lw=1., color=dcolor['para3'], 
                         label='base target')
    lns4 = ax_pair1.plot(df_res['Base calc'], color='k', lw=1., ls=':', label='base calc')
    lns = lns1+lns2+lns3+lns4
    labs = [l.get_label() for l in lns]
    ax_pair.legend(lns, labs, loc='center left', fontsize=fs*0.75)

    # axis layout
    ax_pair.tick_params(labelsize=fs)
    ax_pair1.tick_params(labelsize=fs)

    ax_pair.spines['right'].set_visible(True)
    ax_pair1.spines['right'].set_visible(True)
    for axis in ['top','bottom','left','right']:
        ax_pair.spines[axis].set_linewidth(.75)
        ax_pair1.spines[axis].set_linewidth(.75) 
    fig.tight_layout()
    return fig, ax_pair, ax_pair1


def plot_totalParameter(total_para_grp, total_para, df_res, ylim=None, fgs=(5, 3.5)):
    # get the current unit
    unit = df_res.filter(like='target').filter(like='Acid').columns[0].split('_')[1].split(' ')[0]

    fig, ax_sum = plt.subplots(figsize=fgs)
    ax_sum.set_xlabel('Time, sec', fontsize=fs)
    ax_sum.set_ylabel(total_para_grp + ', ' + unit, fontsize=fs)
    ax_sum.set_title('Total parameter signal · target vs response', fontsize=fs)

    # total parameter response curve
    df_target = df_res.filter(like='target').filter(like='Sum')
    ax_sum.plot(df_target, color=dcolor['total para'], label=total_para +' target')
    ax_sum.plot(df_res['Sum calc'], color='k', ls='-.', lw=1., label=total_para+' calc')
    ax_sum.legend(fontsize=fs*0.75)

    # axis layout
    if ylim:
        if isinstance(ylim, tuple):
            ax_sum.set_ylim(ylim[0], ylim[1])
        elif '%' in ylim:
            s = float(ylim.split('%')[0].strip())/100
            ax_sum.set_ylim(df_target.min()[0]*(1-s), df_target.max()[0]*(1+s))
        elif 'pc' in ylim:
            s = float(ylim.split('p')[0].strip())/100
            ax_sum.set_ylim(df_target.min()[0]*(1-s), df_target.max()[0]*(1+s))

    ax_sum.tick_params(labelsize=fs)
    for axis in ['top','bottom','left','right']:
        ax_sum.spines[axis].set_linewidth(.75)
    ax_sum.xaxis.set_major_locator(MultipleLocator(25))
    ax_sum.xaxis.set_minor_locator(MultipleLocator(5))
    fig.tight_layout()
    return fig, ax_sum


# ---------------------------------------------------------------------------------------------------------------
def _linEdit2list(line):
    ls = list()
    if '.' in line and ',' in line:
        ls = [float(i) for i in line.split(',')]
    elif '.' in line:
        # double-check only one '.' in str
        if len(re.findall(r"\.", line)) == 1:
            ls.append(float(line))
        else:
            if len(line.split(' ')) != 1:
                ls = [float(i) for i in line.split(' ')]
            else:
                raise ValueError('Typo error. Please check your pH entry: ' + line)
    elif ',' in line:
        ls = [float(i) for i in line.split(',') if len(i.strip()) != 0]
    else:
        ls.append(float(line))
    return ls


def align_concentrations(ls_ph, ls_cSum_gmL):
    if len(ls_cSum_gmL) == len(ls_ph):
        pass
    elif len(ls_cSum_gmL) == 1 and len(ls_ph) > 1:
        ls_cSum_gmL = ls_cSum_gmL * len(ls_ph)
    elif len(ls_cSum_gmL) > 1 and len(ls_ph) == 1:
        ls_ph = ls_ph * len(ls_cSum_gmL)
    else:
        if max(len(ls_ph), len(ls_cSum_gmL)) % min(len(ls_ph), len(ls_cSum_gmL)) == 0:
            print('double the length of each other')
            if len(ls_ph) < len(ls_cSum_gmL):
                ls_ph = ls_ph * int(len(ls_cSum_gmL) / 2)
            else:
                ls_cSum_gmL = ls_cSum_gmL * int(len(ls_ph) / 2)
        else:
            if len(ls_ph) < len(ls_cSum_gmL):
                ls_ph.append(ls_ph[-1])
            else:
                ls_cSum_gmL.append(ls_cSum_gmL[-1])
    return ls_ph, ls_cSum_gmL


def _getAnalyteSystem(para2_name):
    para2_grp, para3_grp, total_para_grp, para2, para3, total_para, analyte = None, None, None, None, None, None, None
    # If an entirely new system is designed, add it as a tuple. The first parameter is the measured, the second the corresponding analyte and the third the total/sum parameter.
    if isinstance(para2_name, tuple):
        para2_grp, para2 = para2_name[0], para2_name[0]
        para3_grp, para3 = para2_name[1], para2_name[1]
        total_para_grp, total_para = para2_name[2], para2_name[2]  
        analyte = 'acid' if '+' in para2 or '-' in para3 else 'base'

    # otherwise use the pre-defined systems TAN or TDS
    else:
        if para2_name:
            # H2S/HS- known system
            if para2_name in systems['TDS']['base']:
                para2_grp, para3_grp, total_para_grp = 'HS-', 'H2S', 'TDS'
                para2, para3, total_para = 'HS$^-$', 'H$_2$S', '$TDS$'
                analyte = 'base'
            elif para2_name in systems['TDS']['acid']:
                para2_grp, para3_grp, total_para_grp = 'H2S',  'HS-', 'TDS'
                para2, para3, total_para = 'H$_2$S', 'HS$^-$', '$TDS$'
                analyte = 'acid'

            # NH3/NH4+ known system
            if para2_name in systems['TAN']['base']:
                para2_grp, para3_grp, total_para_grp = 'NH3', 'NH4+', 'TAN'
                para2, para3, total_para = 'NH$_3$', 'NH$_4^+$', '$TAN$'
                analyte = 'base'
            elif para2_name in systems['TAN']['acid']:
                para2_grp, para3_grp, total_para_grp = 'NH4+',  'NH3', 'TAN'
                para2, para3, total_para = 'NH$_4^+$', 'NH$_3$', '$TAN$'
                analyte = 'acid'

            # if unknown use "base / acid / sum" (the first one is the measured, the second one the calculated)
            if '/' in para2_name:
                system = para2_name.split('/')
                para2_grp, para3_grp, total_para_grp = system[0], system[1], system[-1]
                para2, para3, total_para = system[0], system[1], system[-1]
                analyte = system[0]
        else:
            para2_grp, para3_grp, total_para_grp = 'Base', 'Acid', 'Sum Parmeter'
            para2, para3, total_para = 'base', 'acid', 'sum parmeter'
            analyte = 'sum parmeter'
    return para2_grp, para3_grp, total_para_grp, para2, para3, total_para, analyte  


def adjust_tunit(para):
    if para[1] in dic_tunits2sec.keys():
        return (para[0]*dic_tunits2sec[para[1]], 's')


def parameter_prep(t_plateau, ls_ph, ls_total, pKa):
    # target fluctuation of analytes
    target_ph = target_fluctuation(ls_conc=ls_ph, tstart=0, tstop=t_plateau * 2, nP=1, analyte='pH')
    target_sum = target_fluctuation(ls_conc=ls_total, tstart=0, tstop=t_plateau * 2, nP=1, analyte='Sum')

    # target concentration para2/para3 based on target pH (para2... base, para3... acid)
    para2 = pd.DataFrame(henderson_Sum4Base(pKa=pKa, pH=target_ph['signal pH'].to_numpy(),
                                            c_sum=target_sum['signal Sum'].to_numpy()), index=target_ph.index,
                                            columns=['signal Base'])
    para3 = pd.DataFrame(henderson_Sum4Acid(pKa=pKa, pH=target_ph['signal pH'].to_numpy(),
                                            c_sum=target_sum['signal Sum'].to_numpy()), index=target_ph.index,
                                            columns=['signal Acid'])

    # calculate target concentration of individual NHx parameters
    df_res = pd.concat([target_ph, target_sum, para2, para3], axis=1)
    return df_res


def check_parameter(analyte, paraSum_conc_edit, ph_edit, para2_pka_edit, tsteady_edit, ph_t90_edit, para2_t90_edit):
    # check that units are corresponding
    if (tsteady_edit[1] == ph_t90_edit[1] == para2_t90_edit[1]) is False:
        tsteady_edit = adjust_tunit(tsteady_edit)
        ph_t90_edit, para2_t90_edit = adjust_tunit(ph_t90_edit), adjust_tunit(para2_t90_edit) 

    # get sum concentration(s) and pH value(s) and align list of fluctuation points
    ls_cSum_gmL = _linEdit2list(line=paraSum_conc_edit)
    ls_ph = _linEdit2list(line=ph_edit)
    ls_ph, ls_cSum_gmL = align_concentrations(ls_ph=ls_ph, ls_cSum_gmL=ls_cSum_gmL)
    
    # determine all required target concentrations including the sensor response
    df_target = parameter_prep(ls_ph=ls_ph, ls_total=ls_cSum_gmL, pKa=float(para2_pka_edit), t_plateau=float(tsteady_edit[0]))

    # calculate target concentrations of sensor 2
    if analyte in ['base', 'Base']:
        df_target_sens2 = pd.DataFrame(df_target['signal Base'])
    elif analyte in ['acid', 'Acid']:
        df_target_sens2 = pd.DataFrame(df_target['signal Acid'])
    else:
        df_target_sens2 = pd.DataFrame(df_target['signal Base'])
    ls_target_sens2 = list(df_target_sens2[df_target_sens2.columns[0]].drop_duplicates().values)

    # collect all relevant parameter
    sensor_ph = dict({'E0': ph_E0, 't90': float(ph_t90_edit[0]), 'resolution': ph_res, 'start pH': ls_ph[0], 
                      'time steps': tsteps, 'pH target': df_target['signal pH'], 'set values': ls_ph})

    sensor_para2 = dict({'analyte': analyte, 'time steps': tsteps, 'pKa': float(para2_pka_edit), 'E': sens_E, 
                         'resolution': sens_res, 't90': float(para2_t90_edit[0]), 
                         'Base target': df_target['signal Base'], 'Acid target': df_target['signal Acid'], 
                         'Sum target': df_target['signal Sum'], 'conc_calib': conc_calib,
                         'set values': ls_cSum_gmL, 'set sensor2': ls_target_sens2})
    para_meas = dict({'plateau time': float(tsteady_edit[0])})
    return sensor_ph, sensor_para2, para_meas


def _alignSensorSettings(df_target, sensor_ph, sensor2):
    # individual sensor calibration
    df_sigpH_mV, cplateaupH = calibSensor_pH(target_ph=df_target, sensor_ph=sensor_ph)
    df_sigPara_mV, cplateauPara, para_para = calibSensor_para(target_para=df_target, sensor2=sensor2)

    df_ = pd.concat([df_sigpH_mV, df_sigPara_mV], axis=1)
    df_.columns = ['Potential mV pH', 'Potential mV Sum']
    df_res = pd.concat([df_target, df_], axis=1)
    return df_res, cplateaupH, cplateauPara, para_para


def _all_plateaus_inlist(ls):
    ls_ = list()
    for xi in [x for i,x in enumerate(np.diff(ls)) if x != 0]:
        ls_ += (list(np.where(np.diff(ls) == xi)[0]))

    pos_jump = [0] + sorted(list(dict.fromkeys(ls_))) + [len(ls)]
    plateau = [np.mean(ls[pos_jump[en] + 1:pos_jump[en+1]]) for en in range(len(pos_jump)-1)]

    return plateau


# -----------------------------------------------------------------------------------------------------
# functions related to sensor response 
def _gompertz_curve_v1(x, t90, tau, pstart, s_diff, slope='increase'):
    """ Sensor response curve:
    :param x:       np.array; time range
    :param t90:     float; response time of the sensor in s
    :param tau:     float; resolution of the sensor
    :param pstart:  float; starting signal of the response curve
    :param s_diff:  float; difference in signal between end signal and start signal of the response curve
    :param slope:   str; 'increase' or 'decline' describing whether the response increases or decreases
    :return:
    """
    b = np.log(-np.log(tau))
    k = 1 / t90 * (b - np.log(np.log(1 / 0.9)))

    if slope == 'increase':
        y = s_diff * np.exp(-1 * np.exp(b - k*x)) + pstart
    elif slope == 'decline':
        y = pstart - s_diff * np.exp(-1 * np.exp(b - k*x))
    else:
        raise ValueError('Define slope of the response curve either as increase or decline')
    return y


def _sensor_response(func, cplateau, para_response, tplateau, cstart):
    # set the time range
    start, steps = 0, 0.05
    num = int((tplateau - start) / 0.05 + 1)
    sens_time = np.linspace(start, tplateau, num=num)

    # simulate for each signal level a respective sensor response
    c_apparent = cstart
    dsig, df_sensor = dict(), None
    for en, ctarget in enumerate(cplateau):
        sdiff = ctarget - c_apparent
        if en == 0:
            xnew = sens_time
        else:
            xnew = np.linspace((tplateau + steps) * en, (tplateau + steps) * en + tplateau, num=num)

        # actual sensor response function
        df_sig = pd.DataFrame(func(x=sens_time, s_diff=sdiff, pstart=c_apparent, **para_response), index=xnew,
                              columns=['potential mV'])

        # prepare for the next signal change in line
        start, c_apparent = df_sig['potential mV'].index[-1], df_sig['potential mV'].to_numpy()[-1]
        dsig[en] = df_sig
        df_sensor = pd.concat(dsig)
        df_sensor.index = [i[1] for i in df_sensor.index]

    return df_sensor


# ----------------------------------------------------------------
# corresponding acid/base equilibrium
def henderson_base(pH, c_acid, pKa=9.25):
    """ Returns the base concentration depending on the pH and the concentrations of acid """
    return c_acid * 10 ** (pH - pKa)


def henderson_acid(pH, c_base, pKa=9.25):
    """ Returns the acid concentration depending on the pH and the concentrations of the base """
    return c_base * 10 ** (pKa - pH)


def henderson_Sum4Base(pH, c_sum, pKa=9.25):
    """
    Returns the base concentration depending on the pH and the SUM concentration
    :param pH:
    :param c_sum:
    :param pKa:
    :return:
    """
    return c_sum / (1 + 10**(pKa - pH))


def henderson_Sum4Acid(pH, c_sum, pKa=9.25):
    """
    Returns the acid concentration depending on the pH and the SUM concentration
    :param pH:
    :param c_sum:
    :param pKa:
    :return:
    """
    return c_sum / (1 + 10**(pH - pKa))


# -----------------------------------------------------------------------------------------------------
# electrochemical sensor calculations (related to pH sensor)
def Nernst_equation(ls_ph, E0, T=25):
    """ convert apparent pH into electrical potential. The logarithmic of the activity equals the pH.
    :param ls_ph:
    :param E0:
    :param T:
    :return:
    """
    faq = s_nernst * R * (T + Tconvert) / (n * F)
    return 1000 * (E0 - faq * ls_ph)


def Nernst_equation_invert(E, E0, T=25):
    return (E0 - E / 1000) * n * F / (2.303 * R * (T + 273.15))


def calibSensor_para(target_para, sensor2):
    # sensor 2 - targeted concentration fluctuation in mV
    para = _calibration_electroSens(E=sensor2['E'], conc=sensor2['conc_calib'])
    
    df_sigSens2_mV = pd.DataFrame(target_para.filter(like='Sum', axis=1) * para['slope'] + para['E0'])
    df_sigSens2_mV.columns = ['potential mV Sum']

    # what are the specific signal levels
    cplateau = [t * para['slope'] + para['E0'] for t in sensor2['set values']]
    return df_sigSens2_mV, cplateau, para


def _calibration_electroSens(E, conc):
    # Nernst equation with slope -59mV: E = E0 - 59mV/z * log10(c)
    E0 = E + (59/n * np.log(conc))
    para = dict({'slope': 59, 'E0': E0})
    return para


def calibSensor_pH(target_ph, sensor_ph):
    # pH sensor - targeted pH fluctuation in mV
    df_sigpH_mV = pd.DataFrame(Nernst_equation(ls_ph=target_ph['target pH'].to_numpy(), T=temp_std,
                                               E0=sensor_ph['E0']), index=target_ph.index, columns=['potential mV'])
    # what are the specific signal levels
    cplateau = [Nernst_equation(ls_ph=p, T=temp_std, E0=sensor_ph['E0']) for p in sensor_ph['set values']]
    return df_sigpH_mV, cplateau


def para2Sensorcalc(df_res, analyte, cplateauSum, df_pHcalc, sensor2, para_meas, func_response, para_response):
    if analyte == 'base' or analyte == 'Base':
        cBase = [df_res[df_res['Potential mV Sum'] == c].filter(like='Base', axis=1).to_numpy() for c in cplateauSum]
        ls_cBase = list()
        for li in cBase:
            ls_cBase = np.append(ls_cBase, li)
        # find all plateaus even though they repeat themselves
        ls_cBase = _all_plateaus_inlist(ls_cBase)[:len(cplateauSum)]

        df_base_calc = _sensor_response(func=func_response, cplateau=ls_cBase, para_response=para_response,
                                        tplateau=para_meas['plateau time'], cstart=ls_cBase[0])
        df_base_calc.columns = ['Base calc']

        # calculate acid species from Henderson-Hasselbalch and sensed pH
        df_acid_calc = pd.DataFrame(henderson_acid(pH=df_pHcalc['pH calc'].to_numpy(), pKa=sensor2['pKa'],
                                                   c_base=df_base_calc['Base calc'].to_numpy()), columns=['Acid calc'],
                                    index=[round(i, 2) for i in df_pHcalc.index])
    elif analyte == 'acid' or analyte == 'Acid':
        cAcid = [df_res[df_res['Potential mV Sum'] == c].filter(like='Acid', axis=1).to_numpy() for c in cplateauSum]
        ls_cAcid = list()
        for li in cAcid:
            ls_cAcid = np.append(ls_cAcid, li)
        # find all plateaus even though they repeat themselves
        ls_cAcid = _all_plateaus_inlist(ls_cAcid)[:len(cplateauSum)]

        df_acid_calc = _sensor_response(func=func_response, cplateau=ls_cAcid, para_response=para_response, 
                                        tplateau=para_meas['plateau time'], cstart=ls_cAcid[0])

        df_acid_calc.columns = ['Acid calc']

        # calculate base species from Henderson-Hasselbalch and sensed pH
        df_base_calc = pd.DataFrame(henderson_base(pH=df_pHcalc['pH calc'].to_numpy(), pKa=sensor2['pKa'],
                                                   c_acid=df_acid_calc['Acid calc'].to_numpy()), columns=['Base calc'],
                                    index=[round(i, 2) for i in df_pHcalc.index])

    df_sum_calc = pd.DataFrame(df_acid_calc['Acid calc'].to_numpy() + df_base_calc['Base calc'].to_numpy(),
                               index=[round(i, 2) for i in df_pHcalc.index], columns=['Sum calc'])

    # include sensor calculations to result dataframe
    xnew = [round(i, 2) for i in df_sum_calc.index]
    df_acid_calc.index, df_base_calc.index, df_sum_calc.index = xnew, xnew, xnew
    df_res = pd.concat([df_res, df_acid_calc, df_base_calc, df_sum_calc], axis=1).sort_index(axis=1).dropna()

    return df_res


# ------------------------------------------------------------------------------------------------------
# preparing target sensor response and actual sensor response 
def _target_fluctuation(ls_conc, tstart, tstop, analyte, nP=1, steps=0.05):
    """
    Function for periodic block-change of concentration
    :param ls_conc:   list of concentrations (given as floats)
    :param tstart:    start time of the period (in s)
    :param tstop:     stop time of the period (in s)
    :param nP:        frequency; number of cycles/period
    :param analyte:   define the analyte of interest; given as string
    :param steps:     for xdata resolution; given as float
    :return:
    """
    N = ((tstop - tstart) / steps) + 1
    P = N / nP
    D = P / 2

    # when ls_conc has more than 1 entry -> step function
    ls_sig = list(map(lambda n: list((np.arange(N) % P < D) * (ls_conc[n] - ls_conc[n + 1]) + ls_conc[n + 1]),
                      range(len(ls_conc) - 1)))

    x0 = np.arange(0, tstop + steps, steps)
    dsig = dict()
    for i in range(len(ls_sig)):
        if i == 0:
            x = x0
        else:
            x = np.arange((x0[-1]/2+steps)*i, x0[-1]/2*i + steps*(i+1) + x0[-1], steps)
        # additional check to prevent crashing the software
        if len(ls_sig[i]) != len(x):
            x = x[:len(ls_sig[i])]
        dsig[i] = pd.DataFrame(ls_sig[i], columns=['signal ' + analyte], index=x)

    df_target = pd.concat(dsig, axis=0)
    trange = [i[1] for i in df_target.index]
    df_target.index = trange

    # remove duplicates from dataframe
    xnew = [round(x, 2) for x in df_target.index]
    df_target.index = xnew
    df_target = df_target.groupby(df_target.index).mean()

    # finalize output
    df_target.index.name = 'time/s'

    return df_target


def target_fluctuation(ls_conc, tstart, tstop, analyte, nP=1, steps=0.05):
    """
    Function for periodic block-change of concentration
    :param ls_conc:   list of concentrations (given as floats)
    :param tstart:    start time of the period (in s)
    :param tstop:     stop time of the period (in s)
    :param nP:        frequency; number of cycles/period
    :param analyte:   define the analyte of interest; given as string
    :param steps:     for xdata resolution; given as float
    :return:
    """
    if len(ls_conc) > 1:
        df_target = _target_fluctuation(ls_conc=ls_conc, tstart=tstart, tstop=tstop, analyte=analyte, nP=nP,
                                        steps=steps)
    else:
        x0 = np.arange(0, tstop/2 + steps, steps)
        df_target = pd.DataFrame(ls_conc*len(x0), columns=['signal ' + analyte], index=x0)
    return df_target


def pH_sensor(func, cplateau, sensor_ph, para_meas, para_response):
    # include sensor response
    cstart = Nernst_equation(ls_ph=sensor_ph['start pH'], E0=sensor_ph['E0'], T=temp_std)
    df_pHrec = _sensor_response(func=func, cplateau=cplateau, para_response=para_response, tplateau=para_meas['plateau time'], cstart=cstart)
    
    # -------------------------------------
    # re-calculate pH from potential
    df_recalc = pd.DataFrame(Nernst_equation_invert(E=df_pHrec['potential mV'], E0=sensor_ph['E0'], T=temp_std))
    df_recalc.columns = ['pH calc']
    df_recalc.index = [round(i, 2) for i in df_recalc.index]

    return df_pHrec, df_recalc


def calc_integral(df_res, _integral_counter, ls_xcoords, sensor_ph, sensor_para2):
    _integral_counter += 1
    result = dict()
    try:
        ls_xcoords
    except:
        ls_xcoords = list()

    # calculate integral for calculated/target total parameter using simpson's rule for discrete data
    # subtract calculated Sum by target Sum
    if len(ls_xcoords) == 2:
        df_sum = pd.concat([df_res.filter(like='Sum calc'), df_res.filter(like='Sum').filter(like='target')], axis=1)
        df4int = df_sum.loc[min(ls_xcoords):max(ls_xcoords)]
        dfint_calc = integrate.simps(df4int['Sum calc'].to_numpy(), x=df4int.index)

        # target function assumed to have an infinite fast response (~step function)
        search_key = df_sum.filter(like='target').columns[0]
        grp = df4int.groupby(search_key)
        dfint_target = np.sum([(k * (grp.groups[k][-1] - grp.groups[k][0])) for k in grp.groups.keys()])

        # Integration range (s), target pH, t90 pH, target Sum, t90 meas, integral target Sum, integral Sum, error
        result = list([(round(min(ls_xcoords), 2), round(max(ls_xcoords), 2)),
                        df_res['target pH'].loc[min(ls_xcoords):max(ls_xcoords)].mean(),
                        sensor_ph['t90'], df_res[search_key].mean(), sensor_para2['t90'], 
                        dfint_calc, dfint_target, dfint_calc - dfint_target])
    return result


def calculate_integralError(ls_xcoords, df_res, sensor_ph, sensor_para2, para2_name, paraSum_conc):
    if isinstance(ls_xcoords, list):
        # ls_xcoords, average target pH, t90 pH,  average target Sum, t90 sensor,  calculated integral, target integral,
        # difference between target and calculated integral == error
        dres_all, dic_out = dict(), dict()
        for en, xcoords in enumerate(ls_xcoords):
            # determine the actuall error compared to the target integral. The calculation is based on simpson's integral for sample data
            dres = calc_integral(df_res=df_res, _integral_counter=0, ls_xcoords=xcoords, sensor_ph=sensor_ph, 
                                sensor_para2=sensor_para2)
            dres_all[en] = dres

            # combine calculations in a joint DataFrame for better overview
            if isinstance(dres, list):
                out = pd.DataFrame([str(dres[0])]+dres[1:], index=['integral range', 'pH_target', 't90_pH', 'SUM_target, {}'.format(paraSum_conc[1]), 
                                                                   't90_{}'.format(para2_name), 'integral_calc', 'integral_target', 'error'])
            dic_out[en] = out

        # if multiple integral ranges shall be calculated, combine them (again) in a joint DataFrame 
        df_out_all = pd.concat(dic_out, axis=1)
        df_out_all.columns = df_out_all.T.index.levels[0]
        return df_out_all
    elif isinstance(ls_xcoords, tuple):
        # determine the actuall error compared to the target integral. The calculation is based on simpson's integral for sample data
        dres = calc_integral(df_res=df_res, _integral_counter=0, ls_xcoords=ls_xcoords, sensor_ph=sensor_ph, sensor_para2=sensor_para2)

        # combine calculations in a joint DataFrame for better overview
        if isinstance(dres, list):
            out = pd.DataFrame([str(dres[0])]+dres[1:], index=['integral range', 'pH_target', 't90_pH', 'SUM_target, {}'.format(paraSum_conc[1]), 't90_{}'.format(para2_name),
                                                                'integral_calc', 'integral_target', 'error'])
        return out
    else:
        print('Warning! Nothing calculated - please provide either a list of tuples or one single tuple defining the integration ranges.') 
        return None
    

# ----------------------------------------------------------------------------
# main function
def run_simulation(sensor_ph, sensor_para2, para_meas, paraSum_unit):
    global ls_lines, _integral_counter
    ls_lines, ls_xcoords, vlines_list = [], [], []

    # target concentrations
    df_res = pd.concat([sensor_ph['pH target'], sensor_para2['Base target'], sensor_para2['Acid target'], 
                        sensor_para2['Sum target']], axis=1)
    ## adjust the unit
    df_res.columns = ['target pH', 'target_{} Base'.format(paraSum_unit), 'target_{} Acid'.format(paraSum_unit), 
                      'target_{} Sum'.format(paraSum_unit)]
    
    # individual sensor - reduce by individual plotting (different function)
    [df_res, cplateaupH, cplateauTotal, 
     para_para] = _alignSensorSettings(sensor_ph=sensor_ph, sensor2=sensor_para2, df_target=df_res)
    df_res.index = [round(u, 2) for u in df_res.index]
    df_res = df_res.dropna()

    # pH sensor response
    [df_pHrec, df_pHcalc] = pH_sensor(cplateau=cplateaupH, sensor_ph=sensor_ph, para_meas=para_meas)
    df_pHrec.index = [round(i, 2) for i in df_pHrec.index]

    # include sensor response and drift to result dataframe
    df_res = pd.concat([df_res, df_pHrec.loc[df_res.index], df_pHcalc.loc[df_res.index]], axis=1).sort_index(axis=1)

    # para2 sensor - calculate individual sensor response as well as total parameter
    df_res = para2Sensorcalc(df_res=df_res, analyte=sensor_para2['analyte'], cplateauSum=cplateauTotal, 
                             df_pHcalc=df_pHcalc, sensor2=sensor_para2, para_meas=para_meas)
    df_res = df_res.sort_index(axis=1)

    return df_res, ls_lines, ls_xcoords, vlines_list


def run_simulation_divResponses(sensor_ph, sensor_para2, para_meas, paraSum_unit, func_resp_pH, presponse_ph, func_resp_meas, presponse_meas):
    global ls_lines, _integral_counter, df_res
    ls_lines, ls_xcoords, vlines_list = [], [], []
        
    # target concentrations
    df_res = pd.concat([sensor_ph['pH target'], sensor_para2['Base target'], sensor_para2['Acid target'], 
                        sensor_para2['Sum target']], axis=1)
    ## adjust the unit
    df_res.columns = ['target pH', 'target_{} Base'.format(paraSum_unit), 'target_{} Acid'.format(paraSum_unit), 
                        'target_{} Sum'.format(paraSum_unit)]

    # individual sensor - reduce by individual plotting (different function)
    [df_res, cplateaupH, cplateauTotal, 
    para_para] = _alignSensorSettings(sensor_ph=sensor_ph, sensor2=sensor_para2, df_target=df_res)
    df_res.index = [round(u, 2) for u in df_res.index]
    df_res = df_res.dropna()

    # ----------------------------------------------------
    # pH sensor response
    [df_pHrec, df_pHcalc] = pH_sensor(func=func_resp_pH, cplateau=cplateaupH, sensor_ph=sensor_ph, para_meas=para_meas, 
                                      para_response=presponse_ph)
    df_pHrec.index = [round(i, 2) for i in df_pHrec.index]

    # include sensor response and drift to result dataframe
    df_res = pd.concat([df_res, df_pHrec.loc[df_res.index], df_pHcalc.loc[df_res.index]], axis=1).sort_index(axis=1)

    # para2 sensor - calculate individual sensor response as well as total parameter
    df_res = para2Sensorcalc(func_response=func_resp_meas, df_res=df_res, analyte=sensor_para2['analyte'], cplateauSum=cplateauTotal, 
                            df_pHcalc=df_pHcalc, sensor2=sensor_para2, para_meas=para_meas, para_response=presponse_meas)
    df_res = df_res.sort_index(axis=1)

    return df_res, ls_lines, ls_xcoords, vlines_list

