__author__ = 'Silvia E Zieger'
__project__ = 'sensor response in silico'

"""Copyright 2021. All rights reserved.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable 
for any damages arising from the use of this software.
Permission is granted to anyone to use this software within the scope of evaluating mutli-analyte sensing. No permission
is granted to use the software for commercial applications, and alter it or redistribute it.

This notice may not be removed or altered from any distribution.
"""

import matplotlib
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import numpy as np
import pandas as pd
import scipy as sp

# --------------------------------------------------------------------------------------------------------------------
# global parameter
n = 1
R = 8.314         # [J/mol*K]
F = 96485         # [C/mol]
s_nernst = 2.303
Tconvert = 273.15

dcolor = dict({'pH': '#1CC49A', 'sig pH': '#4B5258', 'NH3': '#196E94', 'sig NH3': '#314945',
               'NH4': '#DCA744', 'sig NH4': '#89621A', 'TAN': '#A86349'})
ls = dict({'target': '-.', 'simulation': ':'})
fs = 11


# --------------------------------------------------------------------------------------------------------------------
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


def lin_regression(cnh3_min, cnh3_max, sigNH3_bgd, sigNH3_max, conc_step=0.1):
    cnh3_cali = np.linspace(cnh3_min, cnh3_max, num=int((cnh3_max - cnh3_min) / conc_step + 1))
    ynh3_cali = np.linspace(sigNH3_bgd, sigNH3_max, num=len(cnh3_cali))

    arg = sp.stats.linregress(x=cnh3_cali, y=ynh3_cali)
    return arg


def _calibration_nhx(Emax, cmax):
    # version1: anhx_min, anhx_max, anhx_step, sigNHx_bgd, sigNHx_max):
    # linear regression for all pH values
    #nhx_cali = lin_regression(cnh3_min=anhx_min, cnh3_max=anhx_max, conc_step=anhx_step, sigNH3_bgd=sigNHx_bgd,
    #                          sigNH3_max=sigNHx_max)
    # Nernst equation with slope -59mV: E = E0 - 59mV/z * log10(c)
    E0 = Emax + (59/n * np.log(cmax))
    para = dict({'slope': 59, 'E0': E0})
    return para


def henderson_nh3(pH, c_nh4, pKa=9.25):
    """ Returns the NH3 concentration depending on the pH and the concentrations of NH4+ """
    return c_nh4 * 10 ** (pH - pKa)


def henderson_nh4(pH, c_nh3, pKa=9.25):
    """ Returns the NH3 concentration depending on the pH and the concentrations of NH4+ """
    return c_nh3 * 10 ** (pKa - pH)


def henderson_TAN4NH3(pH, c_tan, pKa=9.25):
    """
    Returns the NH3 concentration depending on the pH and the TAN concentration
    :param pH:
    :param c_tan:
    :param pKa:
    :return:
    """
    return c_tan / (1 + 10**(pKa - pH))


def henderson_TAN4NH4(pH, c_tan, pKa=9.25):
    """
    Returns the NH3 concentration depending on the pH and the TAN concentration
    :param pH:
    :param c_tan:
    :param pKa:
    :return:
    """
    return c_tan / (1 + 10**(pH - pKa))


# --------------------------------------------------------------------------------------------------------------------
def _target_fluctuation(ls_conc, tstart, tstop, analyte, nP=1, steps=0.05):
    """
    Function for periodic block-change of concentration
    :param cnh4_low:  minimal ammonium concentration
    :param cnh4_high: maximal ammonium concentration
    :param tstart:    start time of the period (in s)
    :param tstop:     stop time of the period (in s)
    :param nP:        frequency; number of cycles/period
    :param N:         sample count; a multitude of the time period given
    :param D:         width of pulse; usually half of the period P
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
    :param cnh4_low:  minimal ammonium concentration
    :param cnh4_high: maximal ammonium concentration
    :param tstart:    start time of the period (in s)
    :param tstop:     stop time of the period (in s)
    :param nP:        frequency; number of cycles/period
    :param N:         sample count; a multitude of the time period given
    :param D:         width of pulse; usually half of the period P
    :return:
    """

    if len(ls_conc) > 1:
        df_target = _target_fluctuation(ls_conc=ls_conc, tstart=tstart, tstop=tstop, analyte=analyte, nP=nP,
                                        steps=steps)
    else:
        x0 = np.arange(0, tstop/2 + steps, steps)
        df_target = pd.DataFrame(ls_conc*len(x0), columns=['signal ' + analyte], index=x0)
    return df_target


def _target_concentration(ls_ph, ls_cNH3_ppm, ls_cNH4_ppm, t_plateau):
    # always pH
    [target_ph, trange, D] = _target_fluctuation(ls_conc=ls_ph, tstart=0, tstop=t_plateau * 2, nP=1)

    # define analyte - whether NH3 or NH4+
    if np.all([n == None for n in ls_cNH3_ppm]) and not np.all([n == None for n in ls_cNH4_ppm]):
        analyte, ls_nhx = 'NH4', ls_cNH4_ppm
    elif not np.all([n == None for n in ls_cNH3_ppm]) and np.all([n == None for n in ls_cNH4_ppm]):
        analyte, ls_nhx = 'NH3', ls_cNH3_ppm
    else:
        print('ERROR - define a concentration of NH3 / NH4 you want to study')
        print('predefined NH3: 100ppm')
        ls_nhx = 100
        analyte = 'NH3'

    # specify NHx concentrations - consider all options
    ls_nh = tuple()
    for en, n in enumerate(ls_nhx):
        if n != None:
            tn = (n,)
            ls_nh += tn
    if len(ls_nh) == 1:
        ls_nh = ls_nh[0]
    [target_nhx, trange, D] = _target_fluctuation(ls_conc=ls_nh, tstart=0, tstop=trange[-1], nP=1)

    return D, target_ph, target_nhx, analyte


def signal_drift(xtime, df_mV, psensor):
    df_drift = pd.DataFrame(xtime * psensor['drift'] + df_mV['potential mV'].to_numpy())
    df_drift.columns, df_drift.index = ['potential mV'], xtime
    return df_drift


def _sensor_response(cplateau, psensor, tplateau):
    start, steps = 0, 0.05
    num = int((tplateau - start) / 0.05 + 1)
    sens_time = np.linspace(start, tplateau, num=num)

    # simulate for each signal level a respective sensor response
    c_apparent = psensor['background signal']
    dsig = dict()
    for en, ctarget in enumerate(cplateau):
        sdiff = ctarget - c_apparent
        if en == 0:
            xnew = sens_time
        else:
            xnew = np.linspace((tplateau + steps) * en, (tplateau + steps) * en + tplateau, num=num)
        df_sig = pd.DataFrame(_gompertz_curve_v1(x=sens_time, t90=psensor['t90'], s_diff=sdiff, pstart=c_apparent,
                                                 tau=psensor['resolution'], slope='increase'), index=xnew,
                              columns=['potential mV'])

        start, c_apparent = df_sig['potential mV'].index[-1], df_sig['potential mV'].to_numpy()[-1]
        dsig[en] = df_sig
        df_sensor = pd.concat(dsig)
        df_sensor.index = [i[1] for i in df_sensor.index]

    # include sensor drift if applicable
    df_drift = signal_drift(xtime=df_sensor.index, df_mV=df_sensor, psensor=psensor)

    return df_sensor, df_drift


# --------------------------------------------------------------------------------------------------------------------
def calibSensor_pH(target_ph, sensor_ph, para_meas):
    # pH sensor - targeted pH fluctuation in mV
    df_sigpH_mV = pd.DataFrame(Nernst_equation(ls_ph=target_ph['target pH'].to_numpy(), T=para_meas['temperature'],
                                               E0=sensor_ph['E0']), index=target_ph.index, columns=['potential mV'])
    # what are the specific signal levels
    cplateau = [Nernst_equation(ls_ph=p, T=para_meas['temperature'], E0=sensor_ph['E0'])
                for p in sensor_ph['set values']]
    return df_sigpH_mV, cplateau


def calibSensor_NHx(target_nhx, sensor_nh3, analyte):
    # NHx sensor - targeted concentration fluctuation in mV
    para_nhx = _calibration_nhx(Emax=sensor_nh3['signal max ' + analyte], cmax=sensor_nh3['calib ' + analyte])
    df_sigNHx_mV = pd.DataFrame(target_nhx['target_ppm ' + analyte] * para_nhx['slope'] + para_nhx['E0'])
    df_sigNHx_mV.columns = ['potential mV ' + analyte]

    # what are the specific signal levels
    cplateau = [t * para_nhx['slope'] + para_nhx['E0'] for t in sensor_nh3['set values']]
    return df_sigNHx_mV, cplateau, para_nhx


def _alignSensorSettings(df_target, sensor_ph, sensor_nh3, para_meas):
    # individual sensor calibration
    df_sigpH_mV, cplateaupH = calibSensor_pH(target_ph=df_target, sensor_ph=sensor_ph, para_meas=para_meas)
    df_sigNHx_mV, cplateauNHx, para_nhx = calibSensor_NHx(target_nhx=df_target, sensor_nh3=sensor_nh3,
                                                          analyte=sensor_nh3['analyte'])
    df_ = pd.concat([df_sigpH_mV, df_sigNHx_mV], axis=1)
    df_.columns = ['Potential mV pH', 'Potential mV ' + sensor_nh3['analyte']]
    df_res = pd.concat([df_target, df_], axis=1)
    return df_res, cplateaupH, cplateauNHx, para_nhx


def pH_sensor(cplateau, sensor_ph, para_meas):
    # include sensor response
    df_pHrec, df_pHdrift = _sensor_response(cplateau=cplateau, psensor=sensor_ph, tplateau=para_meas['plateau time'])
    df_pHdrift.columns = ['Drift pH mV']

    # -------------------------------------
    # re-calculate pH from potential
    df_recalc = pd.DataFrame(Nernst_equation_invert(E=df_pHdrift['Drift pH mV'], E0=sensor_ph['E0'],
                                                    T=para_meas['temperature']))
    df_recalc.columns = ['pH calc']
    ind_new = [round(i, 2) for i in df_recalc.index]
    df_recalc.index = ind_new

    return df_pHrec, df_pHdrift, df_recalc


def NHx_sensor(cplateau, sensor_nh3, para_nhx, para_meas):
    # update concentration in case one analyte changes more often than the other one
    # include sensor response
    df_nhrec, df_nhdrift = _sensor_response(cplateau=cplateau, psensor=sensor_nh3, tplateau=para_meas['plateau time'])
    df_nhdrift.columns = ['Drift mV ' + sensor_nh3['analyte']]

    # -------------------------------------
    # re-calculate NHx from potential
    df_recalc = pd.DataFrame((df_nhdrift - para_nhx['E0']) / para_nhx['slope'])
    df_recalc.columns = [sensor_nh3['analyte'] + ' calc']
    ind_new = [round(i, 2) for i in df_recalc.index]
    df_recalc.index = ind_new

    return df_nhrec, df_nhdrift, df_recalc


def _other_analyte(analyte, sensor_nh3, df_target, df_calc):
    if analyte == 'NH3':
        # based on target parameter
        c_nhx = henderson_nh4(pKa=sensor_nh3['pKa'], pH=df_target['target pH'].to_numpy(),
                              c_nh3=df_target['target_ppm ' + analyte].to_numpy())
        df_target.loc[:, 'NH4'] = c_nhx

        # based on sensor response
        c_nhx = henderson_nh4(pH=df_calc['pH calc'].to_numpy(), c_nh3=df_calc[analyte + ' calc'].to_numpy(),
                              pKa=sensor_nh3['pKa'])

        df_calc.loc[:, 'NH4 calc'] = c_nhx
        df_calc.loc[:, 'TAN calc'] = df_calc['NH4 calc'].to_numpy() + df_calc['NH3 calc'].to_numpy()
    elif analyte == 'NH4':
        # based on target parameter
        c_nhx = henderson_nh3(pKa=sensor_nh3['pKa'], pH=df_target['target pH'].to_numpy(),
                              c_nh4=df_target['target_ppm ' + analyte].to_numpy())
        df_target.loc[:, 'NH3'] = c_nhx

        # based on sensor response
        c_nhx = henderson_nh3(pH=df_calc['pH calc'].to_numpy(), c_nh4=df_calc[analyte + ' calc'].to_numpy(),
                              pKa=sensor_nh3['pKa'])
        df_calc.loc[:, 'NH3 calc'] = c_nhx
        df_calc.loc[:, 'TAN calc'] = df_calc['NH4 calc'].to_numpy() + df_calc['NH3 calc'].to_numpy()
    else:
        df_target, df_calc = None, None
    return df_target, df_calc


def _tan_calculation(df_target, df_calc):
    # targeted
    df_target['TAN'] = np.sum(df_target.filter(like='NH'), axis=1)
    df_target = df_target.dropna()

    # based on sensor response
    df_calc['TAN'] = np.sum(df_calc.filter(like='NH'), axis=1)
    return df_target, df_calc


def individualAnalytes(analyte, sensor_nh3, df):
    # prep data for same time range
    df_prep = df.dropna()

    # calculate individual analytes based on pH and TAN
    if analyte == 'NH3':
        # NH3 was measured
        df_nh3 = pd.DataFrame(henderson_TAN4NH3(pH=df_prep['pH'].to_numpy(), c_tan=df_prep['NH4 calc'].to_numpy(),
                                                pKa=sensor_nh3['pKa']), index=df_prep.index, columns=['NH3'])
        # NH4 as the difference
        df_nh4 = pd.DataFrame(df_prep['NH4 calc'] - df_nh3['NH3'], columns=['NH4'])
    elif analyte == 'NH4':
        # NH4 was measured
        df_nh4 = pd.DataFrame(henderson_TAN4NH4(pH=df_prep['pH'].to_numpy(), c_tan=df_prep['TAN'].to_numpy(),
                                                pKa=sensor_nh3['pKa']), index=df_prep.index, columns=['NH4'])
        # NH3 as the difference
        df_nh3 = pd.DataFrame(df_prep['TAN'] - df_nh4['NH4'], columns=['NH3'])
    else:
        df_nh3, df_nh4 = None, None

    df_calc = pd.concat([df_prep, df_nh3, df_nh4], axis=1)
    return df_calc


# --------------------------------------------------------------------------------------------------------------------
def move_figure(xnew, ynew):
    mngr = plt.get_current_fig_manager()
    geom = mngr.window.geometry()
    x, y, dx, dy = geom.getRect()
    mngr.window.setGeometry(xnew, ynew, dx, dy)


def save_report(para_meas, sensor_ph, sensor_nh3, df_res):
    df_p = pd.DataFrame(np.zeros(shape=(len(para_meas.values()), 2)))
    df_p[0] = list(para_meas.keys())
    df_p[1] = para_meas.values()
    df_p.loc[-1, :] = ['parameter', 'values']
    df_p = df_p.sort_index()
    df_p.columns = ['parameter', 'values']
    df_p.index = ['general'] * len(df_p.index)

    df_ph = pd.DataFrame(np.zeros(shape=(len(sensor_ph.values()), 2)))
    df_ph[0] = list(sensor_ph.keys())
    df_ph[1] = sensor_ph.values()
    ph_target = list(dict.fromkeys(sensor_ph['pH target'].to_list()))
    df_ph.loc[7] = ['pH target', ph_target]
    df_ph.columns = ['parameter', 'values']
    df_ph.index = ['ph'] * len(df_ph.index)

    df_nh3 = pd.DataFrame(np.zeros(shape=(len(sensor_nh3.values()), 2)))
    df_nh3[0] = list(sensor_nh3.keys())
    df_nh3[1] = sensor_nh3.values()
    nhx_target = list(dict.fromkeys(sensor_nh3['{} target'.format(sensor_nh3['analyte'])].to_list()))
    df_nh3.loc[6] = ['{} target'.format(sensor_nh3['analyte']), nhx_target]
    df_nh3.columns = ['parameter', 'values']
    df_nh3.index = [sensor_nh3['analyte']] * len(df_nh3.index)
    df_para = pd.concat([df_p, df_ph, df_nh3])

    # ..................................................................
    # results
    df_calc = df_res.filter(like='calc') # pd.concat([df_, df_res_]).sort_index().T.sort_index().T
    xnew = [int(i) for i in df_calc.index]
    df_calc.index = xnew
    df_calc = df_calc.groupby(df_calc.index).mean()
    header_res = pd.DataFrame(df_calc.columns, columns=['Time / s'], index=df_calc.columns).T

    df_out = pd.concat([header_res, df_calc])
    df_para.columns = [0, 1]
    df_out.columns = np.arange(0, len(df_out.columns))

    output = pd.concat([df_para, df_out])

    return output


def load_data(file):
    file_ = open(file, 'r')
    count = 0
    ls_lines = list()
    while True:
        count += 1
        line = file_.readline()
        # if line is empty end of file is reached
        if not line:
            break
        ls_lines.append(line.strip().split('\t'))
    file_.close()

    # ............................................................
    ls_general = list()
    for l in ls_lines:
        if l[0] == 'general':
            ls_general.append(l[1:])
    df_general = pd.DataFrame(ls_general).T.set_index(0).T.set_index('parameter')

    ls_ph = list()
    for l in ls_lines:
        if l[0] == 'ph':
            ls_ph.append(l[1:])
    df_ph = pd.DataFrame(ls_ph, columns=['parameter', 'values'])
    df_ph = df_ph.set_index('parameter')

    ls_nh3 = list()
    for l in ls_lines:
        if 'NH' in l[0]:
            ls_nh3.append(l[1:])
    df_nh3 = pd.DataFrame(ls_nh3, columns=['parameter', 'values'])
    df_nh3 = df_nh3.set_index('parameter')

    return df_general, df_ph, df_nh3
