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


def _calibration_nh3(anh3_min, anh3_max, anh3_step, sigNH3_bgd, sigNH3_max):
    # linear regression for all pH values
    nh3_cali = lin_regression(cnh3_min=anh3_min, cnh3_max=anh3_max, conc_step=anh3_step, sigNH3_bgd=sigNH3_bgd,
                              sigNH3_max=sigNH3_max)
    return nh3_cali


def henderson_nh3(pH, c_nh4, pKa=9.25):
    """ Returns the NH3 concentration depending on the pH and the concentrations of NH4+ """
    return c_nh4 * 10 ** (pH - pKa)


def henderson_nh4(pH, c_nh3, pKa=9.25):
    """ Returns the NH3 concentration depending on the pH and the concentrations of NH4+ """
    return c_nh3 * 10 ** (pKa - pH)


# --------------------------------------------------------------------------------------------------------------------
def _target_fluctuation(ls_conc, tstart, tstop, nP=1):
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

    N = ((tstop - tstart) / 0.05) + 1
    P = N / nP
    D = P / 2

    # when ls_conc has more than 1 entry -> step function
    if isinstance(ls_conc, tuple):
        ls_sig = list(map(lambda n: list((np.arange(N) % P < D) * (ls_conc[n] - ls_conc[n + 1]) + ls_conc[n + 1]),
                          range(len(ls_conc) - 1)))

        dsig = dict()
        for i in range(len(ls_sig)):
            x = np.linspace(tstop / 2 * i, tstop + tstop / 2 * i, num=int(N))
            dsig[i] = pd.DataFrame(ls_sig[i], columns=['signal'], index=x)

        df_target = pd.concat(dsig, axis=0)
        trange = [i[1] for i in df_target.index]
        df_target.index = trange
    else:
        trange = np.linspace(0, tstop, num=int(N))
        df_target = pd.DataFrame([ls_conc] * len(trange), index=trange, columns=['signal'])

    df_target.index.name = 'time/s'

    return df_target, trange, D


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


def _sensor_response(df_mV, psensor):
    # what are the specific signal levels
    conc_plateau = list(dict.fromkeys(df_mV['potential mV'].to_numpy()))

    # timing - define sensor time (for simulation of response)
    tplateau = df_mV[df_mV == conc_plateau[0]].dropna().index[-1]
    start, steps = 0, 0.05
    num = int((tplateau - start) / 0.05 + 1)
    sens_time = np.linspace(start, tplateau, num=num)

    # simulate for each signal level a respective sensor response
    c_apparent = psensor['background signal']
    dsig = dict()
    for en, ctarget in enumerate(conc_plateau):
        sdiff = ctarget - c_apparent
        xnew = np.linspace(start + steps * en, tplateau * (en + 1) + steps * en, num=num)

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
def pH_sensor(target_ph, sensor_ph, para_meas):
    # pH sensor - targeted pH fluctuation in mV
    df_sigpH_mV = pd.DataFrame(Nernst_equation(ls_ph=target_ph['signal'].to_numpy(), T=para_meas['temperature'],
                                               E0=sensor_ph['E0']), index=target_ph.index, columns=['potential mV'])

    # -------------------------------------
    # include sensor response
    df_pHrec, df_pHdrift = _sensor_response(df_mV=df_sigpH_mV, psensor=sensor_ph)

    # -------------------------------------
    # re-calculate pH from potential
    df_recalc = pd.DataFrame(Nernst_equation_invert(E=df_pHdrift['potential mV'], E0=sensor_ph['E0'],
                                                    T=para_meas['temperature']))
    df_recalc.columns = ['pH']
    ind_new = [round(i, 2) for i in df_recalc.index]
    df_recalc.index = ind_new

    return df_sigpH_mV, df_pHrec, df_pHdrift, df_recalc


def NHx_sensor(analyte, target_nhx, sensor_nh3):
    # NHx sensor - targeted concentration fluctuation in mV
    para_nh3 = _calibration_nh3(anh3_min=sensor_nh3['NHx calibration'][0], anh3_max=sensor_nh3['NHx calibration'][1],
                                anh3_step=5, sigNH3_bgd=sensor_nh3['signal min'], sigNH3_max=sensor_nh3['signal max'])
    df_sig_mV = pd.DataFrame(target_nhx * para_nh3[0] + para_nh3[1])
    df_sig_mV.columns = ['potential mV']

    # -------------------------------------
    # include sensor response
    df_nhrec, df_nhdrift = _sensor_response(df_mV=df_sig_mV, psensor=sensor_nh3)

    # -------------------------------------
    # re-calculate NHx from potential
    df_recalc = pd.DataFrame((df_nhdrift - para_nh3[1]) / para_nh3[0])
    df_recalc.columns = [analyte]
    ind_new = [round(i, 2) for i in df_recalc.index]
    df_recalc.index = ind_new

    return df_sig_mV, df_nhrec, df_nhdrift, df_recalc


def _other_analyte(analyte, sensor_nh3, df_target, df_calc):
    if analyte == 'NH3':
        # based on target parameter
        c_nhx = henderson_nh4(pKa=sensor_nh3['pKa'], pH=df_target['pH'].to_numpy(), c_nh3=df_target[analyte].to_numpy())
        df_target['NH4'] = pd.DataFrame(c_nhx, index=df_target.index)
        # based on sensor response
        c_nhx = henderson_nh4(pH=df_calc['pH'].to_numpy(), c_nh3=df_calc[analyte].to_numpy(), pKa=sensor_nh3['pKa'])
        df_calc['NH4'] = pd.DataFrame(c_nhx, index=df_calc.index)

    else:
        # based on target parameter
        c_nhx = henderson_nh3(pKa=sensor_nh3['pKa'], pH=df_target['pH'].to_numpy(), c_nh4=df_target[analyte].to_numpy())
        df_target['NH3'] = pd.DataFrame(c_nhx, index=df_target.index)

        # based on sensor response
        c_nhx = henderson_nh3(pH=df_calc['pH'].to_numpy(), c_nh4=df_calc[analyte].to_numpy(), pKa=sensor_nh3['pKa'])
        df_calc['NH3'] = pd.DataFrame(c_nhx, index=df_calc.index)
    return df_target, df_calc


def _tan_calculation(df_target, df_calc):
    # targeted
    df_target['TAN'] = np.sum(df_target.filter(like='NH'), axis=1)
    df_target = df_target.dropna()

    # based on sensor response
    df_calc['TAN'] = np.sum(df_calc.filter(like='NH'), axis=1)
    return df_target, df_calc


# --------------------------------------------------------------------------------------------------------------------
def move_figure(xnew, ynew):
    mngr = plt.get_current_fig_manager()
    geom = mngr.window.geometry()
    x, y, dx, dy = geom.getRect()
    mngr.window.setGeometry(xnew, ynew, dx, dy)


def save_report(para_meas, sensor_ph, sensor_nh3, dsens_record, dtarget):
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
    df_ph.columns = ['parameter', 'values']
    df_ph.index = ['ph'] * len(df_ph.index)

    df_nh3 = pd.DataFrame(np.zeros(shape=(len(sensor_nh3.values()), 2)))
    df_nh3[0] = list(sensor_nh3.keys())
    df_nh3[1] = sensor_nh3.values()
    df_nh3.columns = ['parameter', 'values']
    df_nh3.index = ['nh3'] * len(df_nh3.index)

    df_para = pd.concat([df_p, df_ph, df_nh3])
    # ..................................................................
    # results
    df_res_ = pd.DataFrame.from_dict(dsens_record)
    df_ = pd.DataFrame.from_dict(dtarget)
    df_.columns = ['TAN_target', '{}_target / ppm'.format(df_.columns[1].split(' ')[0]), 'pH_target']

    df_res = pd.concat([df_, df_res_]).sort_index().T.sort_index().T
    xnew = [int(i) for i in df_res.index]
    df_res.index = xnew
    df_res = df_res.groupby(df_res.index).mean()
    header_res = pd.DataFrame(df_res.columns, columns=['Time / s'], index=df_res.columns).T

    df_out = pd.concat([header_res, df_res])
    df_para.columns = [0, 1]
    df_out.columns = np.arange(0, len(df_out.columns))

    output = pd.concat([df_para, df_out])

    return output
