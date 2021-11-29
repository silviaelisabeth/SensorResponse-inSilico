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
def _gompertz_curve(trange, t_resp, Sbackground, c_max):
    b = np.log(np.log(c_max / Sbackground))
    k = 1 / t_resp * (b - np.log(np.log(1 / 0.9)))
    y = c_max * np.exp(-1 * np.exp(b - k * trange))
    return y


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


def _sensor_response(t, par_sens, y=0, direction='increase', function='gompertz'):
    """
    Positive S-shaped sensor response along a certain time range. Options to choose from are 1) whether it is an
    increase or a decline in the sensor signal. 2) the function to be used is either an empirical or gompertz.
    :param t:           time range
    :param par_sens:    sensor parameters including maximal concentration, response time, background signal
    :param y:           apparent signal
    :param direction:   increase or decline of the sensor response
    :param function:    empirical or gompertz function
    :return:
    """

    if direction == 'increase':
        if function == 'gompertz' or function == 'Gompertz':
            dfSig = pd.DataFrame(_gompertz_curve(trange=t, c_max=par_sens['max conc'], t_resp=par_sens['response time'],
                                                 Sbackground=par_sens['background signal']), index=t, columns=['signal'])
        elif function == 'empirical':
            k = ((1 / 0.9) ** 2 - 1) * par_sens['response time'] ** 2
            dfSig = pd.DataFrame(t * (par_sens['max conc'] - y) / np.sqrt(k + t ** 2) + y, index=t, columns=['signal'])
        else:
            raise ValueError('choose either gompertz or empirical as function')

    elif direction == 'decline':
        if function == 'gompertz' or function == 'Gompertz':
            df = pd.DataFrame(_gompertz_curve(trange=t, c_max=par_sens['max conc'], t_resp=par_sens['response time'],
                                              Sbackground=par_sens['background signal']), index=t, columns=['signal'])
            dfSig = y - df
        elif function == 'empirical':
            k = ((1 / 0.9) ** 2 - 1) * par_sens['response time'] ** 2
            dfSig = pd.DataFrame((-1 * t * par_sens['max conc'] / np.sqrt(k + t ** 2)) + y, index=t, columns=['signal'])
        else:
            raise ValueError('choose either gompertz or empirical as function')

    else:
        raise ValueError('Define function parameter; use either increase or decline.')

    #t10 = dfSig[dfSig['signal'] >= dfSig['signal'].max() * 0.1].index[0]
    #t90 = dfSig[dfSig['signal'] >= dfSig['signal'].max() * 0.9].index[0]
    return dfSig#, t10, t90


def _pHsens_response(t, sig_max, sig_app, sig_bgd, tresp, direction='increase'):
    """
    Positive S-shaped sensor response along a certain time range. Options to choose from are 1) whether it is an
    increase or a decline in the sensor signal. 2) the function to be used is either an empirical or gompertz.
    :param t:           time range
    :param par_sens:    sensor parameters including maximal concentration, response time, background signal
    :param y:           apparent signal
    :param direction:   increase or decline of the sensor response
    :param function:    empirical or gompertz function
    :return:
    """

    if direction == 'increase':
        dfSig = pd.DataFrame(_gompertz_curve(trange=t, c_max=sig_max, t_resp=tresp, Sbackground=sig_bgd), index=t,
                                             columns=['signal'])
    elif direction == 'decline':
        df = pd.DataFrame(_gompertz_curve(trange=t, c_max=sig_max, t_resp=tresp, Sbackground=sig_bgd), index=t,
                                          columns=['signal'])
        dfSig = sig_app - df
    else:
        raise ValueError('Define function parameter; use either increase or decline.')
    return dfSig


def _target_concentration(c_nh3, tstart, tstop, nP=1, N=None, D=None):
    """
    Function for periodic block-change of concentration
    :param c_nh3:   maximal ammonia concentration
    :param tstart:  start time of the period (in s)
    :param tstop:   stop time of the period (in s)
    :param nP:      frequency; number of cycles/period
    :param N:       sample count; a multitude of the time period given
    :param D:       width of pulse; usually half of the period P
    :return:
    """

    if N is None:
        N = tstop - tstart
    P = N / nP
    if D is None:
        D = P / 2

    sig = (np.arange(N) % P < D) * c_nh3
    trange = np.linspace(tstart, tstop, num=int(N))
    dfconc = pd.DataFrame(sig, columns=['signal'], index=trange)
    dfconc.index.name = 'time/s'

    return dfconc, D, P, trange


def _sampling_time(trange, Wpulse, tstart, tend, sampling_rate=None):
    if sampling_rate is None:
        s = 0
        for en, i in enumerate(trange):
            if i >= Wpulse:
                if s == 0:
                    s = en
        srate = trange[s]
    else:
        srate = sampling_rate
    tsampling = np.arange(tstart, tend+srate, srate)
    return tsampling, srate


def henderson_hasselbalch(c_nh3, c_nh4, pKa=9.25):
    """ Returns the pH depending on the concentrations of NH3 and NH4+ """
    return pKa + np.log(c_nh3 / c_nh4)


def henderson_nh3(pH, c_nh4, pKa=9.25):
    """ Returns the NH3 concentration depending on the pH and the concentrations of NH4+ """
    return c_nh4 * 10 ** (pH - pKa)


def henderson_nh4(pH, c_nh3, pKa=9.25):
    """ Returns the NH3 concentration depending on the pH and the concentrations of NH4+ """
    return c_nh3 * 10 ** (pKa - pH)


def total_ammonia(c_nh4, pH, pKa=9.25):
    """ Returns the total ammonia concentration depending on the pH and the
    concentrations of NH4+ """
    if isinstance(pH, list):
        tan = [c_nh4 * (1 + 10 ** (ph_i - pKa)) for ph_i in ls_pH]
    else:
        tan = c_nh4 * (1 + 10 ** (pH - pKa))
    return tan


def partial_concNH3(c_nh4, pH, pKa=9.25):
    """ Returns the percentage of NH3 concentration depending on the pH """
    if isinstance(pH, list):
        numerator = [c_nh4 * 10 ** (ph_i - pKa) for ph_i in ls_pH]
        denominator = [c_nh4 * (1 + 10 ** (ph_i - pKa)) for ph_i in ls_pH]
        c_nh3 = [n / d for (n, d) in zip(numerator, denominator)]
    else:
        c_nh3 = (c_nh4 * (10 ** (pH - pKa))) / (c_nh4 * (1 + 10 ** (pH - pKa)))
    return c_nh3


def int2pH(intF, int_max, pKa):
    # pre-check - intensity out of range
    f_min = round(int_max * np.exp(-pKa) / (1 + np.exp(-pKa)), 5)          # pH >= 0
    f_max = round(int_max * np.exp(14 - pKa) / (1 + np.exp(14 - pKa)), 5)  # pH <= 0

    ls_ph = list()
    for f in intF:
        if f > f_min and f < f_max:
            c = np.log(f / (int_max - f))
            ls_ph.append(round(pKa + c, 5))
        elif f < f_min:
            ls_ph.append(0)
        else:
            ls_ph.append(14)
    return ls_ph


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


def _tan_simulation(c_nh4, phmin, phmax, step_ph, ph_deci, pKa=9.25):
    pH = np.linspace(phmin, phmax, num=int((phmax-phmin)/step_ph+1))
    # NH3 and NH4 percentage depending on pH
    alpha = partial_concNH3(c_nh4=c_nh4*1e6, pH=pH, pKa=pKa)
    df_alpha = pd.DataFrame(alpha, index=pH, columns=['nh3 %'])
    xnew = [round(i, ph_deci) for i in df_alpha.index]
    df_alpha.index = xnew

    # general description of NH3 vs pH curve
    nh4_pc = 1 - df_alpha
    nh4_pc.columns = ['nh4 %']
    df_alpha['nh4 %'] = nh4_pc
    return df_alpha


def _alpha4ph(df_pH, df_alpha):
    for i in df_pH.index:
        df_pH.loc[i, 'alpha_nh4 [%]'] = df_alpha['nh4 %'].loc[df_pH.loc[i, 'pH']]
        df_pH.loc[i, 'alpha_nh3 [%]'] = df_alpha['nh3 %'].loc[df_pH.loc[i, 'pH']]
    df_pH.index.name = 'time / s'
    return df_pH


# --------------------------------------------------------------------------------------------------------------------
def move_figure(xnew, ynew):
    mngr = plt.get_current_fig_manager()
    geom = mngr.window.geometry()
    x, y, dx, dy = geom.getRect()
    mngr.window.setGeometry(xnew, ynew, dx, dy)


def plot_sensorresponse(dfconc, sens_time0, tstart_sens, step, D, par_meas, par_sens, arg_fig, plotCheck=True):
    if 'figsize' in arg_fig.keys():
        figsize = arg_fig['figsize']
    else:
        figsize = (6, 3)
    fig, ax = plt.subplots(figsize=figsize)
    if 'title' in arg_fig.keys():
        fig.canvas.set_window_title(arg_fig['title'])
    if 'fontsize' in arg_fig.keys():
        fs = arg_fig['fontsize']
    else:
        fs = 10
    ax.set_xlabel('Time / s', fontsize=fs), ax.set_ylabel('Analyte concentration', fontsize=fs)
    ax.set_ylim(-0.05, par_sens['max conc']*1.19)

    # ..............................................................
    # target concentration
    ax.plot(dfconc, ':k', label='target concentration')

    # time stamps when conc is changing
    n, tsample, t = 0, list(), 0
    while t < par_meas['stop']:
        t = par_meas['start'] + n * D
        tsample.append(par_meas['start'] + n * D)
        n += 1
    # ..............................................................
    # sensor response
    # each time I change the concentration, the sensor has to respond to the change
    # 1) sensor is delayed according to the sampling rate
    # 2) response does not equal target concentration as sensor has a response time
    ddata, c_apparent = dict(), 0
    for en, t in enumerate(tsample):
        # closest time step in concentration fluctuation
        tc = min(dfconc.index, key=lambda x: abs(x - t))
        if tc != dfconc.index[-1]:
            # find position in list for previous concentration
            pos = [en for en in range(len(dfconc.index)) if dfconc.index[en] == tc][0]
            tc_1, tc1 = dfconc.index[pos - 1], dfconc.index[pos + 1]

            # update actual sensor measurement time +t
            if en == 0:
                sens_time = sens_time0
            else:
                sens_time = np.linspace(sens_time[-1], tc + D, num=int((D-tstart_sens)/step + 1))
            # sensor signal response
            df = sensor_response(t=tc1, t_1=tc_1, dfconc=dfconc, par_sens=par_sens, sens_time0=sens_time0,
                                 c_apparent=c_apparent)
            df.index = sens_time
            ddata[tc], c_apparent = df, df['signal'].to_numpy()[-1]
            ax.plot(sens_time, df, lw=1.25, color='teal', label='sensor response')

            if plotCheck is True:
                # check sensor response time
                ax.axvline(par_sens['response time'], lw=0.5, color='crimson')
                ax.axhline(0.9 * par_sens['max conc'], lw=0.5, color='crimson')

        for en, t in enumerate(tsample):
            # arrows indicating that the sensor starts to react to the surrounding
            ax.arrow(t, par_sens['max conc']*1.15, 0.0, (par_sens['max conc']-par_sens['max conc']*1.15)/2.5, fc="k",
                     ec="k", head_width=.8, head_length=0.1)

    ax.legend(['target concentration', 'sensor response'], loc=0, fontsize=fs * 0.7)
    plt.tight_layout(), sns.despine(), plt.grid(False)

    # .............................................................................
    # rearrange data dictionary
    tstarts = list(ddata.keys())
    df_data = pd.concat(ddata, axis=0)
    xnew = df_data.index.levels[1]

    ls_time = list(df_data.loc[tstarts[0]].index)
    for i in range(len(tstarts) - 1):
        ls_time.extend(df_data.loc[tstarts[i + 1]].index)

    dfdata = pd.DataFrame(df_data['signal'].to_numpy(), ls_time, columns=['signal'])
    dfdata.index.name = 'time / s'

    return fig, ax, dfdata


def plot_2sumresponse(sens_time0, par_target, par_sens1, par_meas, arg_fig, par_sens2=None, plotCheck=False):
    # individual sensor responses
    plt.ioff()
    if isinstance(par_target['conc1'], pd.DataFrame):
        fig1, ax, dfdata1 = plot_sensorresponse(dfconc=par_target['conc1'], sens_time0=sens_time0, par_meas=par_meas,
                                                par_sens=par_sens1, step=par_meas['steps'], D=par_target['pulse width'],
                                                tstart_sens=par_meas['start sensing'], arg_fig=arg_fig,
                                                plotCheck=plotCheck)
    else:
        fig1, dfdata1 = None, None
    if isinstance(par_target['conc2'], pd.DataFrame) and par_sens2:
        fig2, ax, dfdata2 = plot_sensorresponse(dfconc=par_target['conc2'], sens_time0=sens_time0, par_meas=par_meas,
                                                par_sens=par_sens2, step=par_meas['steps'], D=par_target['pulse width'],
                                                tstart_sens=par_meas['start sensing'], arg_fig=arg_fig,
                                                plotCheck=plotCheck)
    else:
        fig2, dfdata2 = None, None
    plt.close(fig1), plt.close(fig2)

    # ----------------------------------------------------------------------------------------
    figTAN, axTAN = plt.subplots(figsize=arg_fig['figsize'], nrows=2, sharex=True)
    figTAN.canvas.set_window_title(arg_fig['title'])

    axTAN[1].set_xlabel('Time / s', fontsize=arg_fig['fontsize']),
    axTAN[0].set_ylabel('Ind. Analyte concentration', fontsize=arg_fig['fontsize'])
    axTAN[1].set_ylabel('sum concentration', fontsize=arg_fig['fontsize'])

    # individual sensor signal
    axTAN[0].plot(pd.concat([dfdata1, dfdata2], axis=1), lw=1.)
    if isinstance(par_target['conc1'], pd.DataFrame):
        axTAN[0].plot(par_target['conc1'], lw=1., ls=':', color='gray')
    if isinstance(par_target['conc2'], pd.DataFrame):
        axTAN[0].plot(par_target['conc2'], lw=1., ls='--', color='gray')

    # sum parameter compared to target concentration
    tan = pd.DataFrame(pd.concat([dfdata1, dfdata2], axis=1, ignore_index=True)).sum(axis=1)
    dfconc_sum = par_target['conc1'] + par_target['conc2']
    axTAN[1].plot(tan, lw=0., marker='.', ms=2, color='#C51D74')
    axTAN[1].plot(dfconc_sum, lw=1., ls=':', color='k')
    axTAN[0].legend(['sensor1 - response:{:.1f}s'.format(par_sens1['response time']),
                     'sensor2 - response:{:.1f}s'.format(par_sens2['response time']),
                     'target conc1', 'target conc2'], fontsize=arg_fig['fontsize']*0.7)#, loc='center right')
    axTAN[1].legend(['sum conc.', 'target total conc.'], fontsize=arg_fig['fontsize'] * 0.7, loc='center right')

    sns.despine(), plt.tight_layout()

    return figTAN, axTAN, tan, dfconc_sum


def plot_tanSimulation(df_alpha, phmax=14, pKa=9.25, xnew=50, ynew=60, figsize=(5, 3)):
    fig, ax = plt.subplots(figsize=figsize)
    move_figure(xnew=xnew, ynew=ynew)
    fig.canvas.set_window_title('Ammonia vs. Ammonium')
    ax.set_xlabel('pH', fontsize=fs), ax.set_ylabel('alpha [%]', fontsize=fs)

    ax.plot(df_alpha['nh4 %'], color=dcolor['NH4'], lw=1., label='NH$_4^+$')
    ax.plot(df_alpha['nh3 %'], color=dcolor['NH3'], lw=1., label='NH$_3$')

    ax.axvline(pKa, color='k', ls=':', lw=1.)
    ax.axhline(0.5, color='k', ls=':', lw=1.)
    ax.legend(loc='upper center', bbox_to_anchor=(1, 0.9), frameon=True, fancybox=True, fontsize=fs * 0.7)

    ax.set_xlim(-0.5, phmax * 1.05)
    sns.despine(), plt.tight_layout()
    return fig


def calibration_ph(dfSig_calib, xnew=50, ynew=450, figsize=(5, 3)):
    fig, ax = plt.subplots(figsize=figsize)
    move_figure(xnew=xnew, ynew=ynew)

    fig.canvas.set_window_title('Calibration pH sensor (electrochemical)')
    ax.set_xlabel('pH', fontsize=fs), ax.set_ylabel('Potential / mV', fontsize=fs)

    ax.plot(dfSig_calib, lw=1., color='k')
    sns.despine(), plt.tight_layout()
    return fig


def calibration_nh3(conc_nh3, para_nh3, xnew=50, ynew=500, figsize=(5, 3)):
    fig, ax = plt.subplots(figsize=figsize)
    move_figure(xnew=xnew, ynew=ynew)

    fig.canvas.set_window_title('Calibration NH3 sensor (electrochemical)')
    ax.set_xlabel('alpha(NH$_3$) / %', fontsize=fs * 0.9), ax.set_ylabel('Sensor signal / mV', fontsize=fs * 0.9)

    ax.plot(conc_nh3, para_nh3[0] * conc_nh3 + para_nh3[1], lw=1., color='k')

    sns.despine(), plt.tight_layout()
    return fig


def plot_tanModel(dfph_re, dfconc_target, df_record, df_tan_target, df_tan, phmax, xnew=690, ynew=350, figsize=(6,4.5)):
    fig = plt.figure(figsize=figsize)
    move_figure(xnew=xnew, ynew=ynew)
    fig.canvas.set_window_title('NH3 / NH4+ simulation')

    gs = GridSpec(nrows=3, ncols=1)
    ax = fig.add_subplot(gs[:2, 0])
    ax_ = fig.add_subplot(gs[2, 0], sharex=ax)
    ax_ph = ax_.twinx()
    ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False)

    ax_.set_xlabel('Time / s'), ax_ph.set_ylabel('pH', color='gray')
    ax_.set_ylabel('TAN / ppm', color=dcolor['TAN'])

    # top plot
    ax.plot(df_record['nh3 / ppm'], lw=1., color=dcolor['NH3'], label='NH$_3$')
    ax.plot(df_record['nh4 / ppm'], lw=1., color=dcolor['NH4'], label='NH$_4^+$')
    ax.legend(frameon=True, fancybox=True, loc=0, fontsize=fs * 0.7)
    ax.plot(dfconc_target['nh3 / ppm'], lw=1., ls=ls['target'], color='k')
    ax.plot(dfconc_target['nh4 / ppm'], lw=1., ls=ls['target'], color='gray')
    ax.set_ylim(-10, dfconc_target['nh3 / ppm'].max()*1.1)

    # bottom plot
    ax_.plot(df_tan_target, lw=1., ls=ls['target'], color='k')
    ax_.plot(df_tan, lw=1., color=dcolor['TAN'])
    ax_ph.plot(dfph_re, lw=1., ls='-.', color='gray')
    ax_.set_ylim(df_tan['TAN'].min()*0.95, df_tan['TAN'].max()*1.05), ax_ph.set_ylim(-0.5, phmax * 1.05)

    plt.tight_layout()
    return fig


# --------------------------------------------------------------------------------------------------------------------
def sensor_response(t, t_1, dfconc, par_sens, sens_time0, c_apparent):
    if dfconc.loc[t_1, 'signal'] == 0 and dfconc.loc[t, 'signal'] == 0 and c_apparent < par_sens['max conc']:
        # print('case1 - increase', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = _sensor_response(t=sens_time0, par_sens=par_sens, y=0, direction='increase')[0] + c_apparent

    elif dfconc.loc[t_1, 'signal'] == 0 and dfconc.loc[t, 'signal'] == 0 and c_apparent == par_sens['max conc']:
        # print('case2', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = pd.DataFrame([c_apparent]*len(sens_time0), index=sens_time0, columns=['signal'])

    elif dfconc.loc[t_1, 'signal'] == par_sens['max conc'] and dfconc.loc[t, 'signal'] == 0 \
            and c_apparent < par_sens['max conc']/2:
        # print('case3 - decline', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = _sensor_response(t=sens_time0, par_sens=par_sens, direction='decline')[0]

    elif dfconc.loc[t_1, 'signal'] == par_sens['max conc'] and dfconc.loc[t, 'signal'] == 0 \
            and c_apparent > par_sens['max conc']/2:
        # print('case4 - decline', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = _sensor_response(t=sens_time0, par_sens=par_sens, y=c_apparent, direction='decline')[0]

    elif dfconc.loc[t_1, 'signal'] == par_sens['max conc'] and dfconc.loc[t, 'signal'] == 0 \
            and c_apparent == par_sens['max conc']:
        # print('case5', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = pd.DataFrame([c_apparent]*len(sens_time0), index=sens_time0, columns=['signal'])

    elif dfconc.loc[t_1, 'signal'] == 0 and dfconc.loc[t, 'signal'] == par_sens['max conc'] \
            and c_apparent < par_sens['max conc']:
        # print('case6 - increase', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = _sensor_response(t=sens_time0, par_sens=par_sens, direction='increase')[0] +\
             c_apparent

    elif dfconc.loc[t_1, 'signal'] == par_sens['max conc'] and dfconc.loc[t, 'signal'] == par_sens['max conc']\
            and c_apparent == par_sens['max conc']:
        # print('case7', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = pd.DataFrame([c_apparent]*len(sens_time0), index=sens_time0, columns=['signal'])

    elif dfconc.loc[t_1, 'signal'] == par_sens['max conc'] and dfconc.loc[t, 'signal'] == par_sens['max conc']\
            and c_apparent < par_sens['max conc']/2:
        # print('case8 - increase', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = _sensor_response(t=sens_time0, par_sens=par_sens, y=c_apparent, direction='increase')[0]

    elif dfconc.loc[t_1, 'signal'] == par_sens['max conc'] and dfconc.loc[t, 'signal'] == par_sens['max conc']\
            and c_apparent > par_sens['max conc']/2:
        # print('case9 - decline', dfconc.loc[t_1, 'signal'], dfconc.loc[t, 'signal'], c_apparent)
        df = _sensor_response(t=sens_time0, par_sens=par_sens, y=c_apparent, direction='decline')[0]

    else:
        raise ValueError('case undefined at {:.4f}: c(t)={:.4f}, c(t-1)={:.4f}, '
                         'c_apparent={:.4f}'.format(t, dfconc.loc[t, 'signal'], dfconc.loc[t_1, 'signal'], c_apparent))
    return df


def _pHsensor_response(t_plateau, pH_plateau, dfpH_target, t90_pH, ph_res, sig_bgd, step=0.01):
    dsig = dict()
    for n in range(len(pH_plateau)):
        if n == 0:
            y_1 = sig_bgd
        else:
            y_1 = dsig[n - 1]['signal / mV'].to_numpy()[-1]
        sdiff = dfpH_target['target potential / mV'].loc[n * t_plateau] - y_1

        if n == 0:
            # 1st sensor response - always different
            t_sens = np.arange(0, t_plateau - 1, step)
            y0 = _gompertz_curve_v1(x=t_sens, t90=t90_pH, tau=ph_res, pstart=sig_bgd, slope='increase', s_diff=sdiff)
            dfSig = pd.DataFrame(y0, index=t_sens, columns=['signal / mV'])
        else:
            # other sensor responses behave similar
            t_sens_ = np.arange(0, t_plateau + step, step)

            # different cases - decline or increase
            if sdiff < 0:
                y1 = _gompertz_curve_v1(x=t_sens_, t90=t90_pH, tau=ph_res, pstart=y_1, slope='decline',
                                        s_diff=np.abs(sdiff))
            else:
                y1 = _gompertz_curve_v1(x=t_sens_, t90=t90_pH, tau=ph_res, pstart=y_1, slope='increase',
                                        s_diff=np.abs(sdiff))
            t_sens1 = t_sens_ + n * t_plateau + (n - 1) * step - 1
            dfSig = pd.DataFrame(y1, index=t_sens1, columns=['signal / mV'])
        dsig[n] = dfSig

    sens_response = pd.concat(dsig)
    xnew = sens_response.index.levels[1]
    sens_response.index = xnew
    return sens_response


def _potential2pH(sens_response, ph_deci, E0, T=25):
    y = Nernst_equation_invert(E=sens_response['signal / mV'], T=T, E0=E0)
    dfph_re = pd.DataFrame([round(i, ph_deci) for i in y], index=sens_response.index)
    dfph_re.columns = ['pH recalc.']
    return dfph_re


def _calibration_nh3(anh3_min, anh3_max, anh3_step, sigNH3_bgd, sigNH3_max):
    # linear regression for all pH values
    nh3_cali = lin_regression(cnh3_min=anh3_min, cnh3_max=anh3_max, conc_step=anh3_step, sigNH3_bgd=sigNH3_bgd,
                              sigNH3_max=sigNH3_max)
    return nh3_cali


def _alpha_vs_time(df_alpha, dfpH_target):
    # specific target concentration (NH3) over time in percentages
    c_nh3 = pd.DataFrame(df_alpha['nh3 %'].loc[dfpH_target['pH'].to_numpy()])
    c_nh4 = pd.DataFrame(df_alpha['nh4 %'].loc[dfpH_target['pH'].to_numpy()])
    c_nh3.index, c_nh4.index = dfpH_target['pH'].index, dfpH_target['pH'].index

    # DataFrame for NH3/NH4 concentration in percent over time (according to respective pH)
    dfalpha_target = pd.concat([c_nh3, c_nh4, dfpH_target['pH']], axis=1)
    return c_nh3, c_nh4, dfalpha_target


def conv2potent_nh3(dfalpha_target, c_nh3, c_nh4, dfpH_target, para_nh3):
    # target concentration over time in ppm
    dfconc = pd.concat([dfalpha_target['nh3 %'] * c_nh3, dfalpha_target['nh4 %'] * c_nh4,
                        dfpH_target['pH']], axis=1)
    dfconc.columns = ['nh3 / ppm', 'nh4 / ppm', 'pH']

    # target potential over time in mV
    dfpotent = pd.DataFrame(para_nh3[0] * dfconc['nh3 / ppm'] + para_nh3[1])
    dfpotent.columns, dfpotent.index.name = ['nh3 / mV'], 'time / s'
    dfpotent['pH'] = dfpH_target['pH']

    return dfconc, dfpotent


def _nh3sensor_response(t_plateau, pH_plateau, dfph_re, t90_nh3, nh3_res, dfpot_target, step=.01):
    dsig_nh3 = dict()
    for n in range(len(pH_plateau)):
        if n == 0:
            y_1 = dfpot_target.loc[0, 'nh3 / mV']
        else:
            y_1 = dsig_nh3[n - 1]['signal / mV'].to_numpy()[-1]
        sdiff = dfpot_target['nh3 / mV'].loc[n * t_plateau] - y_1

        if n == 0:
            # 1st sensor response - always different
            t_sens = np.arange(0, t_plateau - 1, step)
            # apparent pH
            ph_i = dfph_re.loc[t_plateau]
            y1 = _gompertz_curve_v1(x=t_sens, t90=t90_nh3, tau=nh3_res, slope='decline',  s_diff=sdiff,
                                    pstart=dfpot_target.loc[0, 'nh3 / mV'])
            dfSig = pd.DataFrame(y1, index=t_sens, columns=['signal / mV'])
        else:
            # other sensor responses behave similar
            t_sens_ = np.arange(0, t_plateau + step, step)

            # different cases - decline or increase
            if sdiff < 0:
                y1 = _gompertz_curve_v1(x=t_sens_, t90=t90_nh3, tau=nh3_res, pstart=y_1, slope='decline',
                                        s_diff=np.abs(sdiff))
            else:
                y1 = _gompertz_curve_v1(x=t_sens_, t90=t90_nh3, tau=nh3_res, pstart=y_1, slope='increase',
                                        s_diff=np.abs(sdiff))
            t_sens1 = t_sens_ + n * t_plateau + (n - 1) * step - 1
            dfSig = pd.DataFrame(y1, index=t_sens1, columns=['signal / mV'])
        dsig_nh3[n] = dfSig

    sensNH3_response = pd.concat(dsig_nh3)
    xnew = sensNH3_response.index.levels[1]
    sensNH3_response.index = xnew
    sensNH3_response[dfph_re.columns[0]] = pd.DataFrame(dfph_re, index=xnew)
    return sensNH3_response


def _potent2nh3(sensNH3_response, para_nh3):
    # apparent pH and respective calibration parameter
    df = pd.DataFrame((sensNH3_response['signal / mV'] - para_nh3[1]) / para_nh3[0], index=sensNH3_response.index)
    df.columns = ['nh3 / ppm']
    return df


# --------------------------------------------------------------------------------------------------------------------
def pH_sensor(sensor_ph, para_meas, df_alpha):
    # unpacking dictionaries
    E0, t90_pH, pHres, Tres_ph = sensor_ph['E0'], sensor_ph['t90'], sensor_ph['resolution'], sensor_ph['time steps']
    pHsens, sig_bgd = sensor_ph['sensitivity'], sensor_ph['background signal']
    Temp, t_steady, pH_steps = para_meas['Temperature'], para_meas['Plateau time'], para_meas['pH steps']

    # ------------------------------------------------------------------------------------------------------------------
    # electrochemical pH sensor calibration; source vernier.com/2018/04/06/the-theory-behind-ph-measurements/
    dfSig_calib = pd.DataFrame(Nernst_equation(ls_ph=np.array(pH_steps), T=Temp, E0=E0), index=pH_steps,
                               columns=['signal / mV'])
    dfSig_calib.index.name = 'pH'

    # pH changes to target potential over time
    # create time range
    l = [[ph_i] * t_steady for ph_i in pH_steps]
    ls_pH = l[0]
    [ls_pH.extend(l_i) for l_i in l[1:]]

    # finalize target pH and concentrations
    df_ph = pd.DataFrame(ls_pH, columns=['pH'])
    df_pH = _alpha4ph(df_pH=df_ph, df_alpha=df_alpha)

    dfSig_steps = pd.DataFrame(Nernst_equation(ls_ph=np.array(df_pH['pH'].to_numpy()), T=Temp, E0=E0),
                               index=df_pH['pH'].index, columns=['target potential / mV'])
    dfpH_target = pd.concat([df_pH['pH'], dfSig_steps], axis=1)

    # include sensor response
    sens_response = _pHsensor_response(t_plateau=t_steady, pH_plateau=pH_steps, ph_res=pHres, dfpH_target=dfpH_target,
                                       t90_pH=t90_pH, step=Tres_ph, sig_bgd=sig_bgd)

    # recalculation of pH values
    dfph_re = _potential2pH(sens_response=sens_response, ph_deci=pHsens, E0=E0, T=Temp)

    return dfpH_target, dfSig_calib, dfph_re


def NH3_sensor(df_alpha, dfph_re, dfpH_target, sensor_nh3, para_meas):
    # unpacking dictionaries
    t90_nh3, nh3res, Tsteps_nh3 = sensor_nh3['response time'], sensor_nh3['resolution'], sensor_nh3['time steps']
    Temp, t_steady, pH_steps = para_meas['Temperature'], para_meas['Plateau time'], para_meas['pH steps']

    pHrange, ph_deci, pKa = sensor_nh3['pH range'], sensor_nh3['sensitivity'], sensor_nh3['pKa']
    nh3_range, signal_min, signal_max = sensor_nh3['nh3 range'], sensor_nh3['signal min'], sensor_nh3['signal max']

    # ------------------------------------------------------------------------------------------------------------------
    # linear regression for all pH values
    para_nh3 = _calibration_nh3(anh3_min=nh3_range[0], anh3_max=nh3_range[1], anh3_step=nh3_range[2],
                                sigNH3_bgd=signal_min, sigNH3_max=signal_max)

    # target NH3 signal for NH3 and NH4+ according to the measured pH - specific target concentration over time in %
    nh3_target, nh4_target, dfalpha_target = _alpha_vs_time(df_alpha=df_alpha, dfpH_target=dfpH_target)

    if 'target concentration' in para_meas:
        cnh3_target = para_meas['target concentration']
    else:
        # assume a constant ammonia concentration over time
        cnh3_target = [para_meas['GGW concentration']]*len(dfalpha_target.index)

    [dfconc_target,
     dfpot_target] = conv2potent_nh3(dfalpha_target=dfalpha_target, c_nh3=cnh3_target, c_nh4=cnh3_target,
                                     dfpH_target=dfpH_target, para_nh3=para_nh3)

    # include sensor response - double checked. NH3 starting with ideal situation
    sensNH3_resp = _nh3sensor_response(t_plateau=t_steady, pH_plateau=pH_steps, step=Tsteps_nh3, dfph_re=dfph_re,
                                       t90_nh3=t90_nh3, nh3_res=nh3res, dfpot_target=dfpot_target)

    # recalculation of NH3 concentrations while considering the pH(recalc.)
    dfnh3_calc = _potent2nh3(sensNH3_response=sensNH3_resp, para_nh3=para_nh3)

    # get NH4 concentration via Henderson-Hasselbalch
    dfnh4_calc = pd.DataFrame(henderson_nh4(pKa=pKa, pH=dfph_re['pH recalc.'].to_numpy(),
                                            c_nh3=dfnh3_calc['nh3 / ppm'].to_numpy()), columns=['nh4 / ppm'],
                              index=dfnh3_calc.index)

    # include sampling rate for actually reading out sensor data
    x_smg = np.arange(dfnh3_calc.index[0], dfnh3_calc.index[-1], para_meas['sampling rate'])
    x_smg = [i.round(2) for i in x_smg]

    df_record = pd.concat([dfnh3_calc, dfnh4_calc], axis=1)
    tnew = [round(i, 2) for i in df_record.index]
    df_record.index = tnew
    df_record = df_record.loc[x_smg]

    # Calculation total ammonia nitrogen: TAN = NH3 + NH4+
    df_tan_target = pd.DataFrame(dfconc_target[['nh3 / ppm', 'nh4 / ppm']].sum(axis=1), columns=['TAN'])
    df_tan = pd.DataFrame(df_record.dropna().sum(axis=1), columns=['TAN'])
    return df_tan_target, dfconc_target, para_nh3, df_record, df_tan


def tan_simulation(sensor_ph, para_meas, sensor_nh3, plot='result'):
    # unpacking dictionaries
    pHrange, ph_deci, pKa = sensor_nh3['pH range'], sensor_nh3['sensitivity'], sensor_nh3['pKa']
    nh3_range, signal_min, signal_max = sensor_nh3['nh3 range'], sensor_nh3['signal min'], sensor_nh3['signal max']

    # ------------------------------------------------------------------------------------------------------------------
    # TAN system - pH curve of NH3 and NH4 in percentage
    df_alpha = _tan_simulation(c_nh4=para_meas['GGW concentration'], phmin=pHrange[0], phmax=pHrange[1],
                               step_ph=pHrange[2], ph_deci=ph_deci, pKa=pKa)

    # ------------------------------------------------------------------------------------------------------------------
    # pH sensor modelation
    dfpH_target, dfSig_calib, dfph_re = pH_sensor(sensor_ph=sensor_ph, para_meas=para_meas, df_alpha=df_alpha)

    # ------------------------------------------------------------------------------------------------------------------
    # NH3 / NH4+ sensor modulation
    [df_tan_target, dfconc_target, para_nh3,
     df_record, df_tan] = NH3_sensor(df_alpha=df_alpha, dfph_re=dfph_re, dfpH_target=dfpH_target, sensor_nh3=sensor_nh3,
                                     para_meas=para_meas)

    # ------------------------------------------------------------------------------------------------------------------
    # plotting part
    # Turn interactive plotting off
    plt.ioff()

    # system simulation
    fig = plot_tanSimulation(df_alpha=df_alpha, phmax=14, pKa=9.25)

    # single sensor calibration
    fig1 = calibration_ph(dfSig_calib=dfSig_calib)
    conc_nh3 = np.arange(nh3_range[0], nh3_range[1], step=nh3_range[2])
    fig2 = calibration_nh3(conc_nh3, para_nh3, xnew=50, ynew=500, figsize=(5, 3))

    # final model
    fig3 = plot_tanModel(dfph_re=dfph_re, dfconc_target=dfconc_target, df_record=df_record, df_tan_target=df_tan_target,
                         df_tan=df_tan, phmax=sensor_nh3['pH range'][1])

    # depending on user input close certain figure plots
    if plot == 'calib':
        plt.close(fig3)
    elif plot == 'results':
        plt.close(fig), plt.close(fig1), plt.close(fig2)
    elif plot == 'all':
        pass
    # Display all "open" (non-closed) figures
    plt.show()

    # ------------------------------------------------------------------------------------------------------------------
    # collect for result output
    dic_target = dict({'TAN': df_tan_target, 'NH3 simulation': df_alpha, 'target conc nh3': dfconc_target})
    dic_sens_calib = dict({'pH': dfSig_calib, 'NH3': para_nh3})
    dic_sens_record = dict({'tan': df_tan, 'NH3': df_record, 'pH': dfph_re})
    dic_figures = dict({'tan simulation': fig, 'pH calib': fig1, 'NH3 calib': fig2, 'model': fig3})

    return dic_target, dic_sens_calib, dic_sens_record, dic_figures, df_tan

