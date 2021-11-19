__author__ = 'Silvia E Zieger'
__project__ = 'sensor response in silico'

"""Copyright 2021. All rights reserved.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable 
for any damages arising from the use of this software.
Permission is granted to anyone to use this software within the scope of evaluating mutli-analyte sensing. No permission
is granted to use the software for commercial applications, and alter it or redistribute it.

This notice may not be removed or altered from any distribution.
"""

import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
import pandas as pd


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


# --------------------------------------------------------------------------------------------------------------------
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
    plt.show()

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
    plt.show()

    return figTAN, axTAN, tan, dfconc_sum


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

