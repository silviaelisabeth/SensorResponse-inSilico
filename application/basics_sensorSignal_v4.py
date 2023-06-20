__author__ = 'Silvia E Zieger'
__project__ = 'sensor response in silico'

"""Copyright 2021. All rights reserved.

This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable 
for any damages arising from the use of this software.
Permission is granted to anyone to use this software within the scope of evaluating mutli-analyte sensing. No permission
is granted to use the software for commercial applications, and alter it or redistribute it.
# 
This notice may not be removed or altered from any distribution.
"""

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------------------------------------------------
# global parameter
n = 1
R = 8.314         # [J/mol*K]
F = 96485         # [C/mol]
s_nernst = 2.303
Tconvert = 273.15
temp_std = 25

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


def _calibration_electroSens(E, conc):
    # Nernst equation with slope -59mV: E = E0 - 59mV/z * log10(c)
    E0 = E + (59/n * np.log(conc))
    para = dict({'slope': 59, 'E0': E0})
    return para


def henderson_base(pH, c_acid, pKa=9.25):
    """ Returns the NH3 concentration depending on the pH and the concentrations of NH4+ """
    return c_acid * 10 ** (pH - pKa)


def henderson_acid(pH, c_base, pKa=9.25):
    """ Returns the NH3 concentration depending on the pH and the concentrations of NH4+ """
    return c_base * 10 ** (pKa - pH)


def henderson_Sum4Base(pH, c_sum, pKa=9.25):
    """
    Returns the NH3 concentration depending on the pH and the TAN concentration
    :param pH:
    :param c_sum:
    :param pKa:
    :return:
    """
    return c_sum / (1 + 10**(pKa - pH))


def henderson_Sum4Acid(pH, c_sum, pKa=9.25):
    """
    Returns the NH4+ concentration depending on the pH and the TAN concentration
    :param pH:
    :param c_sum:
    :param pKa:
    :return:
    """
    return c_sum / (1 + 10**(pH - pKa))


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


def _sensor_response(cplateau, psensor, tplateau, cstart):
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

        df_sig = pd.DataFrame(_gompertz_curve_v1(x=sens_time, t90=psensor['t90'], s_diff=sdiff, pstart=c_apparent,
                                                 tau=psensor['resolution'], slope='increase'), index=xnew,
                              columns=['potential mV'])

        start, c_apparent = df_sig['potential mV'].index[-1], df_sig['potential mV'].to_numpy()[-1]
        dsig[en] = df_sig
        df_sensor = pd.concat(dsig)
        df_sensor.index = [i[1] for i in df_sensor.index]

    return df_sensor


def _all_plateaus_inlist(ls):
    ls_ = list()
    for xi in [x for i,x in enumerate(np.diff(ls)) if x != 0]:
        ls_ += (list(np.where(np.diff(ls) == xi)[0]))

    pos_jump = [0] + sorted(list(dict.fromkeys(ls_))) + [len(ls)]
    plateau = [np.mean(ls[pos_jump[en] + 1:pos_jump[en+1]]) for en in range(len(pos_jump)-1)]

    return plateau


# --------------------------------------------------------------------------------------------------------------------
def calibSensor_pH(target_ph, sensor_ph):
    # pH sensor - targeted pH fluctuation in mV
    df_sigpH_mV = pd.DataFrame(Nernst_equation(ls_ph=target_ph['target pH'].to_numpy(), T=temp_std,
                                               E0=sensor_ph['E0']), index=target_ph.index, columns=['potential mV'])
    # what are the specific signal levels
    cplateau = [Nernst_equation(ls_ph=p, T=temp_std, E0=sensor_ph['E0']) for p in sensor_ph['set values']]
    return df_sigpH_mV, cplateau


def calibSensor_para(target_para, sensor2):
    # sensor 2 - targeted concentration fluctuation in mV
    para = _calibration_electroSens(E=sensor2['E'], conc=sensor2['conc_calib'])

    df_sigSens2_mV = pd.DataFrame(target_para['target_mg/L Sum'] * para['slope'] + para['E0'])
    df_sigSens2_mV.columns = ['potential mV Sum']

    # what are the specific signal levels
    cplateau = [t * para['slope'] + para['E0'] for t in sensor2['set values']]
    return df_sigSens2_mV, cplateau, para


def _alignSensorSettings(df_target, sensor_ph, sensor2):
    # individual sensor calibration
    df_sigpH_mV, cplateaupH = calibSensor_pH(target_ph=df_target, sensor_ph=sensor_ph)
    df_sigPara_mV, cplateauPara, para_para = calibSensor_para(target_para=df_target, sensor2=sensor2)

    df_ = pd.concat([df_sigpH_mV, df_sigPara_mV], axis=1)
    df_.columns = ['Potential mV pH', 'Potential mV Sum']
    df_res = pd.concat([df_target, df_], axis=1)
    return df_res, cplateaupH, cplateauPara, para_para


def pH_sensor(cplateau, sensor_ph, para_meas):
    # include sensor response
    cstart = Nernst_equation(ls_ph=sensor_ph['start pH'], E0=sensor_ph['E0'], T=temp_std)
    df_pHrec = _sensor_response(cplateau=cplateau, psensor=sensor_ph, tplateau=para_meas['plateau time'], cstart=cstart)

    # -------------------------------------
    # re-calculate pH from potential
    df_recalc = pd.DataFrame(Nernst_equation_invert(E=df_pHrec['potential mV'], E0=sensor_ph['E0'], T=temp_std))
    df_recalc.columns = ['pH calc']
    df_recalc.index = [round(i, 2) for i in df_recalc.index]

    return df_pHrec, df_recalc


def para2_sensors(cplateau, sensor2, para_meas):
    # update concentration in case one analyte changes more often than the other one
    # include sensor response
    df2_rec = _sensor_response(cplateau=cplateau, psensor=sensor2, tplateau=para_meas['plateau time'],
                               cstart=cplateau[0])
    return df2_rec


def para2Sensorcalc(df_res, analyte, cplateauSum, df_pHcalc, sensor2, para_meas):
    if analyte == 'base' or analyte == 'Base':
        cBase = [df_res[df_res['Potential mV Sum'] == c]['target_mg/L Base'].to_numpy() for c in cplateauSum]
        ls_cBase = list()
        for li in cBase:
            ls_cBase = np.append(ls_cBase, li)
        # find all plateaus even though they repeat themselves
        ls_cBase = _all_plateaus_inlist(ls_cBase)[:len(cplateauSum)]

        df_base_calc = _sensor_response(cplateau=ls_cBase, psensor=sensor2, tplateau=para_meas['plateau time'],
                                        cstart=ls_cBase[0])
        df_base_calc.columns = ['Base calc']

        # calculate acid species from Henderson-Hasselbalch and sensed pH
        df_acid_calc = pd.DataFrame(henderson_acid(pH=df_pHcalc['pH calc'].to_numpy(), pKa=sensor2['pKa'],
                                                   c_base=df_base_calc['Base calc'].to_numpy()), columns=['Acid calc'],
                                    index=[round(i, 2) for i in df_pHcalc.index])
    elif analyte == 'acid' or analyte == 'Acid':
        cAcid = [df_res[df_res['Potential mV Sum'] == c]['target_mg/L Acid'].to_numpy() for c in cplateauSum]

        ls_cAcid = list()
        for li in cAcid:
            ls_cAcid = np.append(ls_cAcid, li)
        # find all plateaus even though they repeat themselves
        ls_cAcid = _all_plateaus_inlist(ls_cAcid)[:len(cplateauSum)]

        df_acid_calc = _sensor_response(cplateau=ls_cAcid, psensor=sensor2, tplateau=para_meas['plateau time'],
                                        cstart=ls_cAcid[0])

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


def save_integral(para_meas, sensor_ph, sensor_para2, dres, total_lbl):
    # combine measurement parameter
    ind = ['plateau time (s)', 'pKa', 'analyte']
    df_meas1 = pd.DataFrame([para_meas['plateau time'], '-- ', 'pH'], index=ind)
    df_meas2 = pd.DataFrame([para_meas['plateau time'], sensor_para2['pKa'], sensor_para2['analyte']], index=ind)
    df_total = pd.DataFrame([para_meas['plateau time'], sensor_para2['pKa'], total_lbl], index=ind)

    df_meas = pd.concat([df_meas1, df_meas2, df_total], axis=1)
    df_meas.loc['run simulation', :] = df_meas.columns
    df_meas.columns = np.arange(len(df_meas.columns))

    # combine sensor settings (pH, sensor 2, and total parameter)
    df_pH = pd.DataFrame([str(sensor_ph['set values'])[1:-1], sensor_ph['t90']])
    df_para2 = pd.DataFrame([str([round(i, 2) for i in sensor_para2['set sensor2']])[1:-1], sensor_para2['t90']])
    df_sum = pd.DataFrame([str(sensor_para2['set values'])[1:-1], '-- '])

    df_settings = pd.concat([df_pH, df_para2, df_sum], axis=1)
    df_settings.columns, df_settings.index = np.arange(len(df_settings.columns)), ['target values', 't90 (s)']

    # stepwise combination of measurement and sensor settings
    df1 = pd.concat([df_meas, df_settings], axis=0)
    df1 = df1.loc[['run simulation', 'analyte', 'plateau time (s)', 'pKa', 'target values', 't90 (s)']]

    # add empty row to df1 for the ease of use in the report
    df1.loc[' ', :] = [' ']*df1.shape[1]

    # ------------------------------------------------------------
    # integration results
    int_out = pd.DataFrame.from_dict(dres)
    if int_out.empty:
        out = df1
    else:
        int_out.index = ['integration range (s)', 'target pH', 't90 pH', 'target Sum (mg/L)', 't90 meas sens', 
                         'integral target Sum', 'integral observed Sum', 'error']
        int_out.loc['peaks', :] = list([i-1 for i in dres.keys()])
        int_out.columns = np.arange(len(int_out.columns))
        int_out = int_out.loc[['peaks', 'integration range (s)', 't90 pH', 'target Sum (mg/L)', 't90 meas sens', 
                         'integral target Sum', 'integral observed Sum', 'error']]

        # ------------------------------------------------------------
        out = pd.concat([df1, int_out])

    return out


def save_report(para_meas, sensor_ph, sensor_para2, dres, df_res, total_lbl):
    # prep integral and settings for output
    out = save_integral(para_meas=para_meas, sensor_ph=sensor_ph, sensor_para2=sensor_para2, total_lbl=total_lbl,
                        dres=dres)

    # sort df_res index
    xold = list(df_res.index)
    xnew = ['header']
    [xnew.append(i) for i in xold]
    df_res.loc['header', :] = df_res.columns
    df_res = df_res.loc[xnew]

    return out, df_res


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
    # find the empty row indicating the break between the general dataframe and the integral
    df1_end = next(x[0] for x in enumerate(ls_lines) if len(x[1]) < 3)
    df_general = pd.DataFrame(ls_lines[:df1_end]).set_index(0)

    # get the information for pH
    col_ph = next(a[0] for a in enumerate(df_general.loc['analyte'].to_numpy()) if a[1] == 'pH') + 1
    df_ph = pd.DataFrame(df_general[col_ph].loc['analyte':])

    # get the information for sensor 2
    col_para2, para2 = col_ph + 1, df_general.loc['analyte'].to_numpy()[col_ph]
    df_para2 = pd.DataFrame(df_general[col_para2].loc['analyte':])

    # get the information for total parameter
    df_total = pd.DataFrame(df_general[col_para2 + 1].loc['analyte':])

    return df_general, df_ph, para2, df_para2, df_total
