__author__ = 'szieger'
__project__ = 'in silico study for sensor response'

import basics_sensorSignal_v4 as bs
import seaborn as sns
import re
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtGui import *
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QMessageBox, QGridLayout,
                             QLabel, QLineEdit, QGroupBox, QFileDialog, QFrame, QDialog, QTableWidget, QTableWidgetItem,
                             QWizard, QWizardPage, QComboBox)
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
import pyqtgraph.exporters
from scipy import integrate
import pandas as pd
import numpy as np
import os

# branch 1 is used to implement the reviewer's comments. 
# - allowing other concentrations -> drop down menu.
#

# .....................................................................................................................
# global parameter
dcolor = dict({'pH': '#1CC49A', 'sig pH': '#4B5258', 'para2': '#196E94', 'sig para2': '#314945', 'para3': '#DCA744',
               'sig para3': '#89621A', 'total para': '#A86349', 'background': '#383e42', 'font': 'white'})
ls = dict({'target': '-.', 'simulation': ':'})
save_type = ['jpg']
fs_font = 12
sns.set_context('paper', font_scale=1.)
pg.setConfigOption('background', dcolor['background'])
pg.setConfigOption('foreground', 'white')

# known acud/base systems
systems = dict({'TDS': dict({'base': ['HS-', 'HS', 'hs-', 'hs'], 'acid': ['H2S', 'h2s']}),
                'TAN': dict({'base': ['NH3', 'nh3'], 'acid': ['NH4+', 'NH4', 'nh4+', 'nh4']})})
ls_acidBase = [('TDS', 'HS-', 'H2S'), ('TAN', 'NH3', 'NH4+')]
acidBase = (None, None, None)
ls_unit = ['ppm', 'mg/L', 'mol/L']
unit = None


# fixed parameter depending on literature research
ph_res = 1e-3              # resolution of the pH sensor
ph_E0 = 0.22               # standard potential of the reference electrode of the pH meter; V
tsteps = 1e-3              # time steps for pH and measuring sensor (theory); seconds

# electrochemical sensor · 2
conc_calib = 1             # reaction quotient for calibrating the measuring sensor
sens_E = 0.09              # potential of the measuring sensor at certain concentration (conc_calib); V
sens_res = 1e-9            # resolution of the measuring sensor

# others
_integral_counter = 0       # helper scalar counting the integration
ls_lines = list()           # list to collect integration boundaries in pyqtgraph
dres = dict()               # simulation results
# get the local path for relative directories
loc_path = os.getcwd() 

wizard_page_index = {"IntroPage": 0, "SimPage": 1}


# .....................................................................................................................
class MagicWizard(QWizard):
    def __init__(self):
        super(MagicWizard, self).__init__()
        self.introPage = IntroPage()
        self.setPage(wizard_page_index["IntroPage"], self.introPage)
        self.simPage = SimPage()
        self.setPage(wizard_page_index["SimPage"], self.simPage)

        # set start page
        self.setStartId(wizard_page_index["IntroPage"])

        # GUI layout
        self.setWindowTitle("In silico simulation of pH-dependent systems")
        # define Wizard style and certain options
        self.setWizardStyle(QWizard.ClassicStyle)
        self.setOptions(QtWidgets.QWizard.NoCancelButtonOnLastPage |
                        QtWidgets.QWizard.IgnoreSubTitles)
        self.setStyleSheet("color: white; background-color: #383e42")

        # add a background image
        path = os.path.join(loc_path + r'/picture/icon.icns')
        pixmap = QtGui.QPixmap(path)
        pixmap = pixmap.scaled(200, 200, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        self.setPixmap(QWizard.BackgroundPixmap, pixmap)


class IntroPage(QWizardPage):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTitle("IntroPage | define parameter for initial plot")

        # create layout
        self.initUI()

        # connect checkbox and load file button with a function
        self.paraPair_label.currentIndexChanged.connect(self.onCurrentIndexChanged)
        self.onCurrentIndexChanged(self.paraPair_label.currentIndex())
        self.paraSum_conc_unit.currentIndexChanged.connect(self.onCurrentIndexChanged_Unit)
        self.onCurrentIndexChanged_Unit(self.paraSum_conc_unit.currentIndex()) 

        self.para1_conc_edit.editingFinished.connect(self.checkRange)
        self.para2_name_edit.editingFinished.connect(self.checkRange)
        self.para2_pKa_edit.editingFinished.connect(self.check_pka)
        self.load_button.clicked.connect(self.load_data)
        self.save_button.clicked.connect(self.save_path)

        # when all conditions are met, enable NEXT button:
        self.fname, self.fsave = QLineEdit(), QLineEdit()
        self.registerField("Data", self.fname)
        self.registerField("Storage path", self.fsave)
        self.registerField("parameter2*", self.para2_name_edit)
        self.registerField("conc. parameter1*", self.para1_conc_edit)
        self.registerField("conc. sum parameter*", self.paraSum_conc_edit)
        self.registerField("t90 parameter1*", self.para1_respT_edit)
        self.registerField("t90 parameter2*", self.para2_respT_edit)
        self.registerField("pKa parameter2*", self.para2_pKa_edit)
        self.registerField("plateau time", self.plateau_time_edit)
        self.registerField("unit*", self.paraSum_conc_unit, "currentText", self.paraSum_conc_unit.currentIndexChanged) 
        self.registerField("acid/base*", self.paraPair_label, "currentText", self.paraPair_label.currentIndexChanged) 

    def initUI(self):
        # path for measurement file (csv)
        self.load_button = QPushButton('Load meas. file', self)
        self.load_button.setFixedWidth(150)
        self.load_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.inputFileLineEdit = QLineEdit(self)
        self.inputFileLineEdit.setValidator(QtGui.QDoubleValidator())
        self.inputFileLineEdit.setFixedWidth(375) 
        self.inputFileLineEdit.setAlignment(Qt.AlignRight)
        self.inputFileLineEdit.setFont(QFont('Helvetica Neue', int(0.7 * fs_font)))

        # directory to store files
        self.save_button = QPushButton('Storage path', self)
        self.save_button.setFixedWidth(150)
        self.save_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.inputSaveLineEdit = QLineEdit(self)
        self.inputSaveLineEdit.setValidator(QtGui.QDoubleValidator())
        self.inputSaveLineEdit.setFixedWidth(375)
        self.inputSaveLineEdit.setAlignment(Qt.AlignRight)
        self.inputSaveLineEdit.setFont(QFont('Helvetica Neue', int(0.7 * fs_font)))

        # create label for parameters | name, concentration, response time, and pKa (if applicable)
        # while the first parameter will always be pH, the second varies
        para1_lbl, para1_name_lbl, para2_lbl, para2_name_lbl = QLabel(self), QLabel(self), QLabel(self), QLabel(self)
        para1_lbl.setText('Sensor · 1')
        para1_name_lbl.setText('Name')
        para2_lbl.setText('Sensor · 2')
        para2_name_lbl.setText('Name')
        para1_name, self.para2_name_edit = QLabel(self), QLineEdit(self)
        para1_name.setText('pH')
        para1_name.setAlignment(Qt.AlignRight)
        self.para2_name_edit.setAlignment(Qt.AlignRight)

        para1_conc_lbl, paraSum_conc_lbl = QLabel(self), QLabel(self)
        para1_conc_unit, self.paraSum_conc_unit = QLabel(self), QComboBox(self)
        para1_conc_lbl.setText('pH value(s)')
        paraSum_conc_lbl.setText('Concentration(s) sum parameter')
        para1_conc_unit.setText('')
        self.paraSum_conc_unit.setView(QtWidgets.QListView())
        [self.paraSum_conc_unit.addItem(k) for k in ['Select unit ...']+ls_unit]

        self.para1_conc_edit, self.paraSum_conc_edit = QLineEdit(self), QLineEdit(self)
        self.para1_conc_edit.setValidator(QRegExpValidator())
        self.para1_conc_edit.setAlignment(Qt.AlignRight)
        self.paraSum_conc_edit.setValidator(QRegExpValidator())
        self.paraSum_conc_edit.setAlignment(Qt.AlignRight)

        para1_respT_lbl, para2_respT_lbl = QLabel(self), QLabel(self)
        para1_respT_lbl.setText('Response time')
        para2_respT_lbl.setText('Response time')
        para1_respT_unit, para2_respT_unit = QLabel(self), QLabel(self)
        para1_respT_unit.setText('s')
        para2_respT_unit.setText('s')
        self.para1_respT_edit, self.para2_respT_edit = QLineEdit(self), QLineEdit(self)
        self.para1_respT_edit.setValidator(QRegExpValidator())
        self.para1_respT_edit.setAlignment(Qt.AlignRight)
        self.para2_respT_edit.setValidator(QRegExpValidator())
        self.para2_respT_edit.setAlignment(Qt.AlignRight)

        para2_pKa_lbl, self.para2_pKa_edit = QLabel(self), QLineEdit(self)
        para2_pKa_lbl.setText('pKa')
        self.para2_pKa_edit.setValidator(QRegExpValidator())
        self.para2_pKa_edit.setAlignment(Qt.AlignRight)

        # create plateau time
        plateau_time_lbl, plateau_time_unit = QLabel(self), QLabel(self)
        plateau_time_lbl.setText('Plateau time')
        plateau_time_unit.setText('s')
        self.plateau_time_edit = QLineEdit(self)
        self.plateau_time_edit.setValidator(QDoubleValidator())
        self.plateau_time_edit.setAlignment(Qt.AlignRight)
        self.plateau_time_edit.setText('75.')

        # additional information for corresponding acid/base pair and sum parameter
        global systems
        paraPair_label = QLabel(self)
        paraPair_label.setText('Corresponding acid/base pair')
        self.paraPair_label = QComboBox(self)
        self.paraPair_label.setView(QtWidgets.QListView())
        self.paraPair_label.addItem('Select acid/base pair ...')
        for k in systems.keys():
            for i in [(k, systems[k]['base'][0], systems[k]['acid'][0])]:        
                self.paraPair_label.addItem(', '.join([str(j.strip())[1:-1] for j in str(i)[1:-1].split(',')]))
        self.paraPair_label.addItem('Add other acid/base pairs ...')

        divider, divider1 = QLabel(self), QLabel(self)
        divider.setText('   ')
        divider1.setText('   ')

        # -----------------------------------------------------------------------
        # creating main window (GUI)
        w = QWidget()
        # create layout grid
        mlayout = QVBoxLayout(w)
        vbox_top, vbox_middle, vbox_bottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        mlayout.addLayout(vbox_top)
        mlayout.addLayout(vbox_middle)
        mlayout.addLayout(vbox_bottom)

        meas_settings = QGroupBox("Parameter Settings for Simulation")
        meas_settings.setStyleSheet("QGroupBox { font-weight: bold; } ")
        grid_load = QGridLayout()
        meas_settings.setFixedHeight(275)
        meas_settings.setFixedWidth(1000)
        meas_settings.setFont(QFont('Helvetica Neue', 12))
        vbox_top.addWidget(meas_settings)
        meas_settings.setLayout(grid_load)

        # include widgets in the layout
        # information for global system settings
        grid_load.addWidget(plateau_time_lbl, 0, 0)
        grid_load.addWidget(self.plateau_time_edit, 0, 2)
        grid_load.addWidget(plateau_time_unit, 0, 3)
        grid_load.addWidget(paraPair_label, 1, 0)
        grid_load.addWidget(self.paraPair_label, 1, 2)

        grid_load.addWidget(paraSum_conc_lbl, 2, 0)
        grid_load.addWidget(self.paraSum_conc_edit, 2, 2)
        grid_load.addWidget(self.paraSum_conc_unit, 2, 3)
        grid_load.addWidget(divider, 3, 2)

        # information for parameter 1
        grid_load.addWidget(para1_lbl, 4, 0)
        grid_load.addWidget(para1_name_lbl, 4, 1)
        grid_load.addWidget(para1_name, 4, 2)
        grid_load.addWidget(para1_conc_lbl, 5, 1)
        grid_load.addWidget(self.para1_conc_edit, 5, 2)
        grid_load.addWidget(para1_conc_unit, 5, 3)
        grid_load.addWidget(para1_respT_lbl, 6, 1)
        grid_load.addWidget(self.para1_respT_edit, 6, 2)
        grid_load.addWidget(para1_respT_unit, 6, 3)
        grid_load.addWidget(divider1, 7, 2)

        # information for parameter 2
        grid_load.addWidget(para2_lbl, 8, 0)
        grid_load.addWidget(para2_name_lbl, 8, 1)
        grid_load.addWidget(self.para2_name_edit, 8, 2)
        grid_load.addWidget(para2_respT_lbl, 9, 1)
        grid_load.addWidget(self.para2_respT_edit, 9, 2)
        grid_load.addWidget(para2_respT_unit, 9, 3)
        grid_load.addWidget(para2_pKa_lbl, 10, 1)
        grid_load.addWidget(self.para2_pKa_edit, 10, 2)

        divide_group = QGroupBox()
        divide_group.setStyleSheet("QGroupBox { border: 0px solid #383e42;}")
        grid_divide = QGridLayout()
        divide_group.setMinimumHeight(75)
        vbox_middle.addWidget(divide_group)
        divide_group.setLayout(grid_divide)

        meas_file = QGroupBox("Define Directories")
        meas_file.setStyleSheet("QGroupBox { font-weight: bold; } ")
        grid_file = QGridLayout()
        meas_file.setFixedHeight(100)
        meas_file.setFixedWidth(1000)
        meas_file.setFont(QFont('Helvetica Neue', 12))
        vbox_bottom.addWidget(meas_file)
        meas_file.setLayout(grid_file)

        # include widgets in the layout
        grid_file.addWidget(self.load_button, 0, 0)
        grid_file.addWidget(self.inputFileLineEdit, 0, 1)
        grid_file.addWidget(self.save_button, 1, 0)
        grid_file.addWidget(self.inputSaveLineEdit, 1, 1)

        self.setLayout(mlayout)
        # margin order left, top, right, bottom
        self.layout().setContentsMargins(50, 50, 170, 50)
        self.layout().addStretch()

    def check_measSensor(self):
        global ls_acidBase
        ls_check = [True for p in ls_acidBase if self.para2_name_edit.text() in p]
        if len(ls_check) == 0:
            # raise Error that the measurement sensor does not fit to the defined system
            msgBox = error_message(text="The defined measurment sensor does not fit the  defined"
                                   "acid/base pair", type='Error')
            msgBox.exec_()

    def check_pka(self):
        if re.match(r'^-?\d+(?:\.\d+)$', self.para2_pKa_edit.text()) is None:
            if re.match(r'^-?\d+(?:)$', self.para2_pKa_edit.text()):
               pka = float(self.para2_pKa_edit.text()) 
            else:
                pka = -2.
        else:
            pka = float(self.para2_pKa_edit.text().replace(',', '.'))

        if pka < 0 or pka > 20.:
            msgBox = error_message(
                text="Double check pKa value.", type='Error')
            returnValue = msgBox.exec()
            if returnValue == QMessageBox.Ok:
                pass

    def checkRange(self):
        if self.para2_name_edit.text():
            self.check_measSensor()

        if self.para1_conc_edit.text() and self.para2_name_edit.text():
            # check whether we work with H2S sensor / TDS project
            condA = ('H2S' in self.para2_name_edit.text()
                     or 'HS-' in self.para2_name_edit.text())
            if condA or 'HS' in self.para2_name_edit.text():
                for ph in [float(i.strip()) for i in self.para1_conc_edit.text().split(',')]:
                    if ph > 9:
                        text = "Given pH value out of simulation range.  The H2S sensor simulation is only valid up " \
                               "to pH 9."
                        msgBox = error_message(text=text, type='Error')
                        returnValue = msgBox.exec()
                        if returnValue == QMessageBox.Ok:
                            pass

    def load_data(self):
        # opens a dialog window in the current path
        fname, filter = QFileDialog.getOpenFileName(self, "Select specific txt file for temperature compensation",
                                                    "", "Text files (*.txt *.csv *xls)")

        if fname:
            self.inputFileLineEdit.setText(fname)
            df_general, df_ph, para2, df_para, df_total = bs.load_data(fname)

            # set parameter to run the simulation
            self.plateau_time_edit.setText(
                df_general.loc['plateau time (s)'].to_numpy()[0])

            # pH sensor
            self.para1_respT_edit.setText(df_ph.loc['t90 (s)', df_ph.columns[0]])
            ph_target = df_ph.loc['target values', df_ph.columns[0]]
            self.para1_conc_edit.setText(ph_target)

            # info for sensor 2 and total parameter concentration
            self.para2_name_edit.setText(para2)
            self.para2_respT_edit.setText(df_para.loc['t90 (s)', df_para.columns[0]])
            self.para2_pKa_edit.setText(df_para.loc['pKa', df_para.columns[0]])
            self.paraSum_conc_edit.setText(
                df_total.loc['target values'].to_numpy()[0])

    def save_path(self):
        fsave = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Select Folder')
        self.fsave.setText(fsave)
        if fsave:
            self.inputSaveLineEdit.setText(fsave)
            self.fsave.setText(fsave)

    def onCurrentIndexChanged_Unit(self, index):
        global unit
        unit = ls_unit[index-1]

    def onCurrentIndexChanged(self, index):
        global acidBase, ls_acidBase
        # if the last item was selected, then open another window to edit acid/base pair. 
        # counting starts with 0, hence lenght of combobox-1 
        if index == self.paraPair_label.count()-1:
            # open a pop up window with options to select what shall be saved
            global wEdit
            wEdit = AcidBasePair(index, self.paraPair_label)

            if wEdit.isVisible() is False:
                wEdit.show() 
        else:
            acidBase = ls_acidBase[index-1]
        
    def nextId(self) -> int:
        return wizard_page_index["SimPage"]


class AcidBasePair(QDialog):
    def __init__(self, index, box):
        super().__init__()
        self.index, self.box, self.indDone = index, box, 0
        # set title and initiate general layout
        self.initUI()

        # -----------------------------------------------------------
        # add action function when done with editing
        self.paraSum_label_edit.editingFinished.connect(self.updateComboBox)
        self.paraAcid_label_edit.editingFinished.connect(self.updateComboBox)
        self.paraBase_label_edit.editingFinished.connect(self.updateComboBox)
        self.close_button.clicked.connect(self.close_window)
        
    def initUI(self):
        self.setWindowTitle("Add Corresp. Acid/Base Pair")
        # define geometry of integration window x, y, width, height
        self.setGeometry(50, 100, 300, 300)

        # DONE button
        self.close_button = QPushButton('DONE', self)
        self.close_button.setFixedWidth(100)
        self.close_button.setEnabled(False)

        # acid, base and sum parameter text window
        paraSum_label = QLabel(self)
        paraSum_label.setText('Sum Parameter')

        self.paraSum_label_edit = QLineEdit(self)
        self.paraSum_label_edit.setAlignment(Qt.AlignRight)

        paraAcid_label = QLabel(self)
        paraAcid_label.setText('Corresponding Acid')

        self.paraAcid_label_edit = QLineEdit(self)
        self.paraAcid_label_edit.setAlignment(Qt.AlignRight)

        paraBase_label = QLabel(self)
        paraBase_label.setText('Corresponding Base')
        self.paraBase_label_edit = QLineEdit(self)
        self.paraBase_label_edit.setAlignment(Qt.AlignRight)

        # creating window layout
        mlayoutP = QVBoxLayout()
        vbox_center = QHBoxLayout()
        mlayoutP.addLayout(vbox_center)

        # panel for integration results
        pair_grp = QGroupBox("Define your System")
        grid_pair = QGridLayout()
        pair_grp.setFont(QFont('Helvetica Neue', 10))
        vbox_center.addWidget(pair_grp)
        pair_grp.setLayout(grid_pair)

        grid_pair.addWidget(paraSum_label, 0, 0)
        grid_pair.addWidget(self.paraSum_label_edit, 0, 1)
        grid_pair.addWidget(paraAcid_label, 1, 0)
        grid_pair.addWidget(self.paraAcid_label_edit, 1, 1)
        grid_pair.addWidget(paraBase_label, 2, 0)
        grid_pair.addWidget(self.paraBase_label_edit, 2, 1)
        grid_pair.addWidget(self.close_button, 4, 1) 

         # add everything to the window layout
        self.setLayout(mlayoutP)

    def updateComboBox(self):
        if self.paraSum_label_edit.text() and self.paraAcid_label_edit.text() and self.paraBase_label_edit.text():
            # collect the parameters and update the global dictonary of provided systems
            global systems, ls_acidBase
            if self.paraSum_label_edit.text() in systems.keys():
                if self.paraAcid_label_edit.text() not in systems[self.paraSum_label_edit.text()]['acid']:
                    systems[self.paraSum_label_edit.text()]['acid'].append(self.paraAcid_label_edit.text())
                if self.paraBase_label_edit.text() not in systems[self.paraSum_label_edit.text()]['base']:
                    systems[self.paraSum_label_edit.text()]['base'].append(self.paraBase_label_edit.text())               
            else:
                systems[self.paraSum_label_edit.text()] = dict({'base': [self.paraBase_label_edit.text()],
                                                                'acid': [self.paraAcid_label_edit.text()]})
            
            ls_acidBase.append((self.paraSum_label_edit.text(), self.paraAcid_label_edit.text(), self.paraBase_label_edit.text()))
            # allow closing the window when all parameters are provided
            self.close_button.setEnabled(True)

    def close_window(self):
        global ls_acidBase, acidBase
        acidBase = ls_acidBase[self.index-1]
        initializeComboBox(box=self.box)
        self.hide()


def initializeComboBox(box):
    global ls_acidBase
    # clear comboBox
    box.clear()
    # add all items to comboBox (without duplicates)
    box.addItem('Select acid/base pair ...')
    for k in list(dict.fromkeys(ls_acidBase)):
        box.addItem(', '.join([str(j.strip())[1:-1] for j in str(k)[1:-1].split(',')]))
    box.addItem('Add other acid/base pairs ...') 
    box.setCurrentIndex(len(ls_acidBase))

    # after updating the comboBox and selecting the defined acid/base pair, close window
    global wEdit
    wEdit.close()


# .....................................................................................................................
class SimPage(QWizardPage):
    def __init__(self, parent=None):
        super(SimPage, self).__init__(parent)
        self.setTitle("in silico Simulation")

        # define certain parameter and gather from intro page
        global dres
        self.df_res = None

        # general layout
        self.initUI()

        # -----------------------------------------------------------
        # for all parameters - connect LineEdit with function
        self.initializePage()
        self.ph_t90_edit.returnPressed.connect(self.print_ph_t90)
        self.para2_t90_edit.returnPressed.connect(self.print_para2_t90)

        # connect buttons in navigation manager with functions
        self.plot_button.clicked.connect(self.check_parameter)
        self.int_button.clicked.connect(self.calc_integral)
        self.clearP_button.clicked.connect(self.clear_parameters)
        self.clearF_button.clicked.connect(self.clear_phsim)
        self.clearF_button.clicked.connect(self.clear_para2timedrive)
        self.clearF_button.clicked.connect(self.clear_sumtimdrive)
        self.save_button.clicked.connect(self.save)
        self.saveR_button.clicked.connect(self.save_report)

    def initUI(self):
        # general parameters
        tsteady_lbl, tsteady_unit = QLabel(self), QLabel(self)
        tsteady_lbl.setText('Plateau time')
        tsteady_unit.setText('s')
        self.tsteady_edit = QLineEdit(self)
        self.tsteady_edit.setValidator(QDoubleValidator())
        self.tsteady_edit.setAlignment(Qt.AlignRight)

        # pH sensor settings
        pH_lbl = QLabel(self)
        pH_lbl.setText('pH value(s)')
        self.ph_edit = QLineEdit(self)
        self.ph_edit.setValidator(QRegExpValidator())
        self.ph_edit.setAlignment(Qt.AlignRight)

        ph_t90_label, ph_t90_unit = QLabel(self), QLabel(self)
        ph_t90_label.setText('Response time')
        ph_t90_unit.setText('s')
        self.ph_t90_edit = QLineEdit(self)
        self.ph_t90_edit.setValidator(QDoubleValidator())
        self.ph_t90_edit.setAlignment(Qt.AlignRight)

        # second sensor settings
        self.para2_lbl, self.para2_unit = QLabel(self), QLabel(self)
        self.paraSum_conc_edit = QLineEdit(self)
        self.paraSum_conc_edit.setFixedWidth(100)
        self.paraSum_conc_edit.setValidator(QRegExpValidator())
        self.paraSum_conc_edit.setAlignment(Qt.AlignRight)

        para2_t90_lbl, para2_t90_unit = QLabel(self), QLabel(self)
        para2_t90_lbl.setText('Response time')
        para2_t90_unit.setText('s')
        self.para2_t90_edit = QLineEdit(self)
        self.para2_t90_edit.setFixedWidth(100)
        self.para2_t90_edit.setValidator(QDoubleValidator())
        self.para2_t90_edit.setAlignment(Qt.AlignRight)

        para2_pka_lbl = QLabel(self)
        para2_pka_lbl.setText('pKa')
        self.para2_pka_edit = QLineEdit(self)
        self.para2_pka_edit.setFixedWidth(100)
        self.para2_pka_edit.setValidator(QDoubleValidator())
        self.para2_pka_edit.setAlignment(Qt.AlignRight)

        # general navigation
        self.plot_button = QPushButton('Plot', self)
        self.plot_button.setFixedWidth(100)
        self.plot_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.int_button = QPushButton('Integral', self)
        self.int_button.setFixedWidth(100)
        self.int_button.setEnabled(False)
        self.int_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.clearP_button = QPushButton('Clear parameter', self)
        self.clearP_button.setFixedWidth(150)
        self.clearP_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.clearF_button = QPushButton('Clear plots', self)
        self.clearF_button.setFixedWidth(100)
        self.clearF_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.save_button = QPushButton('Save all', self)
        self.save_button.setFixedWidth(100)
        self.save_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.saveR_button = QPushButton('Save report', self)
        self.saveR_button.setFixedWidth(100)
        self.saveR_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
  
        # ---------------------------------------------------------------------------------------
        w = QWidget()
        mlayout = QVBoxLayout(w)

        # 1st layer: box with vertical alignment to split in top and bottom
        vbox_top, vbox_middle, vbox_bottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        mlayout.addLayout(vbox_top)
        mlayout.addLayout(vbox_middle)
        mlayout.addLayout(vbox_bottom)

        # 2nd layer: box with vertical alignment to split in left and right (left: parameters, m1: line, right: plots)
        hbox_left, hbox_m1, hbox_right = QVBoxLayout(), QVBoxLayout(), QVBoxLayout()
        vbox_top.addLayout(hbox_left)
        vbox_top.addLayout(hbox_m1)
        vbox_top.addLayout(hbox_right)

        # 3rd layer: left box with parameter settings (top, middle, bottom)
        hbox_ltop, hbox_lmiddle, hbox_lbottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        hbox_ltop.setContentsMargins(5, 10, 10, 5)
        hbox_left.addLayout(hbox_ltop)
        hbox_left.addLayout(hbox_lmiddle)
        hbox_left.addLayout(hbox_lbottom)

        # 4th layer: right box split horizontally (top: individual sensors bottom: Sum)
        hbox_rtop, hbox_rbottom = QHBoxLayout(), QHBoxLayout()
        hbox_right.addLayout(hbox_rtop)
        hbox_right.addLayout(hbox_rbottom)

        # 5th layer: right top box split vertically (left: pH right: Base/Acid)
        hbox_tright, hbox_tleft = QVBoxLayout(), QVBoxLayout()
        hbox_rtop.addLayout(hbox_tright)
        hbox_rtop.addLayout(hbox_tleft)

        # draw additional "line" to separate parameters from plots and to separate navigation from rest
        vline = QFrame()
        vline.setFrameShape(QFrame.VLine | QFrame.Raised)
        vline.setLineWidth(2)
        hbox_m1.addWidget(vline)

        hline = QFrame()
        hline.setFrameShape(QFrame.HLine | QFrame.Raised)
        hline.setLineWidth(2)
        vbox_middle.addWidget(hline)
        
        # left side of main window (-> data treatment)
        vbox_top.setContentsMargins(5, 10, 10, 5)

        # ---------------------------------------------------------------------------------------------------------
        # draw additional "line" to separate load from plot and plot from save
        vline1 = QFrame()
        vline1.setFrameShape(QFrame.VLine | QFrame.Raised)
        vline1.setLineWidth(2)
        vline2 = QFrame()
        vline2.setFrameShape(QFrame.VLine | QFrame.Raised)
        vline2.setLineWidth(2)
        vline3 = QFrame()
        vline3.setFrameShape(QFrame.VLine | QFrame.Raised)
        vline3.setLineWidth(2)

        # -------------------------------------------------------------------------------------------
        # GroupBoxes to structure the layout
        navigation_group = QGroupBox("Navigation Tool Bar")
        navigation_group.setStyleSheet("QGroupBox { font - weight: bold; border: 0px solid gray;}")
        grid_load = QGridLayout()
        navigation_group.setFixedHeight(70)
        vbox_bottom.addWidget(navigation_group)
        navigation_group.setLayout(grid_load)

        grid_load.addWidget(self.plot_button, 0, 1)
        grid_load.addWidget(self.int_button, 0, 2)
        grid_load.addWidget(vline2, 0, 3)
        grid_load.addWidget(self.clearP_button, 0, 4)
        grid_load.addWidget(self.clearF_button, 0, 5)
        grid_load.addWidget(vline3, 0, 6)
        grid_load.addWidget(self.save_button, 0, 7)
        grid_load.addWidget(self.saveR_button, 0, 8)

        # create GroupBox to structure the layout
        general_group = QGroupBox("General Settings")
        general_group.setStyleSheet("QGroupBox { font - weight: bold; border: 0px solid gray;}")
        general_group.setMinimumHeight(300)
        grid_load = QGridLayout()
        grid_load.setSpacing(1)
        grid_load.setVerticalSpacing(1)        

        # add GroupBox to layout and load buttons in GroupBox
        hbox_ltop.addWidget(general_group)
        general_group.setLayout(grid_load)
        grid_load.addWidget(tsteady_lbl, 0, 0)
        grid_load.addWidget(self.tsteady_edit, 0, 1)
        grid_load.addWidget(tsteady_unit, 0, 2)

        general_group.setContentsMargins(1, 2, 10, 1)
        hbox_ltop.addSpacing(10)

        # -----------------------
        # pH Sensor Settings
        phsens_group = QGroupBox("pH Sensor Settings")
        phsens_group.setStyleSheet("QGroupBox { font - weight: bold; border: 0px solid gray;}")
        phsens_group.setFixedHeight(110)
        grid_load = QGridLayout()
        grid_load.setSpacing(5)
        grid_load.setVerticalSpacing(1)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_lmiddle.addWidget(phsens_group)
        phsens_group.setLayout(grid_load)

        grid_load.addWidget(pH_lbl, 0, 0)
        grid_load.addWidget(self.ph_edit, 0, 1)
        grid_load.addWidget(ph_t90_label, 1, 0)
        grid_load.addWidget(self.ph_t90_edit, 1, 1)
        grid_load.addWidget(ph_t90_unit, 1, 2)

        phsens_group.setContentsMargins(1, 15, 15, 1)
        hbox_lmiddle.addSpacing(10)

        # -----------------------
        # Settings of the other sensor
        self.para23sens_group = QGroupBox()
        self.para23sens_group.setStyleSheet("QGroupBox { font - weight: bold; border: 0px solid gray;}")
        self.para23sens_group.setFixedHeight(150)
        grid_load = QGridLayout()
        grid_load.setSpacing(5)
        grid_load.setVerticalSpacing(1)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_lbottom.addWidget(self.para23sens_group)
        self.para23sens_group.setLayout(grid_load)

        grid_load.addWidget(self.para2_lbl, 0, 0)
        grid_load.addWidget(self.paraSum_conc_edit, 0, 1)
        grid_load.addWidget(self.para2_unit, 0, 2)
        grid_load.addWidget(para2_t90_lbl, 1, 0)
        grid_load.addWidget(self.para2_t90_edit, 1, 1)
        grid_load.addWidget(para2_t90_unit, 1, 2)
        grid_load.addWidget(para2_pka_lbl, 2, 0)
        grid_load.addWidget(self.para2_pka_edit, 2, 1)

        self.para23sens_group.setContentsMargins(1, 15, 15, 1)
        hbox_lbottom.addSpacing(10)

        # --------------------------------------------------------------------------------------------------------------
        # pH Simulation
        self.fig_phsim = pg.PlotWidget()
        self.fig_phsim.setLabel('left', '<font>pH value</font>', color=dcolor['font'])
        self.fig_phsim.setLabel('bottom', '<font>Time</font>', units='<font>s</font>', color=dcolor['font'])

        # create GroupBox to structure the layout
        phsim_group = QGroupBox("pH Sensor Simulation")
        phsim_group.setMinimumWidth(450)
        phsim_group.setMinimumHeight(250)
        grid_phsim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_tright.addWidget(phsim_group)
        phsim_group.setLayout(grid_phsim)
        grid_phsim.addWidget(self.fig_phsim)

        # parameter 2-3 simulation
        self.fig_para2sim = pg.PlotWidget()

        # Axis and ViewBoxes
        self.ax_para2sim = self.fig_para2sim.plotItem
        self.ax1_para2sim = pg.ViewBox()
        self.fig_para2sim.scene().addItem(self.ax1_para2sim)
        self.ax_para2sim.getAxis('right').linkToView(self.ax1_para2sim)
        self.ax1_para2sim.setXLink(self.ax_para2sim)

        def updateViews():
            # view has resized; update auxiliary views to match
            self.ax1_para2sim.setGeometry(
                self.ax_para2sim.vb.sceneBoundingRect())
        updateViews()
        self.ax_para2sim.vb.sigResized.connect(updateViews)

        # ------------------------------------
        self.para2sim_group = QGroupBox("Single Parameter Simulation")
        self.para2sim_group.setMinimumWidth(450)
        self.para2sim_group.setMinimumHeight(250)
        grid_para2sim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_tleft.addWidget(self.para2sim_group)
        self.para2sim_group.setLayout(grid_para2sim)
        grid_para2sim.addWidget(self.fig_para2sim)

        # total parameter simulation
        self.fig_totalsim = pg.GraphicsLayoutWidget(show=True)
        self.ax_totalsim = self.fig_totalsim.addPlot(row=1, col=0)
        self.ax_totalsim.setAutoVisible(y=True)
        label = pg.LabelItem(justify='right')
        self.fig_totalsim.addItem(label)

        self.totalsim_group = QGroupBox()
        self.totalsim_group.setMinimumWidth(900)
        self.totalsim_group.setMinimumHeight(300)
        grid_totalsim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_rbottom.addWidget(self.totalsim_group)
        self.totalsim_group.setLayout(grid_totalsim)
        grid_totalsim.addWidget(self.fig_totalsim)

        # -------------------------------------------------------------------------------------------------------------
        self.setLayout(mlayout)

    # ---------------------------------------------------
    def print_ph_t90(self):
        print('pH sensor response: ', self.para1_t90_edit.text(), 's')

        # in case set response time equals 0, return warning
        if float(self.para1_t90_edit.text()) == 0:
            msgBox = error_message(text="Please enter a valid sensor response time for the pH sensor, notably a "
                                        "number > 0")
            returnValue = msgBox.exec()
            if returnValue == QMessageBox.Ok:
                pass

    def print_para2_t90(self):
        print('{} sensor response: '.format(self.para2_name_edit.text()),
              self.para2_t90_edit.text(), 's')
        # in case set response time equals 0, return warning
        if float(self.para2_t90_edit.text()) == 0:
            msgBox = error_message(text="Please enter a valid sensor response time for the total sensor, notably a "
                                        "number > 0")
            returnValue = msgBox.exec()
            if returnValue == QMessageBox.Ok:
                pass

    def initializePage(self):
        self.save_path = self.field("Storage path")
        self.para1_conc, self.para2_conc = self.field("conc. parameter1"), self.field("conc. sum parameter")
        if isinstance(self.field("t90 parameter1"), type(None)):
            self.para1_respT = self.field("t90 parameter1")
        else:
            self.para1_respT = str(self.field("t90 parameter1").replace(',', '.'))
        if isinstance(self.field("t90 parameter2"), type(None)):
            self.para2_respT = self.field("t90 parameter2")
        else:
            self.para2_respT = str(self.field("t90 parameter2").replace(',', '.'))
        if isinstance(self.field("pKa parameter2"), type(None)):
            self.para2_pka = self.field("pKa parameter2")
        else:
            self.para2_pka = str(self.field("pKa parameter2").replace(',', '.'))
        if isinstance(self.field("plateau time"), type(None)):
            self.plateau_time = self.field("plateau time")
        else:
            self.plateau_time = float(self.field("plateau time").replace(',', '.'))
        self.para2_name = self.field("parameter2")
        
        # identify the system according to given parameter
        # para2 will be the measured, para3 the calculated parameter
        global systems, acidBase

        if self.para2_name:
            # H2S/HS- known system
            if self.para2_name in systems['TDS']['base']:
                self.para2_grp, self.para3_grp, self.total_para_grp = 'HS-', 'H2S', 'TDS'
                self.para2, self.para3, self.total_para = 'HS$^-$', 'H$_2$S', '$TDS$'
                self.analyte = 'base'
            elif self.para2_name in systems['TDS']['acid']:
                self.para2_grp, self.para3_grp, self.total_para_grp = 'H2S',  'HS-', 'TDS'
                self.para2, self.para3, self.total_para = 'H$_2$S', 'HS$^-$', '$TDS$'
                self.analyte = 'acid'

            # NH3/NH4+ known system
            if self.para2_name in systems['TAN']['base']:
                self.para2_grp, self.para3_grp, self.total_para_grp = 'NH3', 'NH4+', 'TAN'
                self.para2, self.para3, self.total_para = 'NH$_3$', 'NH$_4^+$', '$TAN$'
                self.analyte = 'base'
            elif self.para2_name in systems['TAN']['acid']:
                self.para2_grp, self.para3_grp, self.total_para_grp = 'NH4+',  'NH3', 'TAN'
                self.para2, self.para3, self.total_para = 'NH$_4^+$', 'NH$_3$', '$TAN$'
                self.analyte = 'acid'

            # if unknown use (sum, acid, base) by getting the next keyword from system-dictionary
            ls_keyNew = list()
            for k in systems.keys():
                if k != 'TAN' and k != 'TDS':
                    ls_keyNew.append(k)

            if len(ls_keyNew) > 0:
                self.total_para, self.total_para_grp = acidBase[0], acidBase[0]
                # define whether it is acid or base: 
                if self.para2_name in systems[ls_keyNew[0]]['acid']:
                    self.analyte = 'acid'
                    self.para2, self.para2_grp = systems[ls_keyNew[0]]['acid'], systems[ls_keyNew[0]]['acid']
                    self.para3, self.para3_grp = systems[ls_keyNew[0]]['base'], systems[ls_keyNew[0]]['base']
                elif self.para2_name in systems[ls_keyNew[0]]['base']:
                    self.analyte = 'base'
                    self.para2, self.para2_grp = systems[ls_keyNew[0]]['base'], systems[ls_keyNew[0]]['base']
                    self.para3, self.para3_grp = systems[ls_keyNew[0]]['acid'], systems[ls_keyNew[0]]['acid']

                if isinstance(self.para2, list):
                    self.para2, self.para2_grp = self.para2[0], self.para2_grp[0]
                if isinstance(self.para3, list):
                    self.para3, self.para3_grp = self.para3[0], self.para3_grp[0]
                    
            # if unknown use "base / acid / sum" (the first one is the measured, the second one the calculated)
            if '/' in self.para2_name:
                system = self.para2_name.split('/')
                self.para2_grp, self.para3_grp, self.total_para_grp = system[0], system[1], system[-1]
                self.para2, self.para3, self.total_para = system[0], system[1], system[-1]
                self.analyte = system[0]
        else:
            self.para2_grp, self.para3_grp, self.total_para_grp = 'Base', 'Acid', 'Sum Parmeter'
            self.para2, self.para3, self.total_para = 'base', 'acid', 'sum parmeter'
            self.analyte = 'sum parmeter'
        self.system_select = (self.para2_grp, self.para3_grp, self.total_para_grp)

        # fill in parameters
        global unit
        self.para2_unit.setText(unit)
        
        # Initialize units for axes
        self.ax_para2sim.setLabel('bottom', '<font>Time</font>', units='<font>s</font>', color=dcolor['font'])
        self.ax_totalsim.setLabel('bottom', '<font>Time</font>', units='<font>s</font>', color=dcolor['font'])

        self.ax_para2sim.setLabel('left', '<font>Measured analyte</font>', units='<font>{}</font>'.format(unit), color=dcolor['font'])
        self.ax_para2sim.setLabel('right', '<font>Determined analyte</font>', units='<font>{}</font>'.format(unit), color=dcolor['font'])
        self.ax_totalsim.setLabel('left', '<font>Sum concentration</font>', units='<font>{}</font>'.format(unit), color=dcolor['font'])

        if self.plateau_time:
            self.tsteady_edit.setText(str(self.plateau_time))
        self.ph_edit.setText(self.para1_conc)
        self.ph_t90_edit.setText(self.para1_respT)

        self.paraSum_conc_edit.setText(self.para2_conc)
        self.para2_lbl.setText('c({})'.format(self.total_para_grp))
        self.para2_t90_edit.setText(self.para2_respT)
        self.para2_pka_edit.setText(self.para2_pka)

        self.totalsim_group.setTitle("{} Simulation".format(acidBase[0]))
        self.para2sim_group.setTitle("{}/{} Simulation".format(acidBase[1], acidBase[2]))
        self.para23sens_group.setTitle("{}/{} Sensor Settings".format(acidBase[1], acidBase[2]))

        return

    # ---------------------------------------------------
    def clear_parameters(self):
        # re-write default parameters
        self.tsteady_edit.setText('75.')
        self.ph_edit.setText('8.4, 10.9')
        self.ph_t90_edit.setText('30.')
        self.paraSum_conc_edit.setText('100.0')
        self.para2_t90_edit.setText('60.')
        self.para2_pka_edit.setText('9.25')

    def clear_phsim(self):
        self.fig_phsim.clear()
        self.fig_phsim.setLabel('bottom', '<font>Time</font>', units='<font>s</font>', color=dcolor['font'])
        self.fig_phsim.setLabel('left', 'pH value', color=dcolor['font'])

    def clear_para2timedrive(self):
        global unit
        self.fig_para2sim.clear()
        self.ax_para2sim = self.fig_para2sim.plotItem
        self.ax1_para2sim.clear()
        self.ax_para2sim.getAxis('right').linkToView(self.ax1_para2sim)
        self.ax1_para2sim.setXLink(self.ax_para2sim)

        # Axes layout
        self.ax_para2sim.setLabel('bottom', '<font>Time</font>', units='<font>s</font>', color=dcolor['font'])
        self.ax_para2sim.setLabel('left', '<font>Measured analyte</font>', units='<font>{}</font>'.format(unit),
                                  color=dcolor['font'])
        self.ax_para2sim.setLabel('right', '<font>Determined analyte</font>', units='<font>{}</font>'.format(unit),
                                  color=dcolor['font'])

    def clear_sumtimdrive(self):
        global unit
        self.fig_totalsim.clear()
        self.ax_totalsim = self.fig_totalsim.addPlot(row=1, col=0)
        self.ax_totalsim.setAutoVisible(y=True)
        label = pg.LabelItem(justify='right')
        self.fig_totalsim.addItem(label)

        # axis label
        self.ax_totalsim.setLabel('bottom', '<font>Time</font>', units='<font>s</font>', color=dcolor['font'])
        if 'TAN' in self.total_para:
            self.ax_totalsim.setLabel('left', '<math>TAN</math>', units='<font>{}</font>'.format(unit), color=dcolor['font'])
        else:
            self.ax_totalsim.setLabel('left', '<font>Sum concentration</font>', units='<font>{}</font>'.format(unit),
                                      color=dcolor['font'])

    def align_concentrations(self, ls_ph, ls_cSum_gmL):
        if len(ls_cSum_gmL) == len(ls_ph):
            pass
        elif len(ls_cSum_gmL) == 1 and len(ls_ph) > 1:
            ls_cSum_gmL = ls_cSum_gmL * len(ls_ph)
        elif len(ls_cSum_gmL) > 1 and len(ls_ph) == 1:
            ls_ph = ls_ph * len(ls_cSum_gmL)
        else:
            if max(len(ls_ph), len(ls_cSum_gmL)) % min(len(ls_ph), len(ls_cSum_gmL)) == 0:
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

    def parameter_prep(self, t_plateau, ls_ph, ls_total, pKa):
        # target fluctuation of analytes
        target_ph = bs.target_fluctuation(ls_conc=ls_ph, tstart=0, tstop=t_plateau * 2, nP=1, analyte='pH')
        target_sum = bs.target_fluctuation(ls_conc=ls_total, tstart=0, tstop=t_plateau * 2, nP=1, analyte='Sum')

        # target concentration para2/para3 based on target pH (para2... base, para3... acid)
        para2 = pd.DataFrame(bs.henderson_Sum4Base(pKa=pKa, pH=target_ph['signal pH'].to_numpy(),
                                                   c_sum=target_sum['signal Sum'].to_numpy()), index=target_ph.index,
                             columns=['signal Base'])
        para3 = pd.DataFrame(bs.henderson_Sum4Acid(pKa=pKa, pH=target_ph['signal pH'].to_numpy(),
                                                   c_sum=target_sum['signal Sum'].to_numpy()), index=target_ph.index,
                             columns=['signal Acid'])

        # calculate target concentration of individual NHx parameters
        df_res = pd.concat([target_ph, target_sum, para2, para3], axis=1)
        return df_res

    # ---------------------------------------------------
    def save(self):
        # opens window in current path. User input to define file name
        fname_save = QFileDialog.getSaveFileName(self, 'Save File')[0]

        if fname_save:
            if '.txt' in fname_save or '.csv' in fname_save:
                pass
            else:
                fname_save = fname_save + '.txt'

            # check whether we have data to save
            try:
                if self.df_res is None:
                    msgBox = error_message(
                        text="Simulate before saving. Simulation might have been unsuccessful.")
                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # save all reports
                    self.save2txt(fname_save)

                    # save figures in separate folder
                    for f in self.dic_figures.keys():
                        # create figure name and path
                        figure_name = fname_save.split('.')[0] + '_Graph-' + '-'.join(f.split(' ')) + '.'

                        for t in save_type:
                            # create an exporter instance, as an argument give it the item you wish to export
                            exporter = pg.exporters.ImageExporter(self.dic_figures[f].getPlotItem())
                            # set export parameters if needed
                            # (note this also affects height parameter)
                            exporter.parameters()['width'] = 1000

                            # save to file
                            if t == 'jpg':
                                myQImage = exporter.export(toBytes=True)
                                myQImage.save(figure_name + t, "JPG", 100)
                            else:
                                exporter.export(figure_name + t)

            except NameError:
                msgBox = error_message(
                    text="Simulate before saving", type='Error')
                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    def save_report(self):
        # opens window in current path. User input to define file name
        fname_save = QFileDialog.getSaveFileName(self, 'Save File')[0]

        if fname_save:
            if '.txt' in fname_save or '.csv' in fname_save:
                pass
            else:
                fname_save = fname_save + '.txt'

            # check whether we have data to save
            try:
                if self.df_res is None:
                    msgBox = error_message(
                        text="Simulate before saving", type='Error')
                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # save the reports
                    self.save2txt(fname_save)

            except NameError:
                msgBox = error_message(
                    text="Simulate before saving", type='Error')
                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    def save2txt(self, fname_save):
        global dres
        # if integral available, convert to DataFrame
        if len(dres.keys()) != 0:
            int_out = pd.DataFrame.from_dict(dres)
            # same columns as dres: 
            # Int range(s), target pH, t90 pH, target Sum, t90 meas, integral target Sum, integral Sum, error
            int_out.index = ['integration range (s)', 'target pH', 't90 pH', 'target Sum', 't90 meas sens', 
                             'integral target Sum', 'integral observed Sum', 'error']
        else:
            int_out = None

        # save output
        output = bs.save_report(para_meas=self.para_meas, sensor_ph=self.sensor_ph, df_res=self.df_res, dres=dres,
                                sensor_para2=self.sensor_para2, total_lbl=self.total_para_grp)

        with open(fname_save, 'w') as f:
            output[0].to_csv(f, sep='\t', decimal='.')
        with open(fname_save, 'a') as f:
            output[1].to_csv(f, sep='\t', decimal='.', header=False)
        if isinstance(int_out, type(None)):
            pass
        else:
            with open(fname_save.split('.')[0] + '_integral.txt', 'w') as fint:
                int_out.to_csv(fint, sep='\t', decimal='.', header=False)

    # ---------------------------------------------------
    def _linEdit2list(self, line):
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
                    raise ValueError(
                        'Typo error. Please check your pH entry: ' + line)
        elif ',' in line:
            ls = [float(i) for i in line.split(',') if len(i.strip()) != 0]
        else:
            ls.append(float(line))
        return ls

    def check_parameter(self):
        # get sum concentration(s) and pH value(s) and align list of fluctuation points
        ls_cSum_gmL = self._linEdit2list(line=self.paraSum_conc_edit.text())
        ls_ph = self._linEdit2list(line=self.ph_edit.text())
        ls_ph, ls_cSum_gmL = self.align_concentrations(ls_ph=ls_ph, ls_cSum_gmL=ls_cSum_gmL)

        # determine all required target concentrations
        df_target = self.parameter_prep(ls_ph=ls_ph, ls_total=ls_cSum_gmL, pKa=float(self.para2_pka_edit.text()),
                                        t_plateau=float(self.tsteady_edit.text()))

        # calculate target concentrations of sensor 2
        if self.analyte in ['base', 'Base']:
            df_target_sens2 = pd.DataFrame(df_target['signal Base'])
        elif self.analyte in ['acid', 'Acid']:
            df_target_sens2 = pd.DataFrame(df_target['signal Acid'])
        else:
            df_target_sens2 = pd.DataFrame(df_target['signal Base'])
        ls_target_sens2 = list(df_target_sens2[df_target_sens2.columns[0]].drop_duplicates().values)

        # collect all relevant parameter
        self.sensor_ph = dict({'E0': ph_E0, 't90': float(self.ph_t90_edit.text()), 'resolution': ph_res,
                               'start pH': ls_ph[0], 'time steps': tsteps, 'pH target': df_target['signal pH'], 'set values': ls_ph})

        self.sensor_para2 = dict({'analyte': self.analyte, 'time steps': tsteps,
                                  'pKa': float(self.para2_pka_edit.text()), 'E': sens_E, 'resolution': sens_res,
                                  't90': float(self.para2_t90_edit.text()), 'Base target': df_target['signal Base'],
                                  'Acid target': df_target['signal Acid'], 'Sum target': df_target['signal Sum'],
                                  'conc_calib': conc_calib, 'set values': ls_cSum_gmL, 'set sensor2': ls_target_sens2,
                                  'system select': self.system_select})
        self.para_meas = dict({'plateau time': float(self.tsteady_edit.text())})

        # -----------------------------------------------------------------------------------------------------------
        # check pH sensor
        try:
            if float(self.sensor_ph['t90']) == 0:
                msgBox = error_message(text="Please enter a valid sensor response time for the pH sensor, in "
                                            "particular a number > 0", type='Error')
                msgBox.exec_()
            else:
                try:
                    # check Base/Acid sensor
                    if float(self.sensor_para2['t90']) == 0:
                        msgBox = error_message(text="Please enter a valid sensor response time for the pH sensor, in "
                                                    "particular a number > 0", type='Error')
                        msgBox.exec_()
                    else:
                        # in case parameter check was successful, start run_simulation
                        self.run_simulation()
                except:
                    pass
        except:
            pass

    def run_simulation(self):
        # clear potential left-over plots
        self.clear_phsim()
        self.clear_para2timedrive()
        self.clear_sumtimdrive()
        global ls_lines, _integral_counter, unit
        ls_lines, self.ls_xcoords, self.vlines_list = [], [], []
        
        # target concentrations
        df_res = pd.concat([self.sensor_ph['pH target'], self.sensor_para2['Base target'], 
                            self.sensor_para2['Acid target'], self.sensor_para2['Sum target']], axis=1)
        df_res.columns = ['target pH', 'target_{} Base'.format(unit),
                          'target_{} Acid'.format(unit), 'target_{} Sum'.format(unit)]

        # individual sensor - reduce by individual plotting (different function)
        [df_res, cplateaupH, cplateauTotal,
         para_para] = bs._alignSensorSettings(sensor_ph=self.sensor_ph, sensor2=self.sensor_para2, df_target=df_res)
        df_res.index = [round(u, 2) for u in df_res.index]
        df_res = df_res.dropna()

        # pH sensor response
        [df_pHrec, df_pHcalc] = bs.pH_sensor(cplateau=cplateaupH, sensor_ph=self.sensor_ph, para_meas=self.para_meas)
        df_pHrec.index = [round(i, 2) for i in df_pHrec.index]

        # include sensor response and drift to result dataframe
        df_res = pd.concat([df_res, df_pHrec.loc[df_res.index], df_pHcalc.loc[df_res.index]], axis=1).sort_index(axis=1)

        # para2 sensor - calculate individual sensor response as well as total parameter
        df_res = bs.para2Sensorcalc(df_res=df_res, analyte=self.sensor_para2['analyte'], cplateauSum=cplateauTotal,
                                    df_pHcalc=df_pHcalc, sensor2=self.sensor_para2, para_meas=self.para_meas)
        self.df_res = df_res.sort_index(axis=1)

        # ---------------------------------------------------------------------------------------------------------
        # plotting part | individual sensor - target vs record
        _ = plot_phsensor(df_res=self.df_res, color=dcolor['font'], fig=self.fig_phsim)
        _ = plot_parasensor(df_res=self.df_res,  analyte=self.sensor_para2['analyte'], para2=self.para2,
                            para3=self.para3, color=dcolor['font'], para_sum=self.total_para, fig=self.fig_para2sim,
                            ax=self.ax_para2sim, ax1=self.ax1_para2sim, labels=self.system_select)

        # final total parameter model
        _ = plot_totalModel(df_res=self.df_res, color=dcolor['font'], para_sum=self.total_para, fig1=self.ax_totalsim,
                            labels=self.system_select)
        self.fig_totalsim.scene().sigMouseClicked.connect(self.mouseMoved)

        # --------------------------------------------------------------
        # re-draw figures for saving with transparent/white background
        fig_pH = plot_phsensor4save(df_res=self.df_res)
        fig_base = plot_paras4save(df_res=self.df_res, analyte=self.sensor_para2['analyte'], para2=self.para2,
                                   para3=self.para3, para_sum=self.total_para, labels=self.system_select)
        fig_total = plot_total4save(df_res=self.df_res, para_sum=self.total_para, labels=self.system_select)

        # --------------------------------------------------------------------------------------------------------------
        # collect for result output (save data)
        self.dic_figures = dict({'pH': fig_pH, 'Total': fig_total, 'Base': fig_base})

    def mouseMoved(self, event):
        global ls_lines
        # allow click only in combination with additionally pressed key – has to be on purpose
        modifiers = QApplication.keyboardModifiers()
        if modifiers != Qt.ControlModifier:  # change selected range
            return
        if len(self.ls_xcoords) > 2:
            self.ls_xcoords = list()
            self.int_button.setEnabled(False)
        else:
            self.int_button.setEnabled(True)

        # collect xdata in list self.ls_xcoords
        vb = self.ax_totalsim.vb
        mousePoint = vb.mapSceneToView(event._scenePos)
        self.ls_xcoords.append(mousePoint.x())

        # add vertical lines in plot
        self.v_line = pg.InfiniteLine(angle=90, movable=True)
        self.v_line.setPos(mousePoint.x())
        self.ax_totalsim.addItem(self.v_line)

        # make vertical line drag-able and update x-value in list
        self.v_line.sigPositionChanged.connect(self.movedBoundary)
        ls_lines.append(self.v_line)
        
        # remove duplicates from list
        ls_lines = list(dict.fromkeys(ls_lines))
        self.ls_xcoords = list(dict.fromkeys(self.ls_xcoords))

        # allow integral calculation as soon as xcoords have been collected
        if len(self.ls_xcoords) == 2:
            self.int_button.setEnabled(True)

    def movedBoundary(self, line):
        pos = line.pos()[0]
        line_pos = min(range(len(self.ls_xcoords)),
                       key=lambda i: abs(self.ls_xcoords[i]-pos))
        self.ls_xcoords[line_pos] = pos

        # allow integral calculation as soon as xcoords have been collected
        if len(self.ls_xcoords) == 2:
            self.int_button.setEnabled(True)

    def calc_integral(self):
        global dres, _integral_counter, ls_lines, unit
        _integral_counter += 1
        try:
            self.ls_xcoords
        except:
            self.ls_xcoords = list()

        # calculate integral for calculated/target total parameter using simpson's rule for discrete data
        # subtract calculated Sum by target Sum
        if len(self.ls_xcoords) == 2:
            self.df4int = self.df_res[['Sum calc', 'target_{} Sum'.format(unit)]].loc[min(self.ls_xcoords):max(self.ls_xcoords)]
            dfint_calc = integrate.simps(self.df4int['Sum calc'].to_numpy(), x=self.df4int.index)
            # target function assumed to have an infinite fast response (~step function)
            grp = self.df4int.groupby('target_{} Sum'.format(unit))
            dfint_target = np.sum([(k * (grp.groups[k][-1] - grp.groups[k][0])) for k in grp.groups.keys()])

            # Integration range (s), target pH, t90 pH, target Sum, t90 meas, integral target Sum, integral Sum, error
            result = list([(round(min(self.ls_xcoords), 2), round(max(self.ls_xcoords)), 2),
                           self.df_res['target pH'].loc[min(self.ls_xcoords):max(self.ls_xcoords)].mean(),
                           self.sensor_ph['t90'], self.df_res['target_{} Sum'.format(unit)].mean(), self.sensor_para2['t90'], 
                           dfint_calc, dfint_target, dfint_calc - dfint_target])

            # add current results to dictionary (where all selected peak information are stored)
            dres[_integral_counter] = result

            # open a pop up window with options to select what shall be saved
            global wInt
            wInt = IntegralWindow(self.para_meas, self.sensor_ph, self.sensor_para2, self.total_para_grp)
            if wInt.isVisible() is False:
                wInt.show()

            # clear integral area in total parameter plot
            for line in ls_lines:
                self.ax_totalsim.removeItem(line)
            [vlines.remove() for vlines in self.vlines_list]

            # update figure plot
            self.vlines_list, ls_lines = list(), list()
            self.ls_xcoords = list()

            # disable Integral button til 2 xcoords are selected
            self.int_button.setEnabled(False)


class IntegralWindow(QDialog):
    def __init__(self, para_meas, sensor_ph, sensor_para2, sum_lbl):
        super().__init__()
        self.sum_lbl = sum_lbl
        self.para_meas, self.sensor_ph, self.sensor_para2 = para_meas, sensor_ph, sensor_para2
        self.initUI()

        # when checkbox selected, save information in registered field
        self.reset_button.clicked.connect(self.reset)
        self.save_button.clicked.connect(self.save_integral)
        self.close_button.clicked.connect(self.close_window)

        # execute function as soon as the window pups up or as soon as integral_button is clicked
        self.res2table()

    def initUI(self):
        self.setWindowTitle("error calculation via Simpson's integration")
        # define geometry of integration window x, y, width, height
        self.setGeometry(50, 100, 700, 300)

        # close window button
        self.close_button = QPushButton('OK', self)
        self.close_button.setFixedWidth(100),  self.close_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.reset_button = QPushButton('Reset', self)
        self.reset_button.setFixedWidth(100), self.reset_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))
        self.save_button = QPushButton('Save', self)
        self.save_button.setFixedWidth(100), self.save_button.setFont(QFont('Helvetica Neue', int(fs_font*0.7)))

        # create table to store data
        self.tab_report = QTableWidget(self)
        self.tab_report.setColumnCount(8), self.tab_report.setRowCount(1)
        self.tab_report.setHorizontalHeaderLabels(['Integration range (s)', 'Target pH', 't90 pH', 'Target ' + self.sum_lbl,
                                                   't90 meas sens', 'Integral target ' + self.sum_lbl, 'Integral observed ' +\
                                                    self.sum_lbl, 'Error'])
        self.tab_report.resizeColumnsToContents()
        self.tab_report.resizeRowsToContents()

        # creating window layout
        mlayout2 = QVBoxLayout()
        vbox2_top, vbox2_middle, vbox2_bottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        mlayout2.addLayout(vbox2_top), mlayout2.addLayout(vbox2_middle), mlayout2.addLayout(vbox2_bottom)

        # panel for integration results
        table_grp = QGroupBox("Results")
        grid_res = QGridLayout()
        table_grp.setFont(QFont('Helvetica Neue', 10))
        vbox2_top.addWidget(table_grp)
        table_grp.setLayout(grid_res)

        grid_res.addWidget(self.tab_report)
        self.tab_report.adjustSize()

        # vertical line
        hline = QFrame()
        hline.setFrameShape(QFrame.HLine | QFrame.Raised)
        hline.setLineWidth(2)
        vbox2_middle.addWidget(hline)

        # navigation panel
        navi_grp = QGroupBox("Navigation")
        grid_navi = QGridLayout()
        navi_grp.setFont(QFont('Helvetica Neue', 10))
        vbox2_bottom.addWidget(navi_grp)
        navi_grp.setLayout(grid_navi)

        # include widgets in the layout
        grid_navi.addWidget(self.close_button, 0, 0)
        grid_navi.addWidget(self.reset_button, 0, 1)
        grid_navi.addWidget(self.save_button, 0, 2)

        # add everything to the window layout
        self.setLayout(mlayout2)

    def res2table(self):
        global dres, _integral_counter
        # get the current row by counting the clicks on the integral button
        x0 = _integral_counter-1

        # check whether row is empty (first column)
        item = self.tab_report.item(x0, 0)
        x = x0 if not item or not item.text() else x0 + 1
        item = self.tab_report.item(x, 0)

        # add the number of rows according to keys in dictionary
        self.tab_report.setRowCount(len(dres.keys()))

        # go through the dictionary and fill up the table (again)
        for c in dres.keys():
            # columns: Integration range (s), target pH, target Sum, integral target Sum, integral Sum, error
            for en in enumerate(dres[c]):
                if en[0] == 0:
                    l = str(round(en[1][0], 2)) + ' - ' + str(round(en[1][1], 2))
                else:
                    l = str(round(en[1], 2))
                item = QTableWidgetItem(l)
                item.setTextAlignment(Qt.AlignRight)

                # item structure: row, table, content
                self.tab_report.setItem(c-1, en[0], item)
        self.tab_report.resizeColumnsToContents(), self.tab_report.resizeRowsToContents()

    def reset(self):
        global dres, _integral_counter

        for c in np.linspace(1, _integral_counter, num=int((_integral_counter-1)/1+1)):
            dres.pop(c)
        _integral_counter = 0

        # when reset - reset counter and drop previously selected area from dres
        self.tab_report.clear()
        self.tab_report.setHorizontalHeaderLabels(['integration range (s)', 'target pH', 'target Sum',
                                                   'integral target Sum', 'integral observed Sum', 'error'])
        self.tab_report.setRowCount(1)

    def save_integral(self):
        # opens window in current path. User input to define file name
        fname_save = QFileDialog.getSaveFileName(self, 'Save File')[0]

        if fname_save:
            if '.txt' in fname_save or '.csv' in fname_save:
                fname_save = fname_save.split('.')[0] + '_integral.txt'
            else:
                fname_save = fname_save + '_integral.txt'

            # check whether we have data to save
            try:
                if self.result is None:
                    msgBox = error_message(
                        text="Simulate before saving. Simulation might have been unsuccessful.")
                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # save output now
                    output = bs.save_integral(para_meas=self.para_meas, sensor_ph=self.sensor_ph, dres=self.result,
                                              sensor_para2=self.sensor_para2, total_lbl=self.sum_lbl)

                    with open(fname_save, 'w') as f:
                        output.to_csv(f, sep='\t', decimal='.', header=False)

            except NameError:
                msgBox = error_message(text="Simulate before saving.")
                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    def close_window(self):
        self.hide()


# .....................................................................................................................
def plot_phsensor(df_res, color, fig=None):
    # plotting with pyqtgraph
    penT = pg.mkPen(color=color, width=1, style=Qt.DashLine)
    penS = pg.mkPen(color=dcolor['pH'], width=2)
    fig.plot(df_res['target pH'].dropna().index, df_res['target pH'].dropna().to_numpy(), pen=penT)
    fig.plot(df_res['pH calc'].dropna().index, df_res['pH calc'].dropna().to_numpy(), pen=penS)

    fig.setXRange(int(-0.5), int(df_res['pH calc'].dropna().index[-1]*1.05))
    fig.setYRange(int(df_res['pH calc'].min()), int(df_res['pH calc'].max()*1.05))
    return fig


def plot_phsensor4save(df_res):
    fig = pg.PlotWidget()
    fig.setBackground('w')

    styles = {'color': 'k', 'font-size': '15px'}
    fig.setLabel('left', '<font>pH value</font>', **styles)
    fig.setLabel('bottom', '<font>Time</font>',
                 units='<font>s</font>', **styles)

    # plotting with pyqtgraph
    penT = pg.mkPen(color=dcolor['background'], width=1, style=Qt.DashLine)
    penS = pg.mkPen(color=dcolor['pH'], width=2)
    fig.plot(df_res['target pH'].dropna().index,
             df_res['target pH'].dropna().to_numpy(), pen=penT)
    fig.plot(df_res['pH calc'].dropna().index,
             df_res['pH calc'].dropna().to_numpy(), pen=penS)

    fig.setXRange(int(-0.5), int(df_res['pH calc'].dropna().index[-1]*1.05))
    fig.setYRange(int(df_res['pH calc'].min()),
                  int(df_res['pH calc'].max()*1.05))

    # adjust color
    fig.getAxis('bottom').setPen(
        pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('bottom').setTextPen(
        pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('left').setPen(
        pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('left').setTextPen(
        pg.mkPen(color=dcolor['background'], width=1.75))

    return fig


def plot_parasensor(df_res, para2, para3, analyte, color, para_sum, labels, fig=None, ax=None, ax1=None):
    global unit
    if analyte == 'base' or analyte == 'Base':
        df_target2 = df_res['target_{} Base'.format(unit)].dropna()
        df_target3 = df_res['target_{} Acid'.format(unit)].dropna()
        df_para2 = df_res['Base calc'].dropna()
        df_para3 = df_res['Acid calc'].dropna()
        print(1544, labels)

        # specify selected system (para2, para3, para_total)
        if 'TAN' in para_sum:
            fig.setLabel('left', '<font>NH<sub>3</sub></font>', units=unit, color=dcolor['para2'])
            fig.setLabel('right', '<font>NH<sub>4</sub><sup>+</sup></font>', units=unit, color=dcolor['para3'])
        elif 'TDS' in para_sum:
            fig.setLabel('left', '<font>HS<sup>-</sup></font>', units=unit, color=dcolor['para2'])
            fig.setLabel('right', '<font>H<sub>2</sub>S</font>', units=unit, color=dcolor['para3'])
        else:
            fig.setLabel('left', '<font>{}</font>'.format(labels[0]), units=unit, color=dcolor['para2'])
            fig.setLabel('right', '<font>{}</font>'.format(labels[1]), units=unit, color=dcolor['para3']) 
    else:
        df_target2 = df_res['target_{} Acid'.format(unit)].dropna() 
        df_target3 = df_res['target_{} Base'.format(unit)].dropna()
        df_para2 = df_res['Acid calc'].dropna()
        df_para3 = df_res['Base calc'].dropna()
        if 'TAN' in para_sum:
            fig.setLabel('right', '<font>NH<sub>3</sub></font>', units=unit, color=dcolor['para2'])
            fig.setLabel('left', '<font>NH<sub>4</sub><sup>+</sup></font>', units=unit, color=dcolor['para3'])
        elif 'TDS' in para_sum:
            fig.setLabel('left', '<font>HS<sup>-</sup></font>', units=unit, color=dcolor['para2'])
            fig.setLabel('right', '<font>H<sub>2</sub>S</font>', units=unit, color=dcolor['para3'])
        else:

            fig.setLabel('left', '<font>{}</font>'.format(labels[0]), units=unit, color=dcolor['para2'])
            fig.setLabel('right', '<font>{}</font>'.format(labels[1]), units=unit, color=dcolor['para3']) 
    # plotting with pyqtgraph
    penT = pg.mkPen(color=color, width=1, style=Qt.DashLine)
    penS = pg.mkPen(color=dcolor['para2'], width=2)
    penT1 = pg.mkPen(color=color, width=1, style=Qt.DotLine)
    penS1 = pg.mkPen(color=dcolor['para3'], width=2)

    # measured analyte
    if isinstance(para2, list):
        para2 = para2[0]
    ax.plot(df_target2.index.to_numpy(), df_target2.to_numpy(), pen=penT, label=para2 + ' target')
    ax.plot(df_para2.index.to_numpy(), df_para2.to_numpy(), pen=penS, label=para2)
    
    # other analyte
    if isinstance(para3, list):
        para3 = para3[0]
    ax1.addItem(pg.PlotCurveItem(x=df_target3.index.to_numpy(), y=df_target3.to_numpy(), pen=penT1,
                                 label=para3 + ' target'))
    ax1.addItem(pg.PlotCurveItem(x=df_para3.index.to_numpy(), y=df_para3.to_numpy(), pen=penS1, label=para3))
    return fig


def plot_paras4save(df_res, para2, para3, analyte, para_sum, labels):
    global unit
    fig = pg.PlotWidget()
    fig.setBackground('w')

    ax = fig.plotItem
    ax1 = pg.ViewBox()
    fig.scene().addItem(ax1)
    ax.getAxis('right').linkToView(ax1)
    ax1.setXLink(ax)

    styles = {'color': 'k', 'font-size': '15px'}
    fig.setLabel('bottom', '<font>Time</font>', units='<font>s</font>', **styles)
    if analyte == 'base' or analyte == 'Base':
        df_target2 = df_res['target_{} Base'.format(unit)].dropna()
        df_para2 = df_res['Base calc'].dropna()
        df_target3 = df_res['target_{} Acid'.format(unit)].dropna()
        df_para3 = df_res['Acid calc'].dropna()

        if 'TAN' in para_sum:
            fig.setLabel('left', '<font>NH<sub>3</sub></font>', units=unit, **styles)
            fig.setLabel('right', '<font>NH<sub>4</sub><sup>+</sup></font>', units=unit, **styles)
        elif 'TDS' in para_sum:
            fig.setLabel('left', '<font>HS<sup>-</sup></font>', units=unit, **styles)
            fig.setLabel('right', '<font>H<sub>2</sub>S</font>', units=unit, **styles)
        else:
            fig.setLabel('left', '<font>{}</font>'.format(labels[0]), units=unit, **styles)
            fig.setLabel('right', '<font>{}</font>'.format(labels[1]), units=unit, **styles) 
    else:
        df_target3 = df_res['target_{} Base'.format(unit)].dropna()
        df_para3 = df_res['Base calc'].dropna()
        df_target2 = df_res['target_{} Acid'.format(unit)].dropna()
        df_para2 = df_res['Acid calc'].dropna()
        if 'TAN' in para_sum:
            fig.setLabel('right', '<font>NH<sub>3</sub></font>', units=unit, **styles)
            fig.setLabel('left', '<font>NH<sub>4</sub><sup>+</sup></font>', units=unit, **styles)
        elif 'TDS' in para_sum:
            fig.setLabel('left', '<font>HS<sup>-</sup></font>', units=unit, **styles)
            fig.setLabel('right', '<font>H<sub>2</sub>S</font>', units=unit, **styles)
        else:
            fig.setLabel('left', '<font>{}</font>'.format(labels[0]), units=unit, **styles)
            fig.setLabel('right', '<font>{}</font>'.format(labels[1]), units=unit, **styles) 

    # plotting with pyqtgraph
    penT = pg.mkPen(color=dcolor['background'], width=1, style=Qt.DashLine)
    penS = pg.mkPen(color=dcolor['para2'], width=2)
    penT1 = pg.mkPen(color=dcolor['background'], width=1, style=Qt.DotLine)
    penS1 = pg.mkPen(color=dcolor['para3'], width=2)

    # measured analyte
    ax.plot(df_target2.index.to_numpy(), df_target2.to_numpy(), pen=penT, label=para2 + ' target')
    ax.plot(df_para2.index.to_numpy(), df_para2.to_numpy(), pen=penS, label=para2)
    # other analyte
    ax1.addItem(pg.PlotCurveItem(x=df_target3.index.to_numpy(), y=df_target3.to_numpy(), pen=penT1, 
                                 label=para3 + ' target'))
    ax1.addItem(pg.PlotCurveItem(x=df_para3.index.to_numpy(), y=df_para3.to_numpy(), pen=penS1, label=para3))

    def updateViews():
        # view has resized; update auxiliary views to match
        ax1.setGeometry(ax.vb.sceneBoundingRect())
    updateViews()
    ax.vb.sigResized.connect(updateViews)

    # adjust color
    fig.getAxis('top').setPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('bottom').setPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('bottom').setTextPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('left').setPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('left').setTextPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('right').setPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('right').setTextPen(pg.mkPen(color=dcolor['background'], width=1.75))
    return fig


def plot_totalModel(df_res, para_sum, color, labels, fig1=None):
    global unit
    # plotting with pyqtgraph
    df_target = df_res['target_{} Sum'.format(unit)].dropna()
    df_sim = df_res['Sum calc'].dropna()

    if 'TAN' in para_sum:
        fig1.setLabel('left', '<math>TAN</math>', units='<font>{}</font>'.format(unit), color=dcolor['font'])
    elif 'TDS' in para_sum:
        fig1.setLabel('left', '<font>TDS</font>', units='<font>{}</font>'.format(unit), color=dcolor['font'])
    else:
        fig1.setLabel('left', '<font>{}</font>'.format(labels[-1]), units='<font>{}</font>'.format(unit), color=dcolor['font'])

    penT = pg.mkPen(color=color, width=1, style=Qt.DashLine)
    penS = pg.mkPen(color=dcolor['total para'], width=2)
    fig1.plot(list(df_target.index), list(df_target.to_numpy()), pen=penT, label=para_sum + ' target')
    fig1.plot(df_sim.index, df_sim.to_numpy(), pen=penS, label=para_sum)

    # set initial view
    fig1.setXRange(int(-0.5), int(df_sim.index[-1]*1.05))
    if df_sim.max()*1.05 < int(df_target.max()*1.85):
        if df_sim.min() > int(df_target.min()*0.15):
            fig1.setYRange(int(df_sim.min()), df_sim.max()*1.05)
        else:
            fig1.setYRange(int(df_target.min()*0.15), df_sim.max()*1.05)
    else: 
        if df_sim.min() > int(df_target.min()*0.15):
            fig1.setYRange(int(df_sim.min()), int(df_target.max()*1.85))
        else:
            fig1.setYRange(int(df_target.min()*0.15), int(df_target.max()*1.85))

    return fig1


def plot_total4save(df_res, para_sum, labels):
    global unit
    fig = pg.PlotWidget()
    fig.setBackground('w')

    styles = {'color': 'k', 'font-size': '15px'}
    fig.setLabel('bottom', '<font>Time</font>', units='<font>s</font>', **styles)
    if 'TAN' in para_sum:
        fig.setLabel('left', '<math>TAN</math>', units='<font>{}</font>'.format(unit), **styles)
    elif 'TDS' in para_sum:
        fig.setLabel('left', '<font>TDS</font>', units='<font>{}</font>'.format(unit), **styles)
    else:
        fig.setLabel('left', '<font>{}</font>'.format(labels[-1]), units='<font>{}</font>'.format(unit), **styles)

    # plotting with pyqtgraph
    df_target = df_res['target_{} Sum'.format(unit)].dropna()
    df_sim = df_res['Sum calc'].dropna()

    penT = pg.mkPen(color=dcolor['background'], width=1, style=Qt.DashLine)
    penS = pg.mkPen(color=dcolor['total para'], width=2)
    fig.plot(list(df_target.index), list(df_target.to_numpy()), pen=penT, label=para_sum + ' target')
    fig.plot(df_sim.index, df_sim.to_numpy(), pen=penS, label=para_sum)

    # set initial view
    fig.setXRange(int(-0.5), int(df_sim.index[-1]*1.05))

    # adjust color
    fig.getAxis('bottom').setPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('bottom').setTextPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('left').setPen(pg.mkPen(color=dcolor['background'], width=1.75))
    fig.getAxis('left').setTextPen(pg.mkPen(color=dcolor['background'], width=1.75))
    return fig


# .....................................................................................................................
def error_message(text, type='Warning'):
    msgBox = QMessageBox()
    msgBox.setIcon(QMessageBox.Information)
    msgBox.setText(text)
    msgBox.setFont(QFont('Helvetica Neue', 11))
    msgBox.setWindowTitle(type)
    msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

    return msgBox


# .....................................................................................................................
if __name__ == '__main__':
    import sys

    app = QtWidgets.QApplication(sys.argv)

    path = os.path.join(loc_path + r'/picture/icon.png')
    app.setWindowIcon(QIcon(path))

    # options available: 'Breeze', 'Oxygen', 'QtCurve', 'Windows', 'Fusion'
    app.setStyle('Fusion')

    # screen Size adjustment
    Wizard = MagicWizard()
    screen = app.primaryScreen()
    rect = screen.availableGeometry()
    Wizard.setMaximumHeight(int(rect.height() * 0.9))
    # position (screen width, screen height), width, height
    Wizard.setGeometry(int(rect.width() * 0.075), int(rect.height() * 0.1), 710, 550)

    # show wizard
    Wizard.show()
    sys.exit(app.exec_())
