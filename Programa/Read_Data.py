# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 11:46:53 2022

@author: Anderson Sampaio
"""

import pandas as pd

class Read_Data(object):
    def __init__(self, path_pwf, caso):
        self.path_pwf = path_pwf
        self.df_DBAR = None
        self.df_DLIN = None
        self.Sbase = 100
        
        # IEEE 24 BARRAS
        if caso == 1:
            self.df_DLIN = pd.read_excel(path_pwf, sheet_name="DLIN_24Barras")      
            self.df_DBAR = pd.read_excel(path_pwf, sheet_name="DBAR_24Barras")  
        # 4 BARRAS - MELHOR VISUALIZAÇÃO DAS VARIÁVEIS NO PROCESSO
        if caso == 2:
            self.df_DLIN = pd.read_excel(path_pwf, sheet_name="DLIN_4Barras")     
            self.df_DBAR = pd.read_excel(path_pwf, sheet_name="DBAR_4Barras")  

