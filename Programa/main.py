# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 11:42:55 2022

@author: Anderson Sampaio
"""

#%% PROGRAMA PRINCIPAL



#%% Leitura dos Dados
from Read_Data import Read_Data

# IMPORTANTE: Para que o programa funcione é preciso que altere o endereço abaixo de acordo com endereço onde
# salvou o arquivo em seu PC
path_pwf = r"G:\Meu Drive\1- DOUTORADO - UFJF\2021.3\Estudo Dirigido\Flow - Newton Raphson - Leitura Excel\CASES.xlsx"

caso = 1  # 1 (24Barras)    2 (4barras)
Dados_Sist = Read_Data(path_pwf, caso)


#%% Fluxo de Potência Newton Raphson

from Flow_Newton_Raphson import Flow
Flupot = Flow(Dados_Sist.Sbase, Dados_Sist.df_DBAR, Dados_Sist.df_DLIN)