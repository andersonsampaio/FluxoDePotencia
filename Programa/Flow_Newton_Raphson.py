# -*- coding: utf-8 -*-
"""
Created on Wed Feb 05 11:42:55 2022

@author: Anderson Sampaio
"""

import numpy as np
import math
from termcolor import colored


class Flow(object):
    def __init__(self, Sbase, df_DBAR, df_DLIN):
        
        #%% Calculo de Ybarra
        
        # Linhas e transformadores 
        n_lin = len(df_DLIN["DE"]) #tamanho da primeira dimensao
        De = df_DLIN["DE"]
        Para = df_DLIN["PARA"]
               
        n_elem1 = int(max(max(De),max(Para)))   
        ybus = np.zeros((n_elem1 ,n_elem1),dtype = "complex_") 
        
        vetG = []
        vetB = []
        vetBsh = []
        vetA = np.ones(n_lin)
        
        default = np.zeros(n_lin) 
        
        for l in range (0,n_lin):
            
            RelTrans = df_DLIN["TAP"].values[l]
            default[l]=RelTrans
            
            K = De[l]
            k = int(K)
            M = Para[l]
            m = int(M)
            
            yserie = (1/(df_DLIN["R"].values[l]/Sbase + 1j*df_DLIN["X"].values[l]/Sbase))
            vetG.append(np.real(yserie))
            vetB.append(np.imag(yserie))
            
            yshunt = (1j*df_DLIN["Bc"].values[l]/(2*Sbase))
            vetBsh.append(df_DLIN["Bc"].values[l]/(2*Sbase))
            
        
        
            if RelTrans == 0.0:
                ybus[k-1][k-1] =  ybus[k-1][k-1] + yserie + yshunt
                ybus[m-1][m-1] =  ybus[m-1][m-1] + yserie + yshunt
                ybus[k-1][m-1] =  ybus[k-1][m-1] - yserie 
                ybus[m-1][k-1] =  ybus[m-1][k-1] - yserie
        
            else:
                a = 1/RelTrans
                A = a * yserie
                B = a * yserie * (a - 1)
                C = yserie * (1 - a)
                
                ybus[k-1][k-1] =  ybus[k-1][k-1] + A + B
                ybus[m-1][m-1] =  ybus[m-1][m-1] + A + C
                ybus[k-1][m-1] =  ybus[k-1][m-1] - A 
                ybus[m-1][k-1] =  ybus[m-1][k-1] - A
            
                vetA[l]=a    
               
        # Bancos de capacitor e reator 
                
        n_barra = len(df_DBAR["Numero"])
        for c1 in range(0, n_barra):
            if df_DBAR["Sh"].values[c1] != 0.0:
                K1 = df_DBAR["Numero"].values[c1]
                k1 = int(K1)
        
                Bbanco =  (1j*df_DBAR["Sh"].values[c1]/Sbase)
            
                ybus[k1-1][k1-1] =  ybus[k1-1][k1-1] + Bbanco
        
        # Ybus = G + j*B
        Gy = ybus.real
        By = ybus.imag
        
        #%% Métodos para calculo das equações de potência
                
        X = np.zeros(2*n_barra)                                 # Criando vetor X no qual será armazendo as incógnitas do Subsistema 1 
                                                                # disposição do vetor coluna X -> X=[TetaK1 Vk1 TetaK2 Vk2 ... TetaKn Vkn]
                                
        for c4 in range (0, n_barra):                           # Salvando incógnitas do Subsistema 1 em X
            X[2*c4] = df_DBAR["Fase"].values[c4]*math.pi/180              # Salva Teta e converte para radianos
            X[2*c4+1] = df_DBAR["Tensao"][c4]  
                
        # Método para cálculo da equação de POTÊNCIA ATIVA     # dada pela equação:  Pk = Somatório[ Vk*Vm*(G_km*cos(Tk-Tm) + B_km*sin(Tk-Tm)) ]   
        def PCalc(barra): 
        
            pcalc = 0                                           # Zerando variável ao inicializar o método
            
            Tk = X[2*barra]                                     # Colhendo o Tk da barra atual
            Vk = X[2*barra+1]                                   # Colhendo a Vk da barra atual
            
        
            for c2 in range(0, n_barra):                        # Laço para criar o Somatório
                Tm = X[2*c2]                                    # Colhendo o Tm da barra atual até n barras para criar o somatório (for)
                Vm = X[2*c2+1]                                  # Colhendo o Vm da barra atual até n barras para criar o somatório (for)
                
                G = np.real(ybus[barra,c2])                     # Colhendo a CONDUTÅNCIA G_km na YBUS
                B = np.imag(ybus[barra,c2])                     # Colhendo a SUSCEPTÂNCIA B_km na YBUS
            
                P = Vk*Vm*(G*np.cos(Tk-Tm) + B*np.sin(Tk-Tm))   # Calculando equação Potência Ativa
                pcalc += P                                      # Fazendo Somatório
        
            return pcalc
         
        # Método para cálculo da equação de POTÊNCIA REATIVA    # dada pela equação:  Qk = Somatório[ Vk*Vm*(G_km*.sin(Tk-Tm) - B_km*cos(Tk-Tm))]      
        def QCalc(barra):
        
            qcalc = 0
                 
            Tk = X[2*barra]
            Vk = X[2*barra+1]  
        
            
            for c2 in range(0, n_barra):
                Tm = X[2*c2]
                Vm = X[2*c2+1]
                
                G = np.real(ybus[barra,c2])
                B = np.imag(ybus[barra,c2])
                
                Q = Vk*Vm*(G*np.sin(Tk-Tm) - B*np.cos(Tk-Tm))   # Calculando equação Potência Reativa
                qcalc += Q                                      # Fazendo Somatório
        
            return qcalc
        
        #%% Calculo do SUBSISTEMA 1
                
        ite = 0
        iteMax = 11
        tol = 1e-5
        
        continua = True
        while continua:
        
        
        # Calculo de delta P e delta Q                              
        
            deltaY = np.zeros(2*n_barra)                            # vetor de Deltas_P e Deltas_Q
            vetPcalc = np.zeros(n_barra)                            # vetor de Ps calculados
            vetQcalc = np.zeros(n_barra)                            # vetor de Qs calculados
            Pesp_ = np.zeros(n_barra)
            Qesp_ = np.zeros(n_barra)
        
        
        # Calculo de delta P e delta Q    
            for c3 in range (0, n_barra):
                
                # Calculo de DELTA P        
                Pcalc = PCalc(c3)
                if df_DBAR["Tipo"].values[c3] == 1 or df_DBAR["Tipo"].values[c3] == 0 : # Apenas barras PV ==1 ou PQ ==0
                    Pesp = df_DBAR["PG"].values[c3] - df_DBAR["PL"].values[c3]          # P_esp = PG - PD
                    deltaP = (Pesp/Sbase) - Pcalc                                       # calculo de Delta Pk = P_esp - Pk(P_calculado)
                    Pesp_[c3] = Pesp/Sbase
                else: 
                    deltaP = 0                                                          # caso seja Barra Swing
                
                # Calculo de DELTA Q
                Qcalc = QCalc(c3)       
                if df_DBAR["Tipo"].values[c3] == 0 :                                    # Apenas barras PQ ==0
                    Qesp = df_DBAR["QG"].values[c3] - df_DBAR["QL"].values[c3]          # Q_esp = QG - QD
                    deltaQ = (Qesp/Sbase) - Qcalc                                       # calculo de Delta Qk = Q_esp - Qk(Q_calculado)
                    Qesp_[c3] = Qesp/Sbase
                else:
                    deltaQ = 0                                                          # caso seja Barra Swing ou Bara PV
                    
                # Salvando delta_P e Delta_Q
                # disposição do vetor coluna delta -> deltaY=[Delta_Pk1 Delta_Qk1 Delta_Pk2 Delta_Qk2 ... Delta_Pkn Delta_Qkn]
                deltaY[2*c3] = deltaP                               
                deltaY[2*c3+1] = deltaQ
                
                # Salvando valores dos P_calculados e Q_calculados
                vetPcalc[c3] = Pcalc
                vetQcalc[c3] = Qcalc
                
                
            if max(abs(deltaY)) > tol and ite < iteMax :            # Teste de convergencia 
                ite += 1                                            # Contador
                
                Jac = np.zeros((2*n_barra,2*n_barra))               # Incia Jacobiana considerando a estrutura BIG NUMBER
                
                for k2 in range (0,n_barra):
                    Tk2 = X[2*k2]
                    Vk2 = X[2*k2+1] 
                    
                    for m2 in range (0,n_barra):
                        Tm2 = X[2*m2]
                        Vm2 = X[2*m2+1]
                        
                        G = np.real(ybus[k2,m2])
                        B = np.imag(ybus[k2,m2])
                        
                        if k2 == m2:                                # Caso seja para a Diagonal Principal
                            BigN = 10**10                           # Valor BIG NUMBER
                            
                            # Equações para Diagonal Principal
                            H = -Vk2*Vk2*B - vetQcalc[k2]
                            N = (vetPcalc[k2] + Vk2*Vk2*G)/Vk2
                            M = vetPcalc[k2]-Vk2*Vk2*G 
                            L = (vetQcalc[k2] - Vk2*Vk2*B)/Vk2
                            
                            if df_DBAR["Tipo"].values[k2] == 2:     # Caso seja uma Barra Swing
                                H = L = BigN                        # Aplica BIG NUMBER em H e L
                            elif df_DBAR["Tipo"].values[k2] == 1:   # Caso seja uma Barra PV
                                L = BigN                            # Aplica BIG NUMBER em L
        
                        else:
                            H = Vk2*Vm2*(G*np.sin(Tk2-Tm2)-B*np.cos(Tk2-Tm2))
                            N = Vk2*(G*np.cos(Tk2-Tm2) + B*np.sin(Tk2-Tm2))
                            M = -Vk2*Vm2*(G*np.cos(Tk2-Tm2)+B*np.sin(Tk2-Tm2))
                            L = Vk2*(G*np.sin(Tk2-Tm2)-B*np.cos(Tk2-Tm2))                        
        
                        #Monta Jocobiana estrutura BIG NUMBER
                        Jac[2*k2,2*m2] = H
                        Jac[2*k2,2*m2+1] = N
                        Jac[2*k2+1,2*m2] = M
                        Jac[2*k2+1,2*m2+1] = L
                        
                     
                 
                inv_Jac = np.linalg.inv(Jac)                        # Inversa da Jacobiana
                deltaX = inv_Jac.dot(deltaY)                        # Calculo do Vetor de Correção
                
                # Correção das Variáveis de estado
                for c4 in range (0,n_barra):                        
                    X[2*c4] += deltaX[2*c4]                         # Correção de TETA -> teta^(h+1)=teta^(h) + delta_teta^(h)
                    X[2*c4+1] += deltaX[2*c4+1]                     # Correção de V -> V^(h+1)=V^(h) + delta_V^(h)
            else:
                continua = False
            
        #%% Calculo do SUBSISTEMA 2
        
        for i in range(0,n_barra):
            if df_DBAR["Tipo"].values[i] == 2 or df_DBAR["Tipo"].values[i] == 1:
                Qesp_[i] = 0
                if df_DBAR["Tipo"].values[i] == 2:
                    Pesp_[i] = 0            
        
        
        # Calcula Pk e Qk para VTETA
        # Calcula Qk para VTETA
        for k2 in range (0,n_barra):
            Tk2 = X[2*k2]
            Vk2 = X[2*k2+1] 
                    
            for m2 in range (0,n_barra):
                Tm2 = X[2*m2]
                Vm2 = X[2*m2+1]
                
                Gkm = Gy[k2,m2]
                Bkm = By[k2,m2]
                
                if df_DBAR["Tipo"].values[k2] == 2 or df_DBAR["Tipo"].values[k2] == 1: # Apenas barras PV ==1 ou Vteta ==2
                    Qesp_[k2] = Qesp_[k2] + Vk2*Vm2*(Gkm*np.sin(Tk2-Tm2)- Bkm*np.cos(Tk2-Tm2))
                    if df_DBAR["Tipo"].values[k2] == 2:
                        Pesp_[k2] = Pesp_[k2] + Vk2*Vm2*(Gkm*np.cos(Tk2-Tm2)+ Bkm*np.sin(Tk2-Tm2))
         
        PD = np.zeros(n_barra)
        QD = np.zeros(n_barra)
        Vk2 = np.zeros(n_barra)
        BSH = np.zeros(n_barra)
        for c7 in range(0,n_barra):
            PD[c7] = df_DBAR["PL"].values[c7]/Sbase
            QD[c7] = df_DBAR["QL"].values[c7]/Sbase
            Vk2[c7] = X[2*c7+1] 
            BSH[c7] = df_DBAR["Sh"].values[c7]/Sbase
        
        # Calculo das GERAÇÕES
        PG = (Pesp_ + PD)*Sbase
        QG = (Qesp_ + QD)*Sbase
        
        # SHUNT DE BARRA
        BUS_SH = (BSH*(Vk2*Vk2))*Sbase        
        
        #%%Cálculo do Fluxo Ativo e Reativo nas linhas
        
        Pkm_saida = np.zeros(n_lin)
        Pmk_saida = np.zeros(n_lin)
        Qkm_saida = np.zeros(n_lin)
        Qmk_saida = np.zeros(n_lin)
        
        
        for c5 in range (0,n_lin):
             
            k3 = int(De[c5]-1) 
            m3 = int(Para[c5]-1)
            
            Tk3 = X[2*k3]
            Vk3 = X[2*k3+1]
            Tm3 = X[2*m3]
            Vm3 = X[2*m3+1]
            
            g = vetG[c5]
            b = vetB[c5]
            bsh = vetBsh[c5]
            a = vetA[c5]
            df = 0 # trafo defasador
            
            Pkm = a**2*g*Vk3**2 - a*Vk3*Vm3*(g*np.cos(Tk3-Tm3 + df) + b*np.sin(Tk3-Tm3 +df)) 
            Qkm = -a**2*(bsh +b)*Vk3**2 + a*Vk3*Vm3*(b*np.cos(Tk3-Tm3 + df) - g*np.sin(Tk3-Tm3 +df))       
          
            Pmk = (1)**2*g*Vm3**2 - (a)*Vm3*Vk3*(g*np.cos(Tm3-Tk3 + df) + b*np.sin(Tm3-Tk3 +df)) 
            Qmk = -(1)**2*(bsh +b)*Vm3**2 + (a)*Vm3*Vk3*(b*np.cos(Tm3-Tk3 + df) - g*np.sin(Tm3-Tk3 +df)) 
            
            Pkm_saida[c5] = Pkm*Sbase
            Pmk_saida[c5] = Pmk*Sbase
            Qkm_saida[c5] = Qkm*Sbase
            Qmk_saida[c5] = Qmk*Sbase
        
        Perdas_Elementos = (Pkm_saida + Pmk_saida) * Sbase  # Perdas em cada elemento
        Perdas_Totais = sum(Perdas_Elementos)  # Perda total no sistema
        
        print(colored('-------------------------------------------------------------------------------------------------------', 'red'))
        print(colored('                                               RESULTADOS                                       ', 'red'))
        print(colored('-------------------------------------------------------------------------------------------------------', 'red'))
        if int(ite) == iteMax:
            print('  Situação:', colored("O FLUXO DIVERGIU!", 'red'))
            print("  Número de Iterações: %i" % ite)
            print(colored('________________________________________________________________________________________', 'blue'))
        else:
            print('  Situação:', colored("O FLUXO CONVERGIU!", 'green'))
            print("  Número de Iterações: %i" % ite)
            print("  Perdas no sistema: %4.4f MW" % Perdas_Totais)
            
            print(colored('_______________________________________________________________________________________________________','red'))
            print(colored('                                         RELATÓRIO DAS BARRAS                                  ','red'))
            print()
            print(colored('  Barra  Tipo  V (pu)     Ângulo (°)    PG (MW)     QG (Mvar)     PD (MW)     QD (Mvar)   SHUNT(MVAr)  ', 'red'))
            for b in range(n_barra):
                print(f"{'%5.0f' % df_DBAR['Numero'].values[b]} {'%3.0f' % df_DBAR['Tipo'].values[b]:>5} {'%3.2f' % X[2*b+1]:>8}"
                      f"{'%5.2f' % float(X[2*b]*180/np.pi):>13} {'%9.2f' % PG[b]:>13} {'%9.2f' % QG[b]:>11}"
                      f"{'%9.2f' % df_DBAR['PL'].values[b]:>13} {'%9.2f' % df_DBAR['QL'].values[b]:>11} {'%9.2f' % BUS_SH[b]:>13}")
            print(colored('=======================================================================================================','red'))
            
            print(colored('___________________________________________________________________','red'))
            print(colored('                FLUXO DE POTÊNCIA NAS LINHAS                                 ','red'))
            print(colored('-------------------------------------------------------------------', 'red'))
            print(colored('     DE   PARA    Pkm(MW)     Qkm(MVAr)     Pmk(MW)     Qmk(MVAr)  ', 'red'))
            for b in range(n_lin):
                print(f"{' %5.0f ' %De[b]} {' %3.0f ' % Para[b]} {'%3.2f' % Pkm_saida[b]:>10} {'%5.2f ' % Qkm_saida[b]:>13}"
                      f"{'%9.2f' % Pmk_saida[b]:>12} {'%9.2f' % Qmk_saida[b]:>12}")
            print(colored('___________________________________________________________________','red'))
            


            
            