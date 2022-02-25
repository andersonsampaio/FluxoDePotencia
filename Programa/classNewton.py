import cmath as cmt #Biblioteca números complexos
import math as mt
import numpy as np #Comando para instalar "pip install numpy"
import matplotlib.pyplot as plt

class Newton:
    def __init__(self): #Método construtor que será executado automaticamente
        self.__Sbase = 100e6
        self.__dados = dict()
        self.__Ligacoes = dict()

        self.__Sesp = dict()

        self.__tensaoPlot = dict()
        self.__angPlot = dict()
        pass

    def setBarras(self, barra, code, tensao, ang, carga, geracao):  #Método Instanciar Barras 
        """  
        Barra -> nº identificação da barra
        codes: 1 -> Tensão e Ângulo; 2 -> PQ; 3 -> PV;   
        ang -> graus
        tensao -> pu
        carga e geracao -> VA
        """
        self.__dados[barra] = {
            'code': code,
            'tensao': tensao,
            'ang': mt.radians(ang),
            'carga': (carga / self.__Sbase),
            'geracao': (geracao / self.__Sbase)
        }

        self.__tensaoPlot[barra] = [tensao] #Plotar tensão
        self.__angPlot[barra] = [ang] #Plotar Ângulo

    def printBarras(self):
        """
        Método apenas exibe informações na tela
        """
        print('\n\n================================= DADOS: =================================')
        print('Sbase = ', self.__Sbase, ' VA')
        for i in self.__dados:
            print(self.__dados[i])
        print('==========================================================================')

    def setSesp(self):
        """
        Método utilizado para calcular a potência especificada em cada barra. Automaticamente
        """

        for i in self.__dados:
            if self.__dados[i]['code'] == 2: #Barra PQ
                self.__Sesp[i] = {
                    'Pesp': np.real(self.__dados.get(i)['geracao'] - self.__dados.get(i)['carga']),
                    'Qesp': float(np.imag(self.__dados.get(i)['geracao'] - self.__dados.get(i)['carga']))
                }
            elif self.__dados[i]['code'] == 3: #Barra PV
                self.__Sesp[i] = {
                    'Pesp': np.real(self.__dados.get(i)['geracao'] - self.__dados.get(i)['carga']),
                    'Qesp': float(np.imag(self.__dados.get(i)['geracao'] - self.__dados.get(i)['carga']))
                }
        print('\n\n================================= Sesp: =================================')
        print(self.__Sesp, ' pu')
        print('==========================================================================')

    def ligacoes(self, barra1, barra2, impedancia=None, admitancia=None):
        """
        As informações devem estar em pu
        """
        if impedancia is None:
            impedancia = 1 / admitancia
        elif admitancia is None:
            admitancia = 1 / impedancia
        else:
            return 'ERRO! A admitância ou a impedância devem ser informadas.'

        self.__Ligacoes[(barra1, barra2)] = {
            'Impedância': impedancia,
            'Admitância': admitancia
        }

    def printLigacoes(self):
        print('\n\n================================= LIGAÇÕES: =================================')
        for i in self.__Ligacoes:
             print('Ligação ', i, '\t',self.__Ligacoes[i])
        print('==========================================================================')
    