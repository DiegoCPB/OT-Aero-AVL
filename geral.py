# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 08:42:38 2015

@author: Diego Chou, Mateus Pereira
"""


print("\nCarregando modulos de 'Geral'...")

try:
    import numpy as np
    print("Modulos de 'Geral' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Geral'\n")
    raise

"""
 Todos os pesos estão em kg. 
 Todas as distâncias em metros.
 Todos as funções retornam angulos em graus.
"""

# Gravidade em São José dos Campos (m/s^2)
# FONTE (http://www.sismetra.cta.br/Apresentacao_VI_Semetra/Os_Efeitos_da_Gravidade_na_Calibracao_de_Balancas.pdf)
gravidade = 9.787902

# Velocidade do som (m/s)
vel_som = 335.68

class Ar(object):
    """
    Ar para a altura de voo.
    """
    def __init__(self):
        self.alt = 1400
        self.mi = 1.5e-5
        self.modelo = 'FAR'
        self.rho_mar = 1.218974884615383
        
    def rho(self):
        alt = self.alt
        rho_mar = self.rho_mar
        
        if self.modelo == 'FAR':
            # Interpolação dos pontos presentes na FAR
            rho = rho_mar-0.00011516438231497067*alt+\
                  3.9116307633987786e-9*alt**2-4.184870325042321e-14*alt**3
                  
        elif self.modelo == 'IAS':
            # Consultar International Standard Atmosphere na Wikipedia
            # Formula válida para a troposfera (até 11000m)
            rho = rho_mar*((288-6.5e-3*alt)/288)**(5.253283302063791)          
            
        return rho
  
      
class Ar_qualquer(Ar):
    """
    Ar para qualquer altura.
    """
    def __init__(self,alt):
        self.alt = alt
        self.mi = 1.5e-5
        self.modelo = 'FAR'
        self.rho_mar = 1.218974884615383


class Placa_plana(object):
    
    def __init__(self, Re):
        self.Re = Re
        
    def Cf(self):
        Re = self.Re
        
        if Re < 5e5: # Estabelecido por Raymer
            print("\nFoi selecionado Cf laminar")
            Cf = 1.328/np.sqrt(Re)
        else:
            print("\nFoi selecionado Cf turbulento")
            Cf = 0.455/(np.log10(Re)**2.58) # O nº de Mach é considerado nulo
            
        return Cf