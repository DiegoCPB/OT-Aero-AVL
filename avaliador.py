# -*- coding: utf-8 -*-
"""
Created on Sun Feb 08 12:32:47 2015

@author: Diego Chou
"""
print("\nCarregando modulos de 'Avaliador'...")
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import subprocess as sp
    import os.path
    print("Modulos de 'Avaliador' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Avaliador'\n")
    raise

import construtor as con
from apoio import issueCmd

"""
Esse arquivo define os avaliadores primários, ou seja,
classes que realizam uma análise para ser iterada com 
o algoritmo genético, afim de alcançar uma configuração ótima, 
ou próxima dela.
"""

class Avaliador2016(con.Construtor2016):
    """
    Essa classe avalia o indivíduo
    """
    trim = 3.0
    
    def __init__(self,name,dz_asas, vel, 
                 x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                 x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat, p):             
        con.Construtor2016.__init__(self,name,dz_asas, vel, 
                                    x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                                    x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat)  
        if p:
            self.plot()
        
        self.ps = sp.Popen(['avl.exe'],stdin=sp.PIPE,stdout=None,stderr=None)
        self.aerodinamica(self.name,self.vel,self.trim)
                                    
    def aerodinamica(self,name,vel,trim):                                                  
        """
        Analise aerodinâmica do avião gerado.
        Essa funçao executa o programa avl.exe
        """
        ps = self.ps
                
        issueCmd(ps,'load Runs/%s.avl' %(name))
        issueCmd(ps,'mass Runs/%s.mass' %(name))
        issueCmd(ps,'mset 1')
        issueCmd(ps,'oper')
        issueCmd(ps,'a')
        issueCmd(ps,'a')
        issueCmd(ps,'%f' %(trim))
        issueCmd(ps,'m')
        issueCmd(ps,'v')
        issueCmd(ps,'%f' %(vel))
        issueCmd(ps,'')
        issueCmd(ps,'x')
        issueCmd(ps,'st')
        issueCmd(ps,'Runs/%s.st' %(name))
        if os.path.isfile('Runs/%s.st' %(name)):
            issueCmd(ps,'O')
        issueCmd(ps,'')
        issueCmd(ps,'quit')
        
        
    def plot(self):
        """
        Função de debug pra verificar em um gráfico 3d se os pontos foram
        gerados nas posições corretas.
        """
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        
        asaf, asat = self.formato()[2:]        
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d([-1.25,1.25])
        ax.set_ylim3d([-1.25,1.25])
        ax.set_zlim3d([-1.25,1.25])
        
        cg = self.pos_cg
        ax.plot([cg[0]],[cg[1]],[cg[2]],'ro', label='CG')        
        
        for i in range(len(asaf)):
            M = asaf[i]
            x = np.append(M[:,0],M[0,0])
            y = np.append(M[:,1],M[0,1])
            z = np.append(M[:,2],M[0,2])
            ax.plot(x,y,z)   
                    
        for i in range(len(asat)):
            M = asat[i]
            x = np.append(M[:,0],M[0,0])
            y = np.append(M[:,1],M[0,1])
            z = np.append(M[:,2],M[0,2])
            ax.plot(x,y,z)
        
        plt.legend(loc='best')
        plt.show()
        
if __name__ == '__main__':
    name = 'A2016'
    dz_asas = 0.1
    vel = 15
    
    #Asa frontal
    x_ba_asaf = -0.5
    c_asaf = 0.35
    ang_asaf = 5.0
    epsilon_asaf = 0.0
    perfilr_asaf = 'E423'
    perfilp_asaf = 'E423'#'MIN ponta2016'
    
    # Asa traseira
    x_bf_asat = 0.8
    c_asat = 0.35
    ang_asat = 3.0
    epsilon_asat = 0.0
    perfilr_asat = 'S1223 MOD2015'
    perfilp_asat = 'S1223 MOD2015'#'MIN ponta2016'
    
    aviao = Avaliador2016(name,dz_asas, vel, 
                          x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                          x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat, True)
