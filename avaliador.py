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
    
    def __init__(self,name,ang_clear,dz_asas, vel, 
                 cr_asaf,ct_asaf,ang_asaf,enflex_asaf,epsilon_asaf,
                 perfilr_asaf, perfilp_asaf,
                 cr_asat,ang_asat,perfilr_asat, perfilp_asat):             
        con.Construtor2016.__init__(self,name,ang_clear,dz_asas, vel, 
                                    cr_asaf,ct_asaf,ang_asaf,enflex_asaf,epsilon_asaf,
                                    perfilr_asaf, perfilp_asaf,
                                    cr_asat,ang_asat,perfilr_asat, perfilp_asat)  
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
    ang_clear = 4.0
    dz_asas = 0.4
    vel = 15
    perfilr_asaf = 'S1223 MOD2015'
    perfilp_asaf = 'MIN ponta2016'
    perfilr_asat = 'E423'
    perfilp_asat = 'MIN ponta2016'
    cr_asaf = 0.5
    ct_asaf = 0.25
    ang_asaf = 5.0
    enflex_asaf = 37.0
    epsilon_asaf = -3.0
    cr_asat = 0.3
    ang_asat = 3.0
    
    # Avião possível
    # ['A2016',5,0.35,15,
    #  0.5,0.25,5.0,35.0,-3.0,
    #  'S1223 MOD2015','MIN ponta2016',
    #  0.3,2.0,'E423','MIN ponta2016']
    
    aviao = Avaliador2016(name,ang_clear,dz_asas, vel, 
                          cr_asaf,ct_asaf,ang_asaf,enflex_asaf,epsilon_asaf,
                          perfilr_asaf, perfilp_asaf,
                          cr_asat,ang_asat,perfilr_asat, perfilp_asat)
                           
    aviao.plot()