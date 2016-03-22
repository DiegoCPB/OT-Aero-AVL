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
from apoio import issueCmd, executarNaPasta

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
    case_alphas = [0.0,5.0,10.0] #Angulos de ataque de análise
    number_trim_cases = 2 # Números de caso em trimagem
    
    def __init__(self,name,dz_asas, vel, x_motor,
                 x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                 x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                 c_ev,perfil_ev,p):
        con.Construtor2016.__init__(self,name,dz_asas, vel, x_motor,
                                    x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                                    x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                                    c_ev,perfil_ev,p)
        if p:
            self.plot()
            
        self.CL_cruzeiro = self.CL_cruzeiro()
        self.open_avl()
        self.initialize()
        self.config_case()
        self.execute_flow()
        self.eig_values()
        self.quit_avl()
        
    def CL_cruzeiro(self):
        """
        CL necessário para equilibrar o avião em velocidade de cruzeiro.
        """
        rho = self.ro_ar
        v = self.vel
        A = self.S_asaf
        g = self.g
        m = self.m_vazio+self.m_carga
        return 2*m*g/(rho*v**2*A)
        
    def open_avl(self):
        """
        Abre o AVL no python.
        """
        self.ps = sp.Popen(['avl.exe'],stdin=sp.PIPE,stdout=None,stderr=None)
    
    def quit_avl(self):
        ps = self.ps
        issueCmd(ps,'quit')
    
    def initialize(self):
        """
        Inicializa o número de casos especificados.
        """
        name = self.name
        ps = self.ps
        issueCmd(ps,'load Runs/%s.avl' %(name))
        issueCmd(ps,'mass Runs/%s.mass' %(name))
        issueCmd(ps,'mset 1')
        issueCmd(ps,'oper')
        for i in range(len(self.case_alphas)+self.number_trim_cases-1):
            issueCmd(ps,'+')
        
    def config_case(self):
        """
        Configura cada caso com o angulo de ataque especificado        
        """
        ps = self.ps
        vel = self.vel
        var = len(self.case_alphas)        
        
        for i in range(var):
            issueCmd(ps,'%d' %(i+1))
            issueCmd(ps,'a')
            issueCmd(ps,'a')
            issueCmd(ps,'%f' %(self.case_alphas[i]))
            issueCmd(ps,'m')
            issueCmd(ps,'v')
            issueCmd(ps,'%f' %(vel))
            issueCmd(ps,'')
            
        for i in range(var,self.number_trim_cases+var):
            issueCmd(ps,'%d' %(i+1))
            issueCmd(ps,'C1')
            issueCmd(ps,'C %f' %(self.CL_cruzeiro))
            issueCmd(ps,'B %f' %(10.0*(i-var)))
            issueCmd(ps,'')
                                    
    def execute_flow(self):                                                  
        """
        Analise aerodinâmica do avião gerado.
        Essa funçao executa o programa avl.exe
        """
        ps = self.ps
        var = len(self.case_alphas)
        name = self.name
        
        # Forças aerodinâmicas para os alfas especificados
        for i in range(var):
            issueCmd(ps,'%d' %(i+1))
            issueCmd(ps,'x')
            issueCmd(ps,'ft')
            issueCmd(ps,'Runs/%s_case%d.forces' %(name,i+1))
            if os.path.isfile('Runs/%s_case%d.forces' %(name,i+1)):
                issueCmd(ps,'O')
                
        # Casos gerados para a análise de autovalores            
        for i in range(var,self.number_trim_cases+var):
            issueCmd(ps,'%d' %(i+1))
            issueCmd(ps,'x')
            if i-var == 0:
                # Variáveis de estabilidade em angulo de trimagem desejado
                issueCmd(ps,'st')
                issueCmd(ps,'Runs/%s_trim.stability' %(name))
                if os.path.isfile('Runs/%s_trim.stability' %(name)):
                    issueCmd(ps,'O')       
        
        # Deleta os casos nao relacionados à trimagem
        for i in range(var):
            issueCmd(ps,'1')
            issueCmd(ps,'-')
        
        # Retorna para o menu principal
        issueCmd(ps,'')
    
    def eig_values(self):
        ps = self.ps
        name = self.name
        issueCmd(ps,'mode')
        issueCmd(ps,'0')
        issueCmd(ps,'n')
        issueCmd(ps,'w')
        issueCmd(ps,'Runs/%s_values.eig' %(name))
        if os.path.isfile('Runs/%s_values.eig' %(name)):
                    issueCmd(ps,'O')
                    
        # Retorna para o menu principal
        issueCmd(ps,'')
        
    @executarNaPasta('Graficos/Geometria')
    def plot(self):
        """
        Função de debug pra verificar em um gráfico 3d se os pontos foram
        gerados nas posições corretas.
        """
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        
        asaf, asat, ev = self.formato()[2:]        
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d([-1.25,1.25])
        ax.set_ylim3d([-1.25,1.25])
        ax.set_zlim3d([-1.25,1.25])
        ax.azim = -135
        
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
        
        for i in range(len(ev)):
            M = ev[i]
            x = np.append(M[:,0],M[0,0])
            y = np.append(M[:,1],M[0,1])
            z = np.append(M[:,2],M[0,2])
            ax.plot(x,y,z)
        
        plt.legend(loc='best')
        plt.savefig("geometria_%s.png" %(self.name), bbox_inches='tight', dpi=200)
        
if __name__ == '__main__':
    #Parametros gerais    
    name = 'A2016'
    dz_asas = 0.2
    vel = 15
    x_motor = -0.6
    
    #Asa frontal
    x_ba_asaf = -0.4
    c_asaf = 0.25
    ang_asaf = 4.0
    epsilon_asaf = 0.0
    perfilr_asaf = perfilp_asaf = 'S1223 MOD2015'
    
    # Asa traseira
    x_bf_asat = 0.3
    c_asat = 0.25
    ang_asat = 4.0
    epsilon_asat = 0.0
    perfilr_asat = perfilp_asat = 'S1223 MOD2015'
    
    # EV
    c_ev = 0.1
    perfil_ev = 'x'
    
    #Plotar gráficos
    plot = False
    
    aviao = Avaliador2016(name,dz_asas, vel, x_motor, 
                          x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                          x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                          c_ev,perfil_ev,plot)