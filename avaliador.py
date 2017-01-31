# -*- coding: utf-8 -*-
"""
Created on Sun Feb 08 12:32:47 2015

@author: Diego Chou
"""
print("\nCarregando modulos de 'Avaliador'...")
try:
    import numpy as np
    import subprocess as sp
    import os
    import psutil
    import time
    print("Modulos de 'Avaliador' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Avaliador'\n")
    raise

import construtor as con
import aerodinamica as aero
import estatica as est
import dinamica as din
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
    alfa_trim = 3.0
    case_alphas = [0.0,3.0,7.0,10.0] #Angulos de ataque de análise
    number_trim_cases = 1 # Números de caso em condição de trimagem
    
    config_m = [0.08,0.2] #Intervalo de valores aceitável para a margem estática
    
    def __init__(self,name,dz_asas,alfa,vel,x_motor,
                 x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                 x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                 c_ev,perfil_ev,p):
        con.Construtor2016.__init__(self,name,dz_asas, vel, x_motor,
                                    x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                                    x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                                    c_ev,perfil_ev,p)
        self.alfa = alfa
        self.CL_cruzeiro = self.CL_cruzeiro()
        
        # Execução do AVL 
        self.open_avl()
        self.initialize()
        self.config_case()
        self.execute_flow()
        self.eig_values()
        self.quit_avl()
        
        # Espera até que a execução do AVL esteja 
        # finalizada e o respectivo processo fechado.
        print('')
        inicial = time.clock()
        while True:
            if psutil.pid_exists(self.ps.pid):
                print("Executando AVL...")
                time.sleep(2)
            else:
                break
            if time.clock()-inicial > 60:
                self.ps.kill()
                raise ValueError("O AVL nao convergiu.")
                
        print("\nGeometria da aeronave:")
        print("         Posicao CG :           %s m" %(self.pos_cg))
        print("         Xcg :                  %f" %((self.pos_cg[0]-self.x_ba_asaf)/self.c_asaf))
        print("        ---------- ASA FRONTAL ---------")
        print("         X do bordo de ataque : %f m" %(self.x_ba_asaf))
        print("         Z do bordo de ataque : %f m" %(self.z_min))
        print("         Area :                 %f m^2" %(self.S_asaf))
        print("         Envergadura :          %f m" %(self.bw_asaf))
        print("         Corda :                %f m" %(self.c_asaf))
        print("         Angulo de incidencia : %f graus" %(self.ang_cr_asaf))
        print("         Angulo de torsao :     %f graus" %(self.epsilon_asaf))
        print("        --------- ASA TRASEIRA ---------")
        print("         X do bordo de ataque : %f m" %(self.x_bf_asat-self.c_asat))
        print("         Z do bordo de ataque : %f m" %(self.z_min+self.dz_asas))
        print("         Area :                 %f m^2" %(self.S_asat))
        print("         Envergadura :          %f m" %(self.bw_asat))
        print("         Corda :                %f m" %(self.c_asat))
        print("         Angulo de incidencia : %f graus" %(self.ang_cr_asat))
        print("         Angulo de torsao :     %f graus" %(self.epsilon_asat))
        print("        ------------ EV (x2) -----------")
        print("         Area :                 %f m^2" %(self.S_ev))
        print("         Envergadura :          %f m" %(self.b_ev))
        print("         Corda :                %f m" %(self.c_ev))
        print("        ------------ MOTOR -------------")
        print("         Posicao :              %s m" %(self.pos_motor))
        
        print("\nInercia da aeronave:")
        print("         Peso vazio :           %s kg" %(1.5*self.m_vazio))
        print("\nTensor de Inercia (kg*m^2) : \n%s" %(self.J))
            
        # Avaliacao de carga paga máxima
        args_aero = [self.name,self.alfa,self.case_alphas,self.m_vazio,\
                     self.S_asaf,self.bw_asaf,self.perfilr_asaf,self.perfilp_asaf,\
                     self.epsilon_asaf,self.vel,self.p]
        self.CPaga, self.alfa_estol,self.mac = aero.aerodinamica(*args_aero)
        
        # Avaliação da estabilidade estática
        args_estatica = [self.name,self.alfa_trim,self.alfa_estol,self.vel,\
                         self.pos_cg,self.config_m,self.mac,self.p]
        self.Xnp,self.fator_estatica = est.estabilidade_estatica(*args_estatica)
        
        try:
            self.fator_dinamica = din.estabilidade_dinamica(self.name,self.p)
        except ValueError, e:
            print("\n%s" %(e))
            self.fator_dinamica = 0.01
            
        print('\nFator de aerodinamica : %f' %(self.CPaga))
        print('Fator de estabilidade estatica : %f' %(self.fator_estatica))
        print('Fator de estabilidade dinamica : %f' %(self.fator_dinamica))

        # Pontuacao da aeronave
        self.pontuacao = self.CPaga*self.fator_estatica*self.fator_dinamica
        print("Pontuacao final : %f" %(self.pontuacao))
        
        if self.p: self.plot()
        
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
        self.ps = sp.Popen(['avl.exe'],stdin=sp.PIPE,stdout=open(os.devnull,'w'))
    
    def quit_avl(self):
        """
        Finaliza o processo do AVL para que ele não fique rodando em 2º plano.
        """
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
        Configura cada caso com a condicao especificada        
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
        xnp = self.Xnp
        ax.plot([cg[0]],[cg[1]],[cg[2]],'ro', label='CG')   
        ax.plot([xnp],[cg[1]],[cg[2]],'bo', label='$X_{np}$')
        
        objetos = [asaf,asat,ev]
        for i in objetos:
            for j in range(len(i)):
                M = i[j]
                x = np.append(M[:,0],M[0,0])
                y = np.append(M[:,1],M[0,1])
                z = np.append(M[:,2],M[0,2])
                ax.plot(x,y,z)   
        
        plt.legend(loc='best')
        plt.savefig("geometria_%s.png" %(self.name), bbox_inches='tight', dpi=200)
        
if __name__ == '__main__':

    #[0.0655481135642505, -0.6169047374302623, -0.19072996243401064, 0.3991774109065558, 
    # 2.603011485765066, -0.5885822268427869, 0.7139722884775602, 0.2885711556933519, 
    # 1.1756835144489113, 2.8652532575807257, 0.2815570703482143]
    
    #Parametros gerais    
    name = 'A2016'
    dz_asas = 0.0655481135642505
    alfa = 0.0
    vel = 20.0
    x_motor = -0.34
    
    #Asa frontal
    x_ba_asaf = -0.19072996243401064
    c_asaf = 0.3991774109065558
    ang_asaf = 2.603011485765066
    epsilon_asaf = -0.5885822268427869
    perfilr_asaf = perfilp_asaf = 'S1223 MOD2015'
    
    # Asa traseira
    x_bf_asat = 0.7139722884775602
    c_asat = 0.2885711556933519
    ang_asat = -1
    epsilon_asat = 2.8652532575807257
    perfilr_asat = perfilp_asat = 'NACA 0011'
    
    # EV
    c_ev = 0.20
    perfil_ev = 'NACA 0011'
    
    #Plotar gráficos
    plot = True
    
    aviao = Avaliador2016(name,dz_asas,alfa,vel, x_motor, 
                          x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                          x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                          c_ev,perfil_ev,plot)
