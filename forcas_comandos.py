# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:17:38 2016

@author: Usuario
"""

print("\nCarregando modulos de 'Forcas Comandos'...")

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import quad
    print("Modulos de 'Forcas Comandos' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Forcas Comandos'\n")
    raise

import apoio
import AeroPy.AeroPy as aeropy
import perfil

class Superficie_comando(perfil.Analise):
    """
    Superfície de comando simples (aileron, leme, profundor, flap, etc.)
    """
    def __init__(self,nome_perfil,alfa,corda,p_c,env,vel):
        """
        nome_perfil : nome do perfil, que deve estar listado na pasta Perfis/
        alfa : angulo de ataque do perfil        
        corda : corda do perfil em que se encontra a superfície de comando (usar valor médio)
        p_c : comprimento da superficie em porcentagem da corda (usar valor médio)
        env : envergadura da superfície de comando
        vel : velocidade de análise
        """
        self.ni = 15e-6
        self.rho = 1.086
        self.aerofolio = perfil.Analise(nome_perfil) 
        self.alfa = alfa
        self.p_c = p_c
        self.corda = corda
        self.env = env
        self.vel = vel
        
    def thin_airfoil_flap_cm(self,delta):
        """
        Calcula o coeficiente de momento na dobradiça da superfície móvel
        pela teoria do perfil fino
        
        delta : deflexao da superficie em graus
        """
        def B(n):
            k = 2.0/np.pi
            
            def f(theta):
                x = 0.5*((1.0-x_inicial)*(1.0+np.cos(theta))+x_inicial)
                return deta_dx(x)*np.cos(n*theta)

            integ = quad(f,0,np.pi)[0]
            return k*integ
            
        alfa = self.alfa*np.pi/180.0
        delta = delta*np.pi/180
        x_inicial = 1-self.p_c
        
        deta_dx = self.aerofolio.linha_media().derivative() #em função de x
    
        B0 = B(0)
        B1 = B(1)
        B2 = B(2)
        
        dCl_dalpha = 2*np.pi
        Cl = dCl_dalpha*(alfa+delta)-np.pi*(B0+B1)
        Cm_hinge = -0.25*(Cl+np.pi*(B1+B2))
        return -Cm_hinge # o Cm de atuacao é o oposto do aerodinamico

    def xfoil_flap_cm(self,delta):
        """
        Calcula o coeficiente de momento na dobradiça da superfície móvel
        pelo Xfoil, utilizando o AeroPy. Deve ser muito mais precisa que a 
        função anterior.
        
        delta : deflexao da superficie em graus
        """
        alpha = self.alfa
        x_hinge = 1-self.p_c
        x,y = self.aerofolio.coords_padrao(100) #.getPoints()
        delta = delta*np.pi/180
        
        @apoio.executarNaPasta("AeroPy")
        def cm_hinge():
            with apoio.esconderPrint():
                try:
                    return aeropy.calculate_flap_moment(x, y, alpha, x_hinge, delta,
                                                        unit_deflection = 'rad')
                except:
                    return np.nan
        return cm_hinge()

    def momento_servo(self,delta, option = 2):
        """
        Calcula o torque que um servo deve ter para segurar a superficie movel,
        a partir dos coeficientes aerodinamicos calculados.
        """
        rho = self.rho
        vel = self.vel
        env = self.env        
        
        if option == 1:
            Cm = self.thin_airfoil_flap_cm(delta)
            corda = self.corda*self.p_c
        elif option == 2:
            Cm = self.xfoil_flap_cm(delta)   
            corda = self.corda
        
        return Cm*0.5*rho*vel**2*corda**2*env
        
    def principal(self,delta_min,delta_max,name,n = 10):
        name = str(name)
        delta = np.linspace(delta_min,delta_max,n)
        xfoil_M = []
        
        for i in delta:
            xfoil_M.append(self.momento_servo(i,2))
            
        @apoio.executarNaPasta('Graficos/Aerodinamica')
        def grafico():
            plt.figure()
            plt.title(r'Torque requerido %s ($v=%.1f$ $m/s$)' %(name.upper(),vel))
            plt.grid('on')
            plt.plot(delta,xfoil_M,label = "XFoil")
            plt.legend(loc='best')
            plt.xlabel(r'$\delta$ ($graus$)')
            plt.ylabel(r'Torque ($Nm$)')
            plt.savefig("Torque_%s.png" %(name.lower()), bbox_inches='tight', dpi=200)
        
        grafico()

        xfoil_M = np.array(xfoil_M)
        xfoil_M = xfoil_M[np.logical_not(np.isnan(xfoil_M))]           
        
        return [min(xfoil_M),max(xfoil_M)]
        
if __name__ == "__main__":
    
    vel = 24.0
    flap = Superficie_comando("S1223 MOD2015",2.86,0.448,0.2,0.3084,vel)    
    aileron = Superficie_comando("S1223 MOD2015",2.86,0.3695,0.2436,0.480,vel)
    profundor = Superficie_comando("NACA 0011",1.0,0.2886,0.40,1.1301,vel)
    
#    print aileron.momento_servo(15,2)    
    
    fl = flap.principal(0,30,"flap")
    ail = aileron.principal(-15,15, "aileron")
    prof = profundor.principal(-20,20, "profundor")
#    
    print("\nTorques requeridos a %.1f m/s : [min, max]" %(vel))
    print("Aileron :   %s Nm" %(ail))
    print("Profundor : %s Nm" %(prof))
    print("Flap :      %s Nm" %(fl))      
    
