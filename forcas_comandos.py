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
        self.rho = 1.086
        self.aerofolio = perfil.Analise(nome_perfil) 
        self.alfa = alfa*np.pi/180.0
        self.p_c = p_c
        self.corda = corda
        self.env = env
        self.vel = vel
        
        self.linha_media = self.aerofolio.linha_media()
    
    def thin_airfoil_coeffs(self,delta):
        """
        Calcula a linha média da superfície móvel 
        delta : deflexao da superficie em graus
        """
        def B(n):
            k = 2.0/np.pi
            
            def f(theta):
                x = 0.5*((1.0-x_inicial)*(1.0+np.cos(theta))+x_inicial)
                return deta_dx(x)*np.cos(n*theta)

            integ = quad(f,0,np.pi)[0]
            return k*integ
            
        alfa = self.alfa
        delta = delta*np.pi/180
        x_inicial = 1-self.p_c
        
        deta_dx = self.linha_media.derivative() #em função de x
    
        B0 = B(0)
        B1 = B(1)
        B2 = B(2)
        
        dCl_dalpha = 2*np.pi
        Cl = dCl_dalpha*(alfa+delta)-np.pi*(B0+B1)
        Cm0 = -0.25*(Cl+np.pi*(B1+B2))
        print Cm0
        return Cl, Cm0
        
    def forces(self,delta):
        Cl, Cm0 = self.thin_airfoil_coeffs(delta)
        rho = self.rho
        vel = self.vel
        corda = self.corda*self.p_c
        env = self.env
        L = Cl*0.5*rho*vel**2*corda*env
        M = Cm0*0.5*rho*vel**2*corda**2*env
        return L, M 
        
if __name__ == "__main__":
    vel = 20.0    
    flap = Superficie_comando("S1223 MOD2015",2.86,0.448,0.2,0.3084,vel)
    M_flap = flap.forces(30.0)[1]
    
    aileron = Superficie_comando("S1223 MOD2015",2.86,0.3695,0.2436,0.480,vel)
    M_aileron = aileron.forces(15.0)[1]
    
    profundor = Superficie_comando("NACA 0011",1.0,0.2886,0.40,1.1301,vel)
    M_profundor = profundor.forces(20.0)[1]
    
    print("\nMomento flap :      %.3f" %(M_flap))
    print("Momento aileron :   %.3f" %(M_aileron))
    print("Momento profundor : %.3f" %(M_profundor))
        
    