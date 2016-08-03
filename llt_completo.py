# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 15:10:21 2016

@author: Diego Chou
"""

print("\nCarregando modulos de 'LLT completo'...")

try:
    import numpy as np
    import matplotlib.pyplot as plt
    print("Modulos de 'LLT' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'LLT completo'\n")
    raise

import apoio
import perfil
import llt

class S_sustentadora_completa(llt.S_sustentadora):
    """
    Essa classe exerce todas as funçoes da classe S_sustentadora,
    com a adição do efeito da fuselagem e flaps na parte retangular da asa.
    """
    def __init__(self, larg_fus, a_fus, cd0_fus,
                 perfil_raiz, perfil_ponta, perfil_flap,
                 alfa, Sw, bw, afil, graus_torcao, vel, n, otimizar = True, p = False):
        self._larg_fus = larg_fus 
        
        #Os valores abaixo devem ser adimensionalizados pela corda raiz da asa
        self.a_fus = a_fus*180/np.pi
        self.cd0_fus = cd0_fus
        
        self._perfil_flap = perfil.Analise(perfil_flap)
        
        # Se a asa apresenta afilamento, o flap se estenderá sobre a seçao retangular da asa
        # Caso contrário, ele terá o comprimento de 1/3 da semi-envergadura
        if afil == 1.0 or otimizar == False:
            self.l_flap = bw/6.0
            
        llt.S_sustentadora.__init__(self, perfil_raiz, perfil_ponta, alfa, Sw, bw, afil,
                                    graus_torcao, vel, n, otimizar, p)
                            
    def _clmax_local(self):
        n = self._n        
        y = self._y
        
        if self.trans != 0.0:
            trans = self.trans
        else: 
            trans = self.l_flap
            
        clmax = np.zeros(n)
        Re_asa = self._Re(self._cordas)
        
        if self._perfil_raiz.name != self._perfil_ponta.name and self._perfil_raiz.name != self._perfil_flap.name:
            for i in range(n):
                if y[i] < trans:
                    flap = self._perfil_flap.interpReMach_Clmax(Re_asa[i], self._mach)[0]
                    clmax[i] = flap
                else:
                    raiz = self._perfil_raiz.interpReMach_Clmax(Re_asa[i], self._mach)[0]
                    ponta = self._perfil_ponta.interpReMach_Clmax(Re_asa[i], self._mach)[0]
                    clmax[i] = raiz + (ponta-raiz)*(y[i]-trans)/(0.5*self.bw-trans)
                    
        elif self._perfil_raiz.name != self._perfil_flap.name:
            for i in range(n):
                if y[i] < trans:
                    clmax[i] = self._perfil_flap.interpReMach_Clmax(Re_asa[i], self._mach)[0]
                else:                    
                    clmax[i] = self._perfil_raiz.interpReMach_Clmax(Re_asa[i], self._mach)[0]        
                    
        else:
            for i in range(n):
                clmax[i] = self._perfil_raiz.interpReMach_Clmax(Re_asa[i], self._mach)[0]  
            
        return clmax
        
    def _mod_llt(self,alpha, _iter = False, _dist = False, _dw = False):
        """
        Cálculo da Linha Sustentadora de Prandtl modificada
        """
        def circulacao():
            """
            Calculo da circulação sobre a asa (normalizada pela velocidade)
            e angulo induzido.
            """
            gamma = np.zeros(n)
            alpha_i = np.zeros(n)
            for i in range(n):
                for j in range(n):
                    gamma[i]+=2*A[j]*np.sin((2*j+1)*theta[i])
                    alpha_i[i]+=(2*j+1)*A[j]*np.sin((2*j+1)*theta[i])/(np.sin(theta[i]))
            gamma *= span
            alpha_i = np.tan(alpha_i/self.vel)
            return gamma,alpha_i
            
        def coeffs():
            """
            Calculo dos coeficientes aerodinâmicos conforme a teoria da linha 
            sustentadora modificada
            """
            ########### Computing CL ###########
            CL = A[0]*np.pi*self.AR
            ###### End of CL computation #######
            
            ########### Computing CD0 ##########
            mach = self._mach
            Sw = self.Sw
            
            dy = self._dy
            area_s = (c[1:]+c[:-1])*dy/2.0
            dim = len(area_s) 
            CD0s = np.zeros(dim)
            alpha_secao = 0.5*(value_alpha_local[:-1]+value_alpha_local[1:])
            y_secao = 0.5*(y[:-1]+y[1:])
            
            if self._perfil_raiz.name != self._perfil_ponta.name and self._perfil_raiz.name != self._perfil_flap.name:
                for i in range(dim):
                    if y_secao[i] <= l_fus:
                        CD0s[i] = self.cd0_fus
                    elif y_secao[i] > l_fus and y_secao[i] < trans:
                        CD0s[i] = self._perfil_flap.interpReMachAlfa(alpha_secao[i], Re_asa[i], mach, 'Cd')  
                    else:
                        secao_raiz = self._perfil_raiz.interpReMachAlfa(alpha_secao[i], Re_asa[i], mach, 'Cd')
                        secao_ponta = self._perfil_ponta.interpReMachAlfa(alpha_secao[i], Re_asa[i], mach, 'Cd')
                        CD0s[i] = secao_raiz+(secao_ponta-secao_raiz)*(y_secao[i]-trans)/(0.5*span-trans)
                        
            elif self._perfil_raiz.name != self._perfil_flap.name:
                for i in range(dim):
                    if y_secao[i] <= l_fus:
                        CD0s[i] = self.cd0_fus
                    elif y_secao[i] > l_fus and y_secao[i] < trans:
                        CD0s[i] = self._perfil_flap.interpReMachAlfa(alpha_secao[i], Re_asa[i], mach, 'Cd')
                    else:                    
                        CD0s[i] = self._perfil_raiz.interpReMachAlfa(alpha_secao[i], Re_asa[i], mach, 'Cd')
                        
            else:
                for i in range(dim):
                    if y_secao[i] <= l_fus:
                        CD0s[i] = self.cd0_fus
                    else:
                        CD0s[i] = self._perfil_raiz.interpReMachAlfa(alpha_secao[i], Re_asa[i], mach, 'Cd')
                
            CD0s_x = CD0s*np.cos(alpha_secao)
            CD0s_y = CD0s*np.sin(alpha_secao)
            
            CD0 = 2*sum(CD0s_x*area_s)/Sw
            CL_CD0 = 2*sum(CD0s_y*area_s)/Sw
            ###### End of CD0 computation ######
            
            ########### Computing CDi ##########
            CDi = 0.0
            if A[0] != 0.0:
                for i in range(n):
                    CDi += (2.0*i+1)*A[i]**2
                CDi *= np.pi*self.AR                    
            ###### End of CDi computation ######

            ############ Miscelanea ############                        
            # Como CD0s tem uma dimensão a menos do que y por contabilizar
            # seções e não os pontos, é necessário fazer uma média para
            # a seção que falta para que a função retorne a distribuição            
            CD0s_x = 0.5*(np.append(CD0s_x,CD0s_x[-1])+np.append(CD0s_x[0],CD0s_x))
            CD0s_y = 0.5*(np.append(CD0s_y,CD0s_y[-1])+np.append(CD0s_y[0],CD0s_y))
            
            # O CL que a função retorna já conta com a contribuição do CDO
            # no eixo y            
            CL -= CL_CD0            
            ############ Miscelanea ############ 
            
            return CD0s_x,CD0s_y,CL,CD0,CDi
        
        if not _iter and not _dist and not _dw:
            print("\nIniciando calculo da Linha Sustentadora")        
        
        alpha = np.radians(alpha)
        y = self._y
        c = self._cordas
        
        if self.trans != 0.0:
            trans = self.trans
        else: 
            trans = self.l_flap
            
        l_fus = self._larg_fus/2.0 
        a_fus = self.a_fus
        n = self._n
        theta = self._theta
        span = self.bw
        
        Re_asa = self._Re(self._cordas)
        
#        if self._p and (not _iter) and (not _dist) and (not _dw):
#            p = True
#        else:
#            p = False
        
        #Referente à raiz da asa
        thetaroot = 0.0
        
        # Referente à ponta da asa
        thetatip = np.radians(self.epsilon)
        
        #ângulos geométricos e dCl/dAlpha ao longo da asa
        #A torção, caso haja, só é aplicada após a transição
        alp = np.zeros(n)
        alpha_geo = np.zeros(n)
        a = np.zeros(n)
        
        if self._perfil_raiz.name != self._perfil_ponta.name and self._perfil_raiz.name != self._perfil_flap.name:
            for i in range(n):
                if y[i] <= l_fus:
                    alp[i] = alpha_geo[i] = alpha+thetaroot
                    a[i] = a_fus
                elif y[i] > l_fus and y[i] < trans:
                    #Referente ao flap da asa
                    dClflap = self._perfil_flap.ajustelinear(Re_asa[i],self._mach).get('dCl')
                    a0flap = dClflap[0]*180/np.pi
                    alpha0flap = np.radians(-dClflap[1]/dClflap[0])
                    
                    alp[i] = alpha+thetaroot-alpha0flap
                    alpha_geo[i] = alpha+thetaroot
                    a[i] = a0flap
                else:
                    #Referente à raiz da asa
                    dClroot = self._perfil_raiz.ajustelinear(Re_asa[i],self._mach).get('dCl')
                    a0root = dClroot[0]*180/np.pi
                    alpha0root = np.radians(-dClroot[1]/dClroot[0])
                    
                    # Referente à ponta da asa
                    dCltip = self._perfil_ponta.ajustelinear(Re_asa[i],self._mach).get('dCl')
                    a0tip = dCltip[0]*180/np.pi
                    alpha0tip = np.radians(-dCltip[1]/dCltip[0])
                
                    alp[i] = alpha+thetaroot-alpha0root+(thetatip-alpha0tip-thetaroot+alpha0root)*(y[i]-trans)/(0.5*span-trans)
                    alpha_geo[i] = alpha+thetaroot+(thetatip-thetaroot)*(y[i]-trans)/(0.5*span-trans)                
                    a[i] = a0root+(a0tip - a0root)*(y[i]-trans)/(0.5*span-trans)
        
        elif self._perfil_raiz.name != self._perfil_flap.name:            
            for i in range(n):
                if y[i] <= l_fus:
                    alp[i] = alpha_geo[i] = alpha+thetaroot
                    a[i] = a_fus
                elif y[i] > l_fus and y[i] < trans:
                    #Referente ao flap da asa
                    dClflap = self._perfil_flap.ajustelinear(Re_asa[i],self._mach).get('dCl')
                    a0flap = dClflap[0]*180/np.pi
                    alpha0flap = np.radians(-dClflap[1]/dClflap[0])
                    
                    alpha_geo[i] = alpha+thetaroot
                    alp[i] = alpha_geo[i]-alpha0flap
                    a[i] = a0flap
                else:
                    dCl = self._perfil_raiz.ajustelinear(Re_asa[i],self._mach).get('dCl')
                    a0 = dCl[0]*180/np.pi
                    alpha0 = np.radians(-dCl[1]/dCl[0])
                    
                    alpha_geo[i] = alpha+thetaroot+(thetatip-thetaroot)*(y[i]-trans)/(0.5*span-trans) 
                    alp[i] = alpha_geo[i]-alpha0
                    a[i] = a0
                    
        else:
            for i in range(n):
                dCl = self._perfil_raiz.ajustelinear(Re_asa[i],self._mach).get('dCl')
                a0 = dCl[0]*180/np.pi
                alpha0 = np.radians(-dCl[1]/dCl[0])
                
                if y[i] <= l_fus:
                    alp[i] = alpha_geo[i] = alpha+thetaroot
                    a[i] = a_fus
                elif y[i] > l_fus and y[i] < trans:
                    alpha_geo[i] = alpha+thetaroot
                    alp[i] = alpha_geo[i]-alpha0
                    a[i] = a0
                else:
                    alpha_geo[i] = alpha+thetaroot+(thetatip-thetaroot)*(y[i]-trans)/(0.5*span-trans)
                    alp[i] = alpha_geo[i]-alpha0
                    a[i] = a0
        
        #Lado direito (Right Hand Side) da equação 
        # Set up 2n x 2n system of equations for A1, A3 , ... A2n-1
        rhs = alp*np.sin(theta)*c*a/(4.0*span)
        
        b = np.zeros((n,n))
        mu = c*a/(4.*span)
        
        for i in range(n):
            l = float(2*i+1)
            b[:,i]=np.sin(l*theta)*(mu*l+np.sin(theta))
        
        if not _iter and not _dist and not _dw:
            print("         Resolvendo %d coeficientes de Fourier" %(n))
        # Solve for the Fourier Coefficients
        A = np.linalg.solve(b,rhs)
        
        if not _iter and not _dist and not _dw:
            print("         Calculando circulacao sobre a asa")
        #Calculos de circulação e angulo induzido
        value_gamma,value_alpha_i = circulacao()
    
        value_alpha_local = alpha_geo-value_alpha_i
        
        if not _iter and not _dist and not _dw:
            print("         Calculando coeficiente aerodinâmicos")
        # Distribuição de cl
        dist_cl = 2.0*value_gamma/self._cordas

        #Distribuições de CD0 projetadas em cada eixo e valores de CL, CD0 e CDi globais
        CD0s_x,CD0s_y,value_CL,value_CD0,value_CDi = coeffs()
        
        if not _iter and not _dist and not _dw:
            print("         Alfa = %.1f graus" %(self.alfa))
            print("         CL =   %.3f" %(value_CL))
            print("         CD =   %.3f" %(value_CD0+value_CDi))
            print("         CD0 =  %.3f" %(value_CD0))
            print("         CDi =  %.3f" %(value_CDi))

        if _dw:
            return value_gamma
            
        elif _dist:
            # Distribuições
            value_cl = dist_cl-CD0s_y #Contabilizada a contribuição do CD0 no eixo y            
            value_cdi = np.sin(value_alpha_i)*dist_cl
            value_cd0 = CD0s_x
            return value_gamma,value_alpha_i,value_alpha_local,\
                   value_cl,value_cd0,value_cdi,\
                   value_CL,value_CD0,value_CDi
        
        elif _iter:
            return value_CL,value_CD0,value_CDi
        
        else:
            value_cl = dist_cl-CD0s_y #Contabilizada a contribuição do CD0 no eixo y 
            return value_alpha_local,value_CL,value_CD0,value_CDi
            
    def _CM(self):
        """
        Funçao que estima o coeficiente de momento da asa por integração dos
        CM's dos perfis.
        """
        l_fus = self._larg_fus/2.0
        dy = self._dy
        area_secao = (self._cordas[1:]+self._cordas[:-1])*dy/2.0
        n = len(area_secao)
        alpha_secao = 0.5*(self._alpha_local[:-1]+self._alpha_local[1:])
        y_secao = 0.5*(self._y[:-1]+self._y[1:])
        Re_asa = self._Re(self._cordas)
        
        if self.trans != 0.0:
            trans = self.trans
        else: 
            trans = self.l_flap
            
        mach = self._mach
        ca = self._p_c
        Cm_secao = np.zeros(n)  
        corda_secao = 0.5*(self._cordas[:-1]+self._cordas[1:])
        
        if self._perfil_raiz.name != self._perfil_ponta.name and self._perfil_raiz.name != self._perfil_flap.name:
            for i in range(n):
                ang = alpha_secao[i]
                
                if y_secao[i] <= l_fus:
                    Cl = Cd = xCp = 0.0
                
                elif y_secao[i] > l_fus and y_secao[i] < trans:
                    Cl = self._perfil_flap.interpReMachAlfa(ang,Re_asa[i],mach,'Cl')
                    Cd = self._perfil_flap.interpReMachAlfa(ang,Re_asa[i],mach,'Cd')
                    xCp = self._perfil_flap.interpXcp(ang,Re_asa[i],mach)
                     
                else:
                    Cl_raiz = self._perfil_raiz.interpReMachAlfa(ang,Re_asa[i],mach,'Cl')
                    Cl_ponta = self._perfil_ponta.interpReMachAlfa(ang,Re_asa[i],mach,'Cl')
                    Cd_raiz = self._perfil_raiz.interpReMachAlfa(ang,Re_asa[i],mach,'Cd')
                    Cd_ponta = self._perfil_ponta.interpReMachAlfa(ang,Re_asa[i],mach,'Cd')
                    xCp_raiz = self._perfil_raiz.interpXcp(ang,Re_asa[i],mach)
                    xCp_ponta = self._perfil_ponta.interpXcp(ang,Re_asa[i],mach)
                    
                    Cl = Cl_raiz + (Cl_ponta-Cl_raiz)*(y_secao[i]-trans)/(0.5*self.bw-trans)
                    Cd = Cd_raiz + (Cd_ponta-Cd_raiz)*(y_secao[i]-trans)/(0.5*self.bw-trans)
                    xCp = xCp_raiz + (xCp_ponta-xCp_raiz)*(y_secao[i]-trans)/(0.5*self.bw-trans)
                                 
                Cm_secao[i] = (ca[i]-xCp)*(Cl*np.cos(ang)-Cd*np.sin(ang))
                
        elif self._perfil_raiz.name != self._perfil_flap.name:
            for i in range(n):
                ang = alpha_secao[i]
                
                if y_secao[i] <= l_fus:
                    Cl = Cd = xCp = 0.0
                
                elif y_secao[i] > l_fus and y_secao[i] < trans:
                    Cl = self._perfil_flap.interpReMachAlfa(ang,Re_asa[i],mach,'Cl')
                    Cd = self._perfil_flap.interpReMachAlfa(ang,Re_asa[i],mach,'Cd')
                    xCp = self._perfil_flap.interpXcp(ang,Re_asa[i],mach)
                     
                else:
                    Cl = self._perfil_raiz.interpReMachAlfa(ang,Re_asa[i],mach,'Cl')
                    Cd = self._perfil_raiz.interpReMachAlfa(ang,Re_asa[i],mach,'Cd')
                    xCp = self._perfil_raiz.interpXcp(ang,Re_asa[i],mach)
                                 
                Cm_secao[i] = (ca[i]-xCp)*(Cl*np.cos(ang)-Cd*np.sin(ang))
                
        else:
            for i in range(n):
                ang = alpha_secao[i]
                if y_secao[i] <= l_fus:
                    Cl = Cd = xCp = 0.0
                else:
                    Cl = self._perfil_raiz.interpReMachAlfa(ang,Re_asa[i],mach,'Cl')
                    Cd = self._perfil_raiz.interpReMachAlfa(ang,Re_asa[i],mach,'Cd')
                    xCp = self._perfil_raiz.interpXcp(ang,Re_asa[i],mach)
            
                Cm_secao[i] = (ca[i]-xCp)*(Cl*np.cos(ang)-Cd*np.sin(ang))
        
        def cm_local():
            plt.figure()
            plt.grid('on')
            plt.title(r'Coeficiente de momento local')
            plt.xlabel(r'Semi-envergadura normalizada ($2y/b$)')
            plt.ylabel(r'$C_m$')
            plt.plot(2.0*y_secao/self.bw,Cm_secao,label=r'$\alpha = %.1f^o$' %(self.alfa))
            plt.ylim(-0.4,0.0)
            plt.legend(loc='best')
            plt.savefig("distCm.png", bbox_inches='tight', dpi=200)
            plt.close()

        def momento_torsor_local():
            plt.figure()
            plt.grid('on')
            plt.title(r'Momento torsor local')
            plt.xlabel(r'Semi-envergadura ($m$)')
            plt.ylabel(r'$M_t$ ($N$)')
            plt.plot(y_secao,0.5*self._ar.rho()*self.vel**2*corda_secao*Cm_secao,label=r'$\alpha = %.1f^o$' %(self.alfa))
            plt.legend(loc='best')
            plt.savefig("distMt.png", bbox_inches='tight', dpi=200)
            plt.close()
        
        if self._p:
            apoio.executarNaPasta(self._pasta_graficos)(cm_local)()
            apoio.executarNaPasta(self._pasta_graficos)(momento_torsor_local)()
        
        CM = 2*sum(Cm_secao*area_secao)/self.Sw
        return CM
        
if __name__ == "__main__":
    Asa1 = S_sustentadora_completa(0.12,0.006,0.0088, #Fuselagem
                                   "S1223 MOD2015", "S1223 MOD2015", "S1223 MOD2015",
                                   22.5, 0.797, 1.997, 0.65,-3.0, 13.7, 150, p = False)
    print Asa1.de_da((1-0.23)*0.190730+1.25*0.425401,0.065)
