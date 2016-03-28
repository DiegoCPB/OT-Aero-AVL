# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 16:41:01 2014

@author: Diego Chou
"""

print("\nCarregando modulos de 'LLT'...")

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import fmin
    from scipy.interpolate import UnivariateSpline
    from scipy.integrate import quad
    print("Modulos de 'LLT' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'LLT'\n")
    raise

import apoio
import geral
import perfil

"""
 As análises aqui feitas têm como referência a teoria da linha sustentadora de Prandtl,
 modificada para contabilizar também o arrasto parasita da asa.
 
 No caso em que há perfis diferentes entre a ponta e a raiz, tal mudança só ocorre
 a partir da transição reto-trapezoidal.
 
 Da mesma forma, caso haja torção (washout), esta só é aplicada a partir da transição 
 reto-trapezoidal.
 
 Todos os pesos estão em kg. 
 Todas as distâncias em metros.
 Todos as funções retornam angulos em graus.
 Porque eu sou o programador e eu faço o que eu quiser!     
"""

class S_sustentadora(object):
    """
    Essa classe define uma superficie sustentadora geral, bem como
    os metodos para ESTIMAR suas propriedades aerodinamicas a partir 
    de seu formato.
    
    Atributos geometricos:
    ---------------------
    * vel - velocidade de analise
    * alfa - angulo de ataque
    * Sw - area de superficie
    * bw - envergadura
    * afil - afilamento
    * epsilon - torcao geometrica
    * AR - razao de aspecto (alongamento)
    * trans - distancia da transicao reto-trapezoidal da raiz
    * raiz - comprimento da corda na raiz
    * ponta - comprimento da corda na ponta
    * mac - comprimento da corda media aerodinamica
    * pos_mac - distancia da MAC da raiz
    * offset - offset da ponta da asa em relacao a raiz
    
    Atributos aerodinamicos (em relacao a area da asa):
    --------------------------------------------------
    * CL - coeficiente de sustentacao
    * CD0 - coeficiente de arrasto parasita
    * CDi - coeficiente de arrasto induzido
    * CM - coeficiente de momento
    * alfa_CL0 - angulo de sustentacao nula
    * alfa_estol - angulo de estol
    * CLmax - CL maximo
    * funcao_CL - polar de sustentacao da asa
    * funcao_CD - polar de arrasto da asa
    * funcao_CDi - polar de arrasto induzido da asa
    * de_da - derivada do downwash em funcao do angulo de ataque 
    """
    def __init__(self, perfil_raiz, perfil_ponta, alfa, Sw, bw, afil,
                 graus_torcao, vel, n, otimizar = True, p = False):
        #Parâmetros gerais da classe
        self._pasta_graficos = 'Graficos/Aerodinamica/%s' %(type(self).__name__) 
        self._p = p
        self._ar = geral.Ar() # altitude padrão da classe Ar
        self.vel = vel # condições de cálculo
        self._mach = vel/geral.vel_som
        
        # Perfis aerodinâmicos
        self._perfil_raiz = perfil.Analise(perfil_raiz)
        self._p_c_raiz = self._perfil_raiz.ac_medio(self._mach)
        self._perfil_ponta = perfil.Analise(perfil_ponta)
        self._p_c_ponta = self._perfil_ponta.ac_medio(self._mach)
        
        #Parâmetros geométricos da asa
        self.alfa = alfa
        self.Sw = Sw 
        self.bw = bw
        
        if afil > 1.0:
            raise ValueError("O afilamento maior que 1.")
        
        self.afil = afil
        self.epsilon = graus_torcao
        self.AR = self.bw**2/self.Sw

        if self.AR < 4.0:
            print("\n######################################################")
            print("# ATENCAO: AR = %.2f < 4.                            #" %(self.AR)) 
            print("#        O LLT retorna valores nao confiaveis.       #")
            print("#        Reveja a geometria do %s.                   #" %(type(self).__name__))
            print("######################################################")
        
        # Discretização da asa
        self._n = n # nº de pontos discretos 
        self._x,self._dx = np.linspace(0, 0.5*self.bw, num=self._n,retstep = True) #linear
        self._theta = np.linspace(1,self._n,self._n)*np.pi/(2.0*self._n) #para o LLT
        self._y = 0.5*self.bw*np.cos(self._theta) #para o LLT
        
        # Duas formas de calcular dy. Se n for muito grande, descomente a de baixo
        self._dy = self._y[:-1]-self._y[1:]
#        self._dy = 0.5*self.bw*np.sin(0.5*(self._theta[1:]+self._theta[:-1]))*np.pi/(2.0*self._n) #módulo da derivada de y
        
        #Parâmetros da otimização
        if otimizar and self.afil != 1.0:
            self.trans = self._otimo()
        else:
            self.trans = 0.0
        self._p_c = self._p_c()
        self._cordas = self._planta(self.trans,llt=True) #para o LLT
        self.raiz, self.ponta = self._cordas[0::len(self._cordas)-1][::-1]
        self.mac, self.pos_mac, self._Re_mac = self._mac()
        self.offset = self._offset(self._cordas[-1],self._cordas[0])
        
        #Coeficientes Aerodinâmicos:
        self._alpha_local,self.CL,self.CD0,self.CDi = self._mod_llt(self.alfa)
        
        #Estol
        self.alfa_estol,self.CLmax = self._estol()
        
        # Polares
        self.alfa_CL0,self.funcao_CL,\
        self.funcao_CD,self.funcao_CDi = self._polares(np.linspace(0,self.alfa_estol,4))
       
        #Coeficiente de momento
        self.CM = self._CM()
        
        if p:
            apoio.executarNaPasta(self._pasta_graficos)(self._superior)()
            self._distribuicoes([self.alfa])
        
    def _Re(self, corda):
        """
        Cálculo do número de Reynolds para uma corda específica
        """
        rho = self._ar.rho()
        mi = self._ar.mi
        vel = self.vel
        Re = corda*rho*vel/mi
        return Re
    
    def _elipt(self,llt = False):
        """
        Asa elíptica
        """
        if llt:
            c = self._y
        else:
            c = self._x
        Sw = self.Sw
        bw = self.bw
        # Vide fórmula da elipse 
        c_raiz = 4*Sw/(np.pi*bw)
        
        if llt:        
            elipt = c_raiz*np.sin(self._theta)
        else:
            elipt = 2*c_raiz*np.sqrt(0.25-(c/bw)**2)
        return elipt

    def _planta(self,trans,llt = False):
        """
        Asa em planta
        """
        if llt:
            c = self._y
        else:
            c = self._x
        Sw = self.Sw
        bw = self.bw
        afil = self.afil
        c_raiz = -Sw/(2*((afil+1)*(trans/2-bw/4)-trans))       
        c_ponta = afil*c_raiz        
        planta = []
        
        for i in c:
            if i <= trans:
                planta.append(c_raiz)
#                print c_raiz
            else:
                cx = cx = c_raiz + (c_ponta-c_raiz)*(i-trans)/(0.5*bw-trans)
                planta.append(cx)
#                print cx  
        return np.array(planta)
    
    def _erro(self, arg):
        """
        Erro geométrico da asa em planta em relação à asa elíptica
        """
        arg = float(arg)
        a = self._elipt()
        b = self._planta(arg)
        
        dist = 0
        for i in range(len(self._y)):
            dist += (a[i]-b[i])**2  # Distância vetorial quadrática
        
        return dist #sum((a-b)**2)
        
    def _otimo(self):
        """
        Cálculo da posição ótima da transição reto-trapezoidal
        """
        x0 = 0.3  # 0.3 é o chute inicial
        print("\nCalculando erro relativo total minimo do formato da asa...\n")  
        ot = fmin(self._erro, x0) # Otimizaçao do scipy pelo método de Nelder-Mead (Simplex)  
#        print trans_ot        
        return ot[0] # O numpy retorna um array de 1 dimensão     
        
    def _offset(self,raiz, ponta):
        """
        Offset da ponta da asa em relação à raiz
        """
        p_c_raiz = self._p_c_raiz
        p_c_ponta = self._p_c_ponta
        offset = p_c_raiz*raiz-p_c_ponta*ponta
        return offset
        
    def _mac(self):
        """
        Cálculo do comprimento, posição e número de reynolds 
        da corda média aerodinamica (MAC)
        """
        planta = self._planta(self.trans)
        
        if self.afil == 1.0:
            mac = planta[0]
            pos_mac = 0.0
        else:
            Sw = self.Sw
            x = self._x
            dx = self._dx
            mac = (2/Sw)*sum((0.5*(planta[1:]+planta[:-1]))**2*dx)

            corda_antes, corda_depois = apoio.intervalo(planta,mac)
            pos_antes = x[np.nonzero(planta==corda_antes)[0][0]]
            pos_depois = x[np.nonzero(planta==corda_depois)[0][0]]
                                   
            pos_mac = apoio.interplinear(corda_antes, pos_antes, corda_depois, pos_depois, mac)
        
        Re_mac = self._Re(mac)
        return mac, pos_mac, Re_mac
        
    def _p_c(self):
        y = self._y
        n = self._n
        trans = self.trans
        pcr = self._p_c_raiz
        pcp = self._p_c_ponta
        p_c = np.zeros(n)
        
        for i in range(n):        
            if y[i] < trans:
                p_c[i] = pcr
            else:
                p_c[i] = pcr + (pcp-pcr)*(y[i]-trans)/(0.5*self.bw-trans)
        
        return p_c
        
    def _clmax_local(self):
        n = self._n        
        y = self._y
        trans = self.trans
        clmax = np.zeros(n)
        
        if self._perfil_raiz.name != self._perfil_ponta.name:
            raiz = self._perfil_raiz.interpReMach_Clmax(self._Re(self._cordas[-1]), self._mach)[0]
            ponta = self._perfil_ponta.interpReMach_Clmax(self._Re(self._cordas[0]), self._mach)[0]
            for i in range(n):
                if y[i] < trans:
                    clmax[i] = raiz
                else:
                    clmax[i] = raiz + (ponta-raiz)*(y[i]-trans)/(0.5*self.bw-trans)
                    
        else:
            Re_asa = self._Re(self._cordas)
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
            
            if self._perfil_raiz.name != self._perfil_ponta.name:
                for i in range(dim):
                    if y_secao[i] <= trans:
                        CD0s[i] = self._perfil_raiz.interpReMachAlfa(alpha_secao[i], Re_raiz, mach, 'Cd')  
                    else:
                        secao_raiz = self._perfil_raiz.interpReMachAlfa(alpha_secao[i], Re_raiz, mach, 'Cd')
                        secao_ponta = self._perfil_ponta.interpReMachAlfa(alpha_secao[i], Re_ponta, mach, 'Cd')
                        CD0s[i] = secao_raiz+(secao_ponta-secao_raiz)*(y_secao[i]-trans)/(0.5*span-trans)
                        
            else:
                Re_asa = self._Re(self._cordas)
                for i in range(dim):
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
        trans = self.trans
        n = self._n
        theta = self._theta
        span = self.bw
        
        if self._p and (not _iter) and (not _dist) and (not _dw):
            p = True
        else:
            p = False
        
        #Referente à raiz da asa
        croot = c[-1]    
        thetaroot = 0.0
        
        # Referente à ponta da asa
        ctip = c[0]
        thetatip = np.radians(self.epsilon)
        
        #ângulos geométricos e dCl/dAlpha ao longo da asa
        #A torção, caso haja, só é aplicada após a transição
        alp = np.zeros(n)
        alpha_geo = np.zeros(n)
        a = np.zeros(n)
        
        if self._perfil_raiz.name != self._perfil_ponta.name:
            #Referente à raiz da asa
            Re_raiz = self._Re(croot) 
            dClroot = self._perfil_raiz.ajustelinear(Re_raiz,self._mach,p).get('dCl')
            a0root = dClroot[0]*180/np.pi
            alpha0root = np.radians(-dClroot[1]/dClroot[0])
            
            # Referente à ponta da asa
            Re_ponta = self._Re(ctip)
            dCltip = self._perfil_ponta.ajustelinear(Re_ponta,self._mach,p).get('dCl')
            a0tip = dCltip[0]*180/np.pi
            alpha0tip = np.radians(-dCltip[1]/dCltip[0])
            
            for i in range(n):
                if y[i] < trans:
                    alp[i] = alpha+thetaroot-alpha0root
                    alpha_geo[i] = alpha+thetaroot
                    a[i] = a0root
                else:
                    alp[i] = alpha+thetaroot-alpha0root+(thetatip-alpha0tip-thetaroot+alpha0root)*(y[i]-trans)/(0.5*span-trans)
                    alpha_geo[i] = alpha+thetaroot+(thetatip-thetaroot)*(y[i]-trans)/(0.5*span-trans)                
                    a[i] = a0root+(a0tip - a0root)*(y[i]-trans)/(0.5*span-trans)
                    
        else:
            Re_asa = self._Re(self._cordas)
            
            for i in range(n):
                dCl = self._perfil_raiz.ajustelinear(Re_asa[i],self._mach).get('dCl')
                a0 = dCl[0]*180/np.pi
                alpha0 = np.radians(-dCl[1]/dCl[0])
                
                alp[i] = alpha+thetaroot-alpha0
                a[i] = a0
                
                if y[i] < trans:
                    alpha_geo[i] = alpha+thetaroot
                else:
                    alpha_geo[i] = alpha+thetaroot+(thetatip-thetaroot)*(y[i]-trans)/(0.5*span-trans)
        
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
            
    def de_da(self,x,z):
        """

        * Cálculo da variacao do coeficiente de "downwash" do plano de simetria
          da asa para a posicao do CA do EH.
        * 'x' e 'z' representam as coordenadas do centro aerodinâmico do EH
          relativas ao mesmo da corda da asa,  com 'x' no sentido do vento
          e 'z' orientado para cima, com z = 0 na altura da asa.
        * As contas aqui aplicadas foram baseadas no NASA Technical Paper 2414,
          'Downwash in the plane of simmetry of an elliptically loaded wing'
        """
        def e(alfa):
            s = 0.5*self.bw
            y = self._y
            gamma = self._mod_llt(alfa,_dw=True) # circulacao normalizada pela velocidade (Gamma/V)
            gamma = UnivariateSpline(y,gamma,k=3) # Spline dos pontos
            dg_dy = gamma.derivative() #Derivada da circulacao 
            
            def de_dy(y):
                """
                Integracao numeria baseada na formula 7 do artigo
                """
                k1 = x/np.sqrt(x**2+y**2+z**2)  
                k2 = y**2/(x**2+z**2)
                k3 = y**2/(y**2+z**2)
                
                de_dy = -(k1*(k2+k3)+k3)*dg_dy(y)/(2*np.pi*y)
                return de_dy
            
            e = np.array(quad(de_dy,0,s))
            return e*57.29577951308232 #em graus

        print('\nCalculando downwash da asa para a posicao (x,z) = (%.3f,%.3f)...' %(x,z))

        e10, err10 = e(10.0)
        e0, err0 = e(0.0)
        de_da = 0.1*(e10-e0)
        err = 0.1*np.sqrt(err10**2+err0**2)
        
        print('         de/dalfa = %.3f' %(de_da))  
        print('         Erro estimado: %e' %(err))
            
        return de_da
            
    def _estol(self,pres=1):
        """
        Função que faz uma estimativa do angulo de estol e CLmax da asa
        em função dos cl locais comparados ao clmax do perfil. O angulo é
        encontrado através de um algoritmo de procura por bisseção. Esse é 
        um processo bastante lento se self._n for muito grande.
        """
        print("\nCalculando angulo de estol da asa...")
        
        f = lambda x: self._mod_llt(x,_dist=True)
        
        contador = 0
        ini = 5
        fim = 25
        mult = 10**pres
        fim *= mult
        ini *= mult
        mult = float(mult)
        
        while ini < fim:
            if contador > 20:
                raise ValueError("Nao houve convergencia para o  estol da asa")
            else:
                contador += 1
                pos = (ini+fim)/2
                array = f(pos/mult)[3]-self._clmax_local()
                
                # Verifica se todos os elementos são negativos
                if  np.all(np.less(array,0.0)):
                    ini = pos + 1
                else:
                    fim = pos - 1
            
        print("         Numero de iteracoes: %i" %(contador))
        print("         Calculo completo.")
        alpha_estol = pos/mult
        CLmax = self._mod_llt(alpha_estol,_iter=True)[0]
        return alpha_estol,CLmax
            
    def _polares(self,lista_alfa = range(10,3)):
        """
        Função responsavel por gerar as polares da asa, tanto de arrasto
        quanto de sustentacao. Para isso são calculados os coeficientes 
        aerodinamicos para alguns angulos. Depois esses pontos sao interpolados.
        """
        lista_alfa = np.array(lista_alfa)
        n_alfa = len(lista_alfa)
        
        CL = np.zeros(n_alfa)
        CD0 = np.zeros(n_alfa)
        CDi = np.zeros(n_alfa)

        for i in range(n_alfa): 
            val = self._mod_llt(lista_alfa[i], _iter = True)
            CL[i] = val[0]
            CD0[i] = val[1]
            CDi[i] = val[2] 
            
        def sinal(b):
                if b >= 0.0:
                    sinal = '+'
                else:
                    sinal = '-'
                return sinal
        
        def sustentacao():
            """
            Plota a polar de sustentacao e retorna valores de dCl/dAlpha
            e alphaCl0 da asa. 
            """
            a,b = np.polyfit(lista_alfa,CL,deg=1)
            alphaCL0 = -b/a 
            funcao_CL = np.poly1d([a,b])
            
            def grafico():
                alphas = np.append(alphaCL0,lista_alfa)
                CLs = np.append(0.0,CL)
                legenda = r'$C_L = %.3f*(\alpha %s %.3f)$' %(a,sinal(-alphaCL0),abs(-alphaCL0))
                
                plt.figure()
                plt.grid('on')
                plt.title(r'Polar de sustentacao da asa')
                plt.xlabel(r'$\alpha$ ($graus$)')
                plt.ylabel(r'$C_L$')
                plt.plot(alphas,CLs,'b-',label = legenda)
                xmin,xmax = plt.xlim()
                plt.xlim(min(alphas),xmax+(xmax-xmin)*0.05)
                plt.legend(loc='best')
                plt.savefig("polarLiftAsa.png", bbox_inches='tight', dpi=200)
                plt.close()
                
            if self._p:
                apoio.executarNaPasta(self._pasta_graficos)(grafico)()
            
            return funcao_CL, alphaCL0
        
        def arrasto():    
            """
            Plota a polar de arrasto e retorna o polinomio de 2º grau
            interpolado dos pontos calculados. 
            """        
            pol_CD = np.polyfit(CL,CD0+CDi,deg=2)
            funcao_CD = np.poly1d(pol_CD)
            
            pol_CDi = np.polyfit(CL,CDi,deg=2)
            funcao_CDi = np.poly1d(pol_CDi)
            
            def grafico():
                CLs = np.linspace(0.0,max(CL))                
                
                if abs(pol_CD[1]) <= 1e-3: 
                    legenda = r'$C_D = %.3f %s %.3f*C_L^2 $' %(pol_CD[2], 
                                                               sinal(pol_CD[0]), abs(pol_CD[0]))
                else:
                    legenda = r'$C_D = %.3f %s %.3f*C_L %s %.3f*C_L^2$' %(pol_CD[2], 
                                                                          sinal(pol_CD[1]), abs(pol_CD[1]),
                                                                          sinal(pol_CD[0]), abs(pol_CD[0]))                
                    
                plt.figure()
                plt.grid('on')
                plt.title(r'Polar de arrasto da asa')
                plt.xlabel(r'$C_L$')
                plt.ylabel(r'$C_D$')
                plt.plot(CLs,funcao_CD(CLs),'b-',label = legenda)
                plt.legend(loc='best')
                plt.savefig("polarDragAsa.png", bbox_inches='tight', dpi=200)
                plt.close()
                
            if self._p:
                apoio.executarNaPasta(self._pasta_graficos)(grafico)()
            
            return funcao_CD, funcao_CDi
        
        funcao_CL, alphaCL0 = sustentacao()
        funcao_CD,funcao_CDi = arrasto()
        
        return alphaCL0, funcao_CL, funcao_CD, funcao_CDi
        
    def _CM(self):
        """
        Funçao que estima o coeficiente de momento da asa por integração dos
        CM's dos perfis.
        """
        dy = self._dy
        area_secao = (self._cordas[1:]+self._cordas[:-1])*dy/2.0
        n = len(area_secao)
        alpha_secao = 0.5*(self._alpha_local[:-1]+self._alpha_local[1:])
        y_secao = 0.5*(self._y[:-1]+self._y[1:])
        trans = self.trans
        mach = self._mach
        ca = self._p_c
        Cm_secao = np.zeros(n)  
        corda_secao = 0.5*(self._cordas[:-1]+self._cordas[1:])
        
        if self._perfil_raiz.name != self._perfil_ponta.name:
            Re_raiz = self._Re(self._cordas[-1])
            Re_ponta = self._Re(self._cordas[0])
            
            for i in range(n):
                ang = alpha_secao[i]
                if y_secao[i] < trans:
                    Cl = self._perfil_raiz.interpReMachAlfa(ang,Re_raiz,mach,'Cl')
                    Cd = self._perfil_raiz.interpReMachAlfa(ang,Re_raiz,mach,'Cd')
                    xCp = self._perfil_raiz.interpXcp(ang,Re_raiz,mach)
                     
                else:
                    Cl_raiz = self._perfil_raiz.interpReMachAlfa(ang,Re_raiz,mach,'Cl')
                    Cl_ponta = self._perfil_ponta.interpReMachAlfa(ang,Re_ponta,mach,'Cl')
                    Cd_raiz = self._perfil_raiz.interpReMachAlfa(ang,Re_raiz,mach,'Cd')
                    Cd_ponta = self._perfil_ponta.interpReMachAlfa(ang,Re_ponta,mach,'Cd')
                    xCp_raiz = self._perfil_raiz.interpXcp(ang,Re_raiz,mach)
                    xCp_ponta = self._perfil_ponta.interpXcp(ang,Re_ponta,mach)
                    
                    Cl = Cl_raiz + (Cl_ponta-Cl_raiz)*(y_secao[i]-trans)/(0.5*self.bw-trans)
                    Cd = Cd_raiz + (Cd_ponta-Cd_raiz)*(y_secao[i]-trans)/(0.5*self.bw-trans)
                    xCp = xCp_raiz + (xCp_ponta-xCp_raiz)*(y_secao[i]-trans)/(0.5*self.bw-trans)
                                 
                Cm_secao[i] = (ca[i]-xCp)*(Cl*np.cos(ang)-Cd*np.sin(ang))
                
        else:
            Re_asa = self._Re(self._cordas)
            
            for i in range(n):
                ang = alpha_secao[i]
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
        
    def _distribuicoes(self,lista_alfa):
        """
        Função que plota as distribuições de diversas grandezas ao longo da asa
        """
        lista_alfa = np.array(lista_alfa)
        n_alfa = len(lista_alfa)
        clmax_local = self._clmax_local()        

        gamma = [0.0 for i in range(n_alfa)]
        alpha_i = [0.0 for i in range(n_alfa)]
        alpha_local = [0.0 for i in range(n_alfa)]
        cl_local = [0.0 for i in range(n_alfa)]
        cd0_local = [0.0 for i in range(n_alfa)]
        cdi_local = [0.0 for i in range(n_alfa)]
        CL = [0.0 for i in range(n_alfa)]
        
        for i in range(n_alfa):        
            val = self._mod_llt(lista_alfa[i],_dist = True)
            gamma[i] = val[0]
            alpha_i[i] = val[1]
            alpha_local[i] = val[2]            
            cl_local[i] = val[3]
            cd0_local[i] = val[4]
            cdi_local[i] = val[5]
            CL[i] = val[6]
            
        def geometrico():
            plt.figure()
            plt.grid('on')
            plt.title(r'Angulo geometrico local')
            plt.xlabel(r'Semi-envergadura normalizada ($2y/b$)')
            plt.ylabel(r'$\alpha_{geo}$ ($graus$)')
            for i in range(n_alfa):
                ageo = (alpha_local[i]+alpha_i[i])*180/np.pi
                plt.plot(2.0*self._y/self.bw,ageo,label=r'$\alpha = %.1f^o$' %(lista_alfa[i]))
            ymin,ymax = plt.ylim()
            plt.ylim(ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.05)
            plt.legend(loc='best')
            plt.savefig("distAgeo.png", bbox_inches='tight', dpi=200)
            plt.close()
        
        def local():
            plt.figure()
            plt.grid('on')
            plt.title(r'Angulo local')
            plt.xlabel(r'Semi-envergadura normalizada ($2y/b$)')
            plt.ylabel(r'$\alpha$ ($graus$)')
            for i in range(n_alfa):
                ai = alpha_local[i]*180/np.pi
                plt.plot(2.0*self._y/self.bw,ai,label=r'$\alpha = %.1f^o$' %(lista_alfa[i]))
            plt.legend(loc='best')
            plt.savefig("distAlocal.png", bbox_inches='tight', dpi=200)
            plt.close()
        
        def induzido():
            plt.figure()
            plt.grid('on')
            plt.title(r'Angulo induzido local')
            plt.xlabel(r'Semi-envergadura normalizada ($2y/b$)')
            plt.ylabel(r'$\alpha_i$ ($graus$)')
            for i in range(n_alfa):
                ai = -alpha_i[i]*180/np.pi
                plt.plot(2.0*self._y/self.bw,ai,label=r'$\alpha = %.1f^o$ ' %(lista_alfa[i]))
            plt.legend(loc='best')
            plt.savefig("distAi.png", bbox_inches='tight', dpi=200)
            plt.close()
        
        def cl():
            plt.figure()
            plt.grid('on')
            plt.title(r'Coeficiente de sustentacao local')
            plt.xlabel(r'Semi-envergadura normalizada ($2y/b$)')
            plt.ylabel(r'$C_l$')
            for i in range(n_alfa):
                cl = cl_local[i]
                plt.plot(2.0*self._y/self.bw,cl,label=r'$\alpha = %.1f^o$' %(lista_alfa[i]))
            plt.plot(2.0*self._y/self.bw,clmax_local, 'k--',label='$C_{l_{max}}$ $2D$')
            ymin,ymax = plt.ylim()
            plt.ylim(ymax=ymax+(ymax-ymin)*0.05)
            plt.legend(loc='best')
            plt.savefig("distCl.png", bbox_inches='tight', dpi=200)
            plt.close()
        
        def lift():
            plt.figure()
            plt.grid('on')
            plt.title(r'Distribuicao de sustentacao')
            plt.xlabel(r'Semi-envergadura ($m$)')
            plt.ylabel(r'$L(N/m)$')
            for i in range(n_alfa):
                elipt = np.sin(self._theta)*2*self._ar.rho()*self.vel**2*self.Sw*CL[i]/(np.pi*self.bw)
                L = 0.5*self._ar.rho()*self.vel**2*self._cordas*cl_local[i]
                plt.plot(self._y,elipt,'k--')
                plt.plot(self._y,L,label=r'$\alpha = %.1f^o$' %(lista_alfa[i]))
            plt.legend(loc='best')
            plt.savefig("distLift.png", bbox_inches='tight', dpi=200)
            plt.close()

        def cd0():
            plt.figure()
            plt.grid('on')
            plt.title(r'Coeficiente de arrasto parasita local')
            plt.xlabel(r'Semi-envergadura normalizada ($2y/b$)')
            plt.ylabel(r'$C_{D_0}$')
            for i in range(n_alfa):
                cd0 = cd0_local[i]
                plt.plot(2.0*self._y/self.bw,cd0,label=r'$\alpha = %.1f^o$' %(lista_alfa[i]))
            plt.legend(loc='best')
            plt.savefig("distCd0.png", bbox_inches='tight', dpi=200)
            plt.close()

        def cdi():
            plt.figure()
            plt.grid('on')
            plt.title(r'Coeficiente de arrasto induzido local')
            plt.xlabel(r'Semi-envergadura normalizada ($2y/b$)')
            plt.ylabel(r'$C_{D_i}$')
            for i in range(n_alfa):
                cdi = cdi_local[i]
                plt.plot(2.0*self._y/self.bw,cdi,label=r'$\alpha = %.1f^o$' %(lista_alfa[i]))
            plt.legend(loc='best')
            plt.savefig("distCdi.png", bbox_inches='tight', dpi=200)
            plt.close()

        def drag():
            plt.figure()
            plt.grid('on')
            plt.title(r'Distribuicao de Arrasto')
            plt.xlabel(r'Semi-envergadura ($m$)')
            plt.ylabel(r'$D(N/m)$')
            for i in range(n_alfa):
                D = 0.5*self._ar.rho()*self.vel**2*self._cordas*(cd0_local[i]+cdi_local[i])
                plt.plot(self._y,D,label=r'$\alpha = %.1f^o$' %(lista_alfa[i]))
            plt.legend(loc='best')
            plt.savefig("distDrag.png", bbox_inches='tight', dpi=200)
            plt.close()
        
        def plotar():
            """
            Rodando todos os plots por aqui, o script só tem que 
            trocar de pasta uma única vez.
            """
            #distribuicões de angulos
            geometrico(); induzido(); local()
            #distribuições de coeficientes
            cl(); cd0(); cdi()
            #distribuições de forças
            lift(); drag()
            
        apoio.executarNaPasta(self._pasta_graficos)(plotar)()
        
    def _superior(self):
        """
        Gráfico da asa em planta comparada à asa elíptica
        """
        # Gera o gráfico da asa em partes
        elipt_pos = []
        elipt_neg = [] 
        a = self._elipt(llt=True)
        planta_pos = []
        planta_neg = []
        b = self._cordas 
        y = self._y
        p_c = self._p_c
                 
        for i in range(len(y)):
            elipt_pos.append(a[i]*p_c[i])
            elipt_neg.append(-a[i]*(1-p_c[i]))               
            planta_pos.append(b[i]*p_c[i])
            planta_neg.append(-b[i]*(1-p_c[i]))
        
        plt.figure()
        plt.grid('on')
        plt.axis('equal')
        plt.title('Geometria da semi-asa')
        plt.xlabel('Semi-envergadura ($m$)')
        plt.ylabel('Corda ($m$)')
        plt.plot(y, elipt_pos, 'r-') 
        plt.plot(y, elipt_neg, 'r-', label='Asa eliptica')
        plt.plot(y, planta_pos, 'b-')
        plt.plot(y, planta_neg, 'b-', label='Asa em planta')
        plt.legend(loc='best')
        plt.fill_between(y, elipt_pos, elipt_neg, color='r', alpha=0.3)
        plt.fill_between(y, planta_pos, planta_neg, color='b', alpha=0.3)             
        plt.savefig("geometriaAsa.png", bbox_inches='tight', dpi=200)
        plt.close()
        
if __name__ == "__main__":
#    Geometria da aeronave:
#         Posicao CG :           [-0.03036109  0.          0.19970691] m
#         Xcg :                  0.401748
#        ---------- ASA FRONTAL ---------
#         X do bordo de ataque : -0.190730 m
#         Z do bordo de ataque : 0.120000 m
#         Area :                 0.797030 m^2
#         Envergadura :          1.996680 m
#         Corda :                0.399177 m
#         Angulo de incidencia : 2.603011 graus
#         Angulo de torsao :     -0.588582 graus
#        --------- ASA TRASEIRA ---------
#         X do bordo de ataque : 0.425401 m
#         Z do bordo de ataque : 0.185548 m
#         Area :                 0.326136 m^2
#         Envergadura :          1.130175 m
#         Corda :                0.288571 m
#         Angulo de incidencia : 1.175684 graus
#         Angulo de torsao :     2.865253 graus
#        -------------- EV --------------
#         Area :                 0.034810 m^2
#         Envergadura :          0.122112 m
#         Corda :                0.281557 m
#        ------------ MOTOR -------------
#         Posicao :              [-0.61690474  0.          0.2       ] m    

#    def f(x):
#        Asa = S_sustentadora("S1223 MOD2015", "S1223 MOD2015", 2.6, 0.797, 1.997, x, 0.0, 13.7, 25)
#        return Asa.CDi
#        
#    print fmin(f,0.5,maxfun=30)
    
    Asa1 = S_sustentadora("S1223 MOD2015", "S1223 MOD2015", 19.4, 0.797, 1.997, 0.65,
                          -3.0, 20.0, 25, p = False)
    
    print Asa1.raiz,Asa1.ponta,Asa1.offset,Asa1.alfa_estol,Asa1.trans,Asa1.mac
