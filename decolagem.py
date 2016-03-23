# -*- coding: utf-8 -*-
"""
Created on Fri Oct 03 23:54:50 2014

@author: Diego Chou, Gustavo de Almeida, Mateus Pereira
"""

print("\nCarregando modulos de 'Decolagem'...")

try:
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import odeint, quad
    print("Modulos de 'Decolagem' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Decolagem'\n")
    raise

import apoio
import geral

"""    
 Lista de variáveis e funções (unidades no SI, salvo indicado):

 Comprimento de Pista                                  c           INPUT
 Aceleraçao da Gravidade                               g           INPUT
 Altitude densidade                                    alt         INPUT
 Área de referência (da asa)                           Sw          INPUT
 Cl da asa em ângulo de corrida                        Cl          INPUT
 Cl máximo da asa                                      Clm         INPUT
 Cd de corrida                                         Cd          INPUT
 Massa do avião SEM carga                              m_a         INPUT
 Coeficiente de atrito com a pista                     mi          INPUT
"""

class Analise(object):
    """
    Análise de decolagem para definiçao da carga paga máxima    
    """
    #Parametros de pista
    mi = 0.04
    c = 70.0
    vel_inicial = 0.0
    
    # Parametros de motor
    motor = 'OS61FX'
    helice = 'APC12.25x3.75'
    alt_teste = 408
    trac_tabela = [40.6949,40.00848,39.32206,36.2822,30.98696,29.61412,27.75098,18.6314]
    v_tabela = [0.0,3.0,6.5,10.0,12.3,15.0,17.9,20.4]
    
    def __init__(self, Sw, m_a, Cl, Clmax, Cd, p = False):
        #Leitura do arquivo de configuração
        self.p = p
        
        self.g = geral.gravidade 
        self.alt = geral.Ar().alt
        self.Sw = Sw
        self.Cl = Cl
        self.Clmax = Clmax
        self.Cd = Cd
        self.m_a = m_a
        self.rho = geral.Ar().rho()
        self.rho_mar = geral.Ar().rho_mar
        
        # Jogar o polinomio para o __init__ da classe garante que ele seja 
        # calculado apenas uma vez, e não a cada vez que a função é chamada.
        self.polinomio_tracao = self.poly_t()
        
    def poly_t(self, n=3):
        """
        Traçao disponível. A função lê uma lista de valores de tração
        e outra da velocidade correspondente e interpola a tração disponível
        por um polinômio de grau n.
        """
        v_tabela = self.v_tabela
        trac_tabela = self.trac_tabela
        fit = np.polyfit(v_tabela,trac_tabela,3)
        fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
        
        @apoio.executarNaPasta('Graficos/Desempenho')
        def grafico():            
            x = np.linspace(v_tabela[0],v_tabela[-1],100)
            
#            T = lambda x : (-(1e-4)*x**4+(2.6e-3)*x**3-(2.65e-2)*x**2-\
#                           (2.754e-1)*x+3.87*self.g)            
            
            plt.figure()
            plt.grid('on')
            plt.title('Curva de tracao %s - %s' %(self.motor,self.helice))
            plt.plot(x, fit_fn(x), label='2015')
#            plt.plot(x, T(x), label='2013') 
            plt.errorbar(v_tabela,trac_tabela, yerr=2.0, fmt='o')
            plt.xlabel(r'$V$ ($m/s$)')
            plt.ylabel(r'$T$ ($N$)')
            plt.legend(loc='best')
            plt.savefig("curva_%s_%s.png" %(self.motor, self.helice), bbox_inches='tight', dpi=200)

        if self.p:
            grafico()
        
        return fit_fn
        
    def t_disp(self, vel): 
        
        rho = self.rho
        
        rho_teste = geral.Ar_qualquer(self.alt_teste).rho()
        
        T = self.polinomio_tracao(vel)*rho/rho_teste       
        
        return T
    
    def lift(self, vel):
        """
        Cálculo da força de sustentaçao
        """
        rho = self.rho
        Sw = self.Sw
        Cl = self.Cl
        L = rho*Sw*vel**2*Cl/2.0     
        return L
        
    def drag(self, vel):
        """
        Cálculo da força de arrato
        """
        rho = self.rho
        Sw = self.Sw
        Cd = self.Cd
        D = rho*Sw*vel**2*Cd/2.0
        return D
                
    def massa(self, CPaga):
        """
        Massa total do aviao
        """
        m = self.m_a
        M_t = m + CPaga
        return M_t
    
    def fat(self, vel, CPaga):  
        """
        Atrito com o chão.
        É necessário rever o coeficiente de atrito com o chão...
        """
        mi = self.mi
        m = self.massa(CPaga)
        g = self.g
        L = self.lift(vel)
        Fat = mi*(m*g-L)
#        print "Fat: ", Fat
#        print "Vel, CPaga :", vel, CPaga        
        return Fat
            
    def accel(self, vel, CPaga):
        """
        Cálculo da aceleração
        """
        T = self.t_disp(vel)
        D = self.drag(vel)
        Fat = self.fat(vel,CPaga)
        m = self.massa(CPaga)
        a1 = (T-D-Fat)/m
        return a1
                        
    def vel_estol(self, CPaga):
        """
        Cálculo da velocidade de estol
        """
        m = self.massa(CPaga)
        g = self.g
        rho = self.rho
        Sw = self.Sw
        Clmax = self.Clmax
        v = np.sqrt(2.0*m*g/(rho*Sw*Clmax))
        return v
                            
    def vel_decol(self, CPaga):
        """
        Velocidade de decolagem. O ângulo de decolagem é o de Clmax. 
        O calculo é feito com uma margem de 20%
        """
        v_estol = self.vel_estol(CPaga)
        v = 1.2*v_estol
        return v

    def vel_corrida(self, CPaga):
        """
        Velocidade de corrida calculada ponto a ponto pelo método diferencial
        """
        # self.c - 0.4*self.vel_estol(CPaga) é a distância de pista menos a distância de rotação (arfagem).
        distance = np.linspace(0, self.c - 0.4*self.vel_estol(CPaga), 101)        
        inicial = self.vel_inicial+0.001 # O valor inicial da velocidade nao pode ser zero.      
        
        def deriv(y, t):
            f = self.accel(y, CPaga)/y
            return f
            
        y = odeint(deriv, inicial, distance) # Solução da EDO dada pelo scipy

        @apoio.executarNaPasta('Graficos/Desempenho')
        def grafico():
            plt.figure()
            plt.grid('on')
            plt.plot(distance, y)
            plt.plot(distance[-1], y[-1], 'bo')
            plt.title('Velocidade de corrida na decolagem')
            plt.xlabel(r'$D$ ($m$)')
            plt.ylabel(r'$V$ ($m/s$)')
            plt.annotate('Ponto de rotacao:\n $V = %.1f$ $m/s$ \n $D = %.1f$ $m$' %(y[-1], distance[-1]),
                         xy = (distance[-1], y[-1]), xytext = (-7, -100),
                         textcoords = 'offset points', ha = 'right', va = 'bottom',
                         bbox = dict(boxstyle = 'round,pad=0.5', fc = 'white'),
                         arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            plt.savefig("velocidadeCorrida.png", bbox_inches='tight', dpi=200) 
            
        grafico()

#        return y[len(distance)-1]        

    # Metodo integral  
    def s_corrida(self, CPaga):
        """
        Distância percorrida na decolagem calculada pelo método integral
        proposto por Raymer
        """
        
        def f(vel):
            mi = self.mi
            g = self.g
            m = self.massa(CPaga)
            T = self.t_disp(vel)
            D = self.drag(vel)
            L = self.lift(vel)
            f = (m*vel)/(T-D-mi*(m*g-L))
            return f
            
        v_estol = self.vel_estol(CPaga)
        s = quad(f, self.vel_inicial, 1.2*v_estol)[0] + 0.4*v_estol
        return s
        
    def c_paga_max(self, pres=3):
        """
        Caso em que vel_decol == vel_corrida no limite da pista
        Algoritmo de Busca Binária (Bisseção) para encontrar de maneira rápida a c_paga_max
        A variável pres representa o número de casas decimais de precisão do cálculo da c_paga_max        
        """
        print("\nCalculando carga paga maxima...")
        
        contador = 0
        ini = 0
        fim = 30
        mult = 10**pres
        fim *= mult
        mult = float(mult)
        
        while ini <= fim:
            contador += 1
            pos = (ini+fim)/2

            if self.s_corrida(pos/mult) > self.c:
                fim = pos - 1
                j = 0

            else:
                ini = pos + 1
                j = 1
                
        print("         Numero de iteracoes: %i" %(contador))
        if j == 1 :
            print("         Calculo completo.")
            return ini/mult

        else :
            print("         Calculo completo.")
            return fim/mult
        
