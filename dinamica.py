# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 12:10:49 2016

@author: Diego Chou
"""

print("\nCarregando modulos de 'Estabilidade Dinamica'...")
try:
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import lti, step, impulse
    print("Modulos de 'Estabilidade Dinamica' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Estabilidade Dinamica'\n")
    raise

from apoio import executarNaPasta

def estabilidade_dinamica(name,p):
    """
    Analise da estabilidade dinamica
    """
    @executarNaPasta('Runs')
    def getValues():
            """
            Funcao que lê os parâmetros de estabilidade dinamica
            de o arquivo especificado.
            """
            eig = []
            filename = '%s_values.eig' %(name)
            f = open(filename, 'r')
            flines = f.readlines()
            for i in range(3,len(flines)):
                words = flines[i].split()
                num = complex(float(words[1]),float(words[2]))
                eig.append(num)
            f.close()
            
            return eig
            
    def define_modes():
        """
        Funçao que recebe os autovalores e discrimina os modos de vibração. 
        """
        # Separacao entre autovalores oscilantes e nao oscilantes
        non_oscilate = []
        oscilate = []  
        for i in eig:
            if np.imag(i) == 0.0:
                non_oscilate.append(float(np.real(i)))
            else:
                oscilate.append(i)
        
        # Constante de amortecimento e frequencia natural
        csi = -np.real(oscilate)/np.absolute(oscilate)
        omega_amort = abs(np.imag(oscilate))
        omega = np.absolute(oscilate)
        oscilate = list(set(zip(csi,omega,omega_amort)))
        oscilate = np.array(sorted(oscilate,key=lambda x: x[0]))
        
        # Tempo de meia vida
        non_oscilate = np.sort(non_oscilate)
        
        if len(oscilate) != 3:
            raise ValueError("Existem apenas %d de 3 modos oscilantes" %(len(oscilate)))
        
        long_val = [oscilate[0],oscilate[2]]
        latdir_val = [non_oscilate[0],non_oscilate[1],oscilate[1]]
        latdir_t = -np.log(2)/non_oscilate
        
        # Calculo dos parametros latero-direcionais
        
        
        return long_val,latdir_val,latdir_t
        
    def longitudinal():
        """
        Funcao que analisa a estabilidade dinamica longitudinal
        """
        phug_csi, phug_omega, phug_omega_amort = long_val[0]
        short_csi,short_omega, short_omega_amort = long_val[1]
    
        def root_locus():
            """
            Gera o grafico de root locus para phugoid e short-period.
            """            
            X1 = -phug_csi*phug_omega
            Y1 = phug_omega*(1-phug_csi**2)**0.5
            X2 = X1
            Y2 = -Y1
            X3 = -short_csi*short_omega
            Y3 = short_omega*(1-short_csi**2)**0.5
            X4 = X3
            Y4 = -Y3
            
            def grafico():
                plt.figure()
                plt.grid('on')
    #            plt.axis('equal')
                plt.plot(X1, Y1, 'bo', label='Polos Phugoid')
                plt.plot(X2, Y2, 'bo')
                plt.plot(X3, Y3, 'ys', label='Polos Short-Period')
                plt.plot(X4, Y4, 'ys')
                plt.title('Root Locus para Phugoid e Short-Period')
                plt.xlabel('Eixo Real')
                plt.ylabel('Eixo Imaginario')
                plt.legend(loc='best')
                plt.savefig("root_long.png", bbox_inches='tight', dpi=200)
            
            grafico()
            
        def phugoid():
            """
            Gera o grafico da variacao de altitude com o tempo.
            """            
            # num and den, can be list or numpy array type        
            num = [1.0] 
            den = [1.0, 2*phug_csi*phug_omega, phug_omega**2]
             
            tf = lti(num, den)
             
            # get t = time, s = unit-step response
            t, s = step(tf)
                 
            # recalculate t and s to get smooth plot
            t, s = step(tf, T = np.linspace(min(t), t[-1], 500))
             
            # get i = impulse
            t, i = impulse(tf, T = np.linspace(min(t), t[-1], 500))        
            
            def grafico():
                plt.figure()
                plt.grid('on')
    #            plt.axis('equal')
                plt.plot(t, i, 'b-')
                plt.title('Phugoid')
                plt.xlabel('Tempo (s)')
                plt.ylabel('Variacao da altitude (feet)')
                plt.legend(loc='best')
                plt.savefig("phugoid.png", bbox_inches='tight', dpi=200)
                
            grafico()                      
            
        def short_period():
            """
            Gera o grafico da variacao do angulo de ataque com o tempo.
            """            
            # num and den, can be list or numpy array type       
            num = [1.0] 
            den = [1.0, 2*short_csi*short_omega, short_omega**2]
             
            tf = lti(num, den)
             
            # get t = time, s = unit-step response
            t, s = step(tf)
             
            # recalculate t and s to get smooth plot
            t, s = step(tf, T = np.linspace(min(t), t[-1], 500))
             
            # get i = impulse
            t, i = impulse(tf, T = np.linspace(min(t), t[-1], 500))        
    
            def grafico():
                plt.figure()
                plt.grid('on')
    #            plt.axis('equal')
                plt.plot(t, i, 'b-')
                plt.title('Short-Period')
                plt.xlabel('Tempo (s)')
                plt.ylabel('Variacao do angulo de ataque (rad)')
                plt.legend(loc='best')
                plt.savefig("short_period.png", bbox_inches='tight', dpi=200)
                
            grafico()   
         
        @executarNaPasta('Graficos/Estabilidade') 
        def plot():
            root_locus()
            phugoid()
            short_period()
        
        def pontuacao():
            '''
            Analise principal.
            Valores para o amortecimento do modo Phugoid maiores que 0 sao considerados bons
            Valores para o amortecimento do modo Short-Period entre 0.2 e 2 sao considerados bons
            O aviao foi classificado em Class III, Category B e Level 2
            '''
            
            print('\nResultados da Estabilidade Dinamica Longitudinal:')
            print('        ----------------Phugoid----------------')
            print('                Razao de amortecimento : %.5f' %(phug_csi))
            print('                Frequencia de oscilacao : %.5f rad/s' %(phug_omega_amort))
            print('        --------------Short-Period-------------')
            print('                Razao de amortecimento : %.5f' %(short_csi))
            print('                Frequencia de oscilacao : %.5f rad/s' %(short_omega_amort))
                
            fator = 1.0
            
            if (phug_csi>=0) and (short_csi>=0.2) and (short_csi<=2.0):
                print('         Os valores sao ACEITAVEIS')
            else:
                print('         Os valores NAO sao aceitaveis')
                if phug_csi < 0:
                    fator *= np.exp(phug_csi)
                if short_csi < 0.2:
                    fator *= np.exp(short_csi-0.2)
                if short_csi > 2:
                    fator *= np.exp(2.0-short_csi)
            
            print('         Fator de pontuacao : %f' %(fator))
            return fator
        
        if p: plot()
        return pontuacao()
            
    def latdir():
        """
        Funcao que analisa a estabilidade dinamica longitudinal
        """
        roll,spiral = latdir_val[:2]
        dutch_csi, dutch_omega, dutch_omega_amort = latdir_val[-1]
        tau,t_meio = latdir_t
    
        def root_locus():
            """
            Gera o grafico de root locus para Spiral, Roll e Dutch Roll.
            """            
            X1 = roll
            Y1 = 0
            X2 = spiral
            Y2 = 0
            X3 = -dutch_omega * dutch_csi
            Y3 = dutch_omega * (1 - dutch_csi**2)**0.5
            X4 = X3
            Y4 = -Y3
    
            @executarNaPasta('Graficos/Estabilidade')        
            def grafico():
                plt.figure()
                plt.grid('on')
    #            plt.axis('equal')
                plt.plot(X1, Y1, 'bo', label='Roll')
                plt.plot(X2, Y2, 'yv', label='Spiral')
                plt.plot(X3, Y3, 'ms', label='Dutch Roll')
                plt.plot(X4, Y4, 'ms')
                plt.title('Root Locus para Roll, Spiral e Dutch-Roll')
                plt.xlabel('Eixo Real')
                plt.ylabel('Eixo Imaginario')
                plt.legend(loc='best')
                plt.savefig("root_latdir.png", bbox_inches='tight', dpi=200)
            
            grafico()
            
        def pontuacao():
            '''
            Analise principal.
            O aviao foi classificado em Class III, Category B e Level 2            
            Valores para o amortecimento do modo Dutch Roll maiores que 0.08 sao considerados bons
            Valores para a frequencia do modo Dutch Roll maiores que 0.4 sao considerados bons
            Valores para amort*freq do modo Dutch Roll maiores que 0.15 sao considerados bons
            Valores para tau do modo Roll menores que 3 sao considerados bons
            Valores para t_meio do modo Spiral maiores que 12 sao considerados bons
            '''    
            print('\nResultados da Estabilidade Dinamica Latero-Direcional:')
            print('        ----------------Dutch-Roll----------------')
            print('                Razao de amortecimento : %.5f' %(dutch_csi))
            print('                Frequencia de oscilacao : %.5f rad/s' %(dutch_omega_amort))
            print('                Amortecimento vezes frequencia : %.5f rad/s' %(dutch_csi*dutch_omega_amort))
            print('        -------------------Roll-------------------')
            print('                Constante de tempo : %.5f s' %(tau))
            print('        ------------------Spiral------------------')
            print('                Tempo para a amplitude cair a metade: %.5f s' %(t_meio))
            
            fator = 1.0
            
            if dutch_csi>=0.08 and dutch_omega_amort>=0.4 and dutch_csi*dutch_omega_amort>=0.15 and t_meio>=12.0 and tau<=3.0:
                print('         Os valores sao ACEITAVEIS')
            else:
                print('         Os valores NAO sao aceitaveis')
                    
                if dutch_csi < 0.08:
                    fator *= np.exp(dutch_csi-0.08)
                    
                if dutch_omega_amort < 0.4:
                    fator *= np.exp(dutch_omega_amort-0.4)
                    
                if dutch_csi*dutch_omega_amort < 0.15:
                    fator *= np.exp(dutch_csi*dutch_omega_amort-0.15)
                    
                #if t_meio < 12.0:
                #    fator *= np.exp((t_meio-12.0)/12.0)
                    
                if tau > 3.0:
                    fator *= np.exp(3-tau)
            
            print('         Fator de pontuacao : %f' %(fator))
            return fator
            
        if p: root_locus()
        return pontuacao()
    
    eig = getValues()
    long_val,latdir_val,latdir_t = define_modes()
    fator_long = longitudinal()
    fator_latdir = latdir()
    fator = fator_long*fator_latdir
    
    return fator
