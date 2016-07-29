# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 12:21:22 2016

@author: Diego Chou
"""

print("\nCarregando modulos de 'Estabilidade Estatica'...")
try:
    import numpy as np
    import matplotlib.pyplot as plt
    print("Modulos de 'Estabilidade Estatica' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Estabilidade Estatica'\n")
    raise
    
from apoio import executarNaPasta
    
def estabilidade_estatica(name,trim_desejado,alfa_estol,vel,
                          pos_cg,config_m,mac,p=False):
    """
    Funcao que retorna a pontuacao da estabilidade estática
    e a posicao do ponto neutro (centro aerodinamico do aviao)
    """
    @executarNaPasta('Runs')
    def getValues():
        """
        Funcao que lê os parâmetros de estabilidade estática
        de o arquivo especificado.
        """
        filename = '%s_trim.stability' %(name)
        f = open(filename, 'r')
        flines = f.readlines()
        for i in range(len(flines)):
            words = flines[i].split()
            if i == 15:
                alpha = float(words[2])
            elif i == 19:                    
                Cltot = float(words[5])
            elif i == 20:
                Cmtot = float(words[5])
            elif i == 21:
                Cntot = float(words[5])
            elif i == 38:
                Clb = float(words[8])*np.pi/180
            elif i == 39:
                Cma = float(words[6])*np.pi/180
            elif i == 40:
                Cnb = float(words[8])*np.pi/180
            elif i == 50:
                Xnp = float(words[4])
            elif i == 52:
                spiral = float(words[6])
        f.close()
        
        return alpha,Xnp,Cmtot,Cma,Cltot,Clb,Cntot,Cnb,spiral
        
    def m_estatica():
        """
        Calculo da margem estatica do aviao.
        """
        Xcg = pos_cg[0]
        return (Xnp-Xcg)/mac
    
    def pontuacao_long():
        """
        Calcula o fator de pontuacao da estabiliade estatica longitudinal.
        """
        print('\nResultados de Estabilidade Estatica Longitudinal:')
        print('         Margem Estatica do aviao : %.5f' %(me))
        print('         Angulo alfa de trimagem :  %.5f graus' %(alfa_trim))
        print('         Cm0 :                      %.5f' %(Cmtot-Cma*alpha))
        print('         Cma :                      %.5f graus^(-1)' %(Cma))        

        max_trim = trim_desejado+2.0

        if alfa_trim > 0.0 and alfa_trim < max_trim:
            f = 1.0
        else:
            if alfa_trim < 0.0:
                f = np.exp(alfa_trim)
            else:
                f = np.exp((max_trim-alfa_trim)/10.0)
        
        if me >= config_m[0] and me <= config_m[1]:
            print('         Os valores sao ACEITAVEIS')
     
        else:
            print('         Os valores NAO sao aceitaveis')
            if me < config_m[0]:
                f *= np.exp(-(config_m[0]-me)/0.05)
            else:
                f *= np.exp(-(me-config_m[1])/0.1)
        print('         Fator de pontuacao : %f' %(f))      
        return f
        
    def pontuacao_lat():
        """
        Calcula o fator de pontuacao da estabilidade estatica lateral.
        """
        print('\nResultados da Estabilidade Estatica Lateral:')
        print('         Clb : %.5f' %(Clb))

        # Valores negativos de Cl_beta sao considerados estaveis
        if Clb < 0.0:
            print '         Aviao lateralmente ESTAVEL'
            f = 1.0
        else:
            print '         Aviao lateralmente INSTAVEL'
            f = np.exp(-Clb)
        print('         Fator de pontuacao : %f' %(f))
        return f
            
    def pontuacao_dir():
        """
        Calcula o fator de pontuacao da estabilidade estatica direcional.
        """
        print('\nResultados da Estabilidade Estatica Direcional:')
        print('         Cnb :  %.5f' %(Cnb))

        # Valores positivos de Cnb sao considerados estaveis
        if Cnb > 0.0:
            print '         Aviao direcionalmente ESTAVEL'
            f = 1.0
        else:
            print '         Aviao direcionalmente INSTAVEL'
            f = np.exp(Cnb)
        print('         Fator de pontuacao : %f' %(f))
        return f
            
    def pontuacao_spiral():
        """
        Calcula o fator de pontuacao da estabilidade espiral.
        """
        print('\nResultados da Estabilidade Espiral:')
        print('         Clb Cnr / Clr Cnb :  %.5f' %(spiral))

        # Valores positivos de spiral sao considerados estaveis
        if spiral > 1.0:
            print '         Modo espiral ESTAVEL'
            f = 1.0
        else:
            print '         Modo espiral INSTAVEL'
            f = np.exp((spiral-1.0)/10.0)
        print('         Fator de pontuacao : %f' %(f))
        return f
            
    def graf_long():
        """
        Gera o grafico do coeficiente de momento em funcao de alfa.
        """
        X = np.linspace(0, alfa_estol)
        Y = Cmtot + Cma*(X-alpha)

        plt.figure()
        plt.grid('on')
        plt.plot(X, Y, 'b-', label=r'$\alpha_{trim}=%.2f$ graus' %(alfa_trim))
        plt.title('Estabilidade Longitudinal da Aeronave')
        plt.xlabel(r'$\alpha$ (graus)')
        plt.ylabel(r'$C_M$')
        plt.legend(loc='best')
        plt.savefig("longitudinal_aviao.png", bbox_inches='tight', dpi=200)
        
    def graf_latdir():
        """
        Gera o grafico de Cl e Cn em funcao de beta.
        """
        X = np.linspace(-10, 10)
        Cl = Cltot + Clb*X
        Cn = Cntot + Cnb*X
        
        # Gera o gráfico do Cd pela área do endplate
        fig, ax1 = plt.subplots()
        plt.title('Estabilidade Latero-direcional da Aeronave')
        ax1.grid('on')
        ax1.plot(X, Cl, 'b-')
        ax1.set_xlabel(r'$\beta$ (graus)')
        ax1.set_ylabel(r'Coeficiente de momento em x ($C_l$)', color='b')   
        for tl in ax1.get_yticklabels():
            tl.set_color('b')
        
        ax2 = plt.twinx()
        ax2.plot(X, Cn, 'r-')
        ax2.set_ylabel(r'Coeficiente de momento em z ($C_n$)', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        plt.savefig("latdir_aviao.png", bbox_inches='tight', dpi=200)        
    
    @executarNaPasta('Graficos/Estabilidade')
    def plot():
        graf_long()
        graf_latdir()
    
    alpha,Xnp,Cmtot,Cma,Cltot,Clb,Cntot,Cnb,spiral = getValues()
    alfa_trim = alpha-Cmtot/Cma
    
    # Analise da estabilidade estática longitudinal, logitudinal e direcional
    me = m_estatica() 
    fator_long = pontuacao_long()
    fator_lat = pontuacao_lat()
    fator_dir = pontuacao_dir()
    fator_spiral = pontuacao_spiral()

    fator_estatica = fator_long*fator_lat*fator_dir*fator_spiral
    
    # Plota os graficos
    if p: plot()
        
    return Xnp,fator_estatica
