# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 23:55:59 2015

@author: Diego Chou
"""

print("\nCarregando modulos de 'Aerodinamica'...")
try:
    import numpy as np
    import matplotlib.pyplot as plt
    print("Modulos de 'Aerodinamica' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Aerodinamica'\n")
    raise

import llt
import apoio
import decolagem as decol

def aerodinamica(name,alfa_decol,lista_alfas,peso_vazio,
                 S_asaf,b_asaf,perfil_raiz_asaf,perfil_ponta_asaf,
                 epsilon_asaf,vel,p=False):                                                  
    """
    Gera as polares de sustentação e arrasto do avião.
    """
    
    def sinal(b):
        if b >= 0.0:
            sinal = '+'
        else:
            sinal = '-'
        return sinal     
        
    def alfa_estol():
        """
        O ângulo de estol é aquele em que a asa frontal estola.
        """
        Asa_frontal = llt.S_sustentadora(perfil_raiz_asaf, perfil_ponta_asaf, 
                                         alfa_decol, S_asaf, b_asaf, 1.0,
                                         epsilon_asaf, vel, 12, False, p)
        return Asa_frontal.alfa_estol, Asa_frontal.mac
    
    @apoio.executarNaPasta('Runs')
    def getCL_CDi():
        """
        Função que retorna os pontos do arquivo do perfil.
        O arquivo deve estar no formato aceito pelo XFLR5
        """
        CL = []
        CDi = []
        
        for i in range(len(lista_alfas)):
            filename = '%s_case%d.forces' %(name,i+1)
            f = open(filename, 'r')
            flines = f.readlines()[23:25] 
            for i in range(len(flines)):
                words = flines[i].split()
                val = float(words[2])
                if i == 0:
                    CL.append(val)
                elif i == 1:                    
                    CDi.append(val)
            f.close()
        
        return np.array(CL),np.array(CDi)
        
    def polar_sustentacao():
        angulos = np.append(lista_alfas[0],alfa_estol)
        funcao_CL = np.poly1d(np.polyfit(lista_alfas, CL, deg=1))
        CLmax = funcao_CL(alfa_estol)
        
        @apoio.executarNaPasta('Graficos/Aerodinamica')
        def grafico():
            plt.figure()
            plt.title('Polar de sustentacao ($S_{ref}=%.2fm^2$)' %(S_asaf))
            plt.grid('on')
            plt.plot(lista_alfas,CL,'k.')
            plt.plot(angulos, funcao_CL(angulos), 'b-', 
                     label = r'$C_L=%.3f*\alpha %s %.3f$'%(funcao_CL[1],sinal(funcao_CL[0]),abs(funcao_CL[0])))
            plt.legend(loc='best')
            plt.xlabel(r'$\alpha$ ($graus$)')
            plt.ylabel(r'$C_L$')
            plt.savefig("PolarSustentacao.png", bbox_inches='tight', dpi=200)
            
        if p: grafico()
        return funcao_CL,CLmax
        
    def polar_arrasto():
        CLs = np.linspace(CL[0],CLmax)
        funcao_CDi = np.poly1d(np.polyfit(CL, CDi, deg=2))
        
        @apoio.executarNaPasta('Graficos/Aerodinamica')
        def grafico():
                if abs(funcao_CDi[1]) <= 1e-3: 
                    string = r'$C_D=%.3f+%.3f*C_L^2$'%(funcao_CDi[0],
                                                       funcao_CDi[2])     
                else:
                    string = r'$C_D=%.3f %s %.3f*C_L+%.3f*C_L^2$'%(funcao_CDi[0],
                                                                   sinal(funcao_CDi[1]),
                                                                   abs(funcao_CDi[1]),funcao_CDi[2])                
                plt.figure()
                plt.grid('on')
                plt.title('Polar de arrasto induzido da aeronave ($S_{ref}=%.2fm^2$)' %(S_asaf))
                plt.plot(CL,CDi,'k.')                
                plt.plot(CLs, funcao_CDi(CLs), 'b-', label = string)
                plt.legend(loc='best')
                plt.xlabel('$C_L$')
                plt.ylabel('$C_{D_i}$')
                plt.savefig("PolarArrastoInduzido.png", bbox_inches='tight', dpi=200)
        
        if p: grafico()        
        return funcao_CDi
    
    def decolagem():
        m_a = 1.5*peso_vazio
        CL_decol = funcao_CL(alfa_decol)
        CDi_decol = funcao_CDi(CL_decol)
        CD_decol = 2*CDi_decol
        dec = decol.Analise(S_asaf, m_a, 
                            CL_decol, CLmax, CD_decol, p)
        CPaga = dec.c_paga_max()
        vd = dec.vel_decol(CPaga)
        
        def mapa(tipo=1,n=50):
            """
            Essa funçao gera um gráfico de curvas de nível para a avaliação do efeito da 
            variaçao dos parametros de CL e CD da aeroanve para a decolagem.
            Essa análise é importante para a decisão de se utilizar ou não dispositivos hipersustentadores.
            """
            if tipo == 1:
                x = np.linspace(0.8*CLmax,1.2*CLmax,n)
                y = np.linspace(0.8*CD_decol,1.2*CD_decol,n)
            elif tipo == 2:
                x = np.linspace(0.8*CL_decol,1.2*CL_decol,n)
                y = np.linspace(0.8*CD_decol,1.2*CD_decol,n)
            elif tipo == 3:
                x = np.linspace(0.8*CLmax,1.2*CLmax,n)
                y = np.linspace(0.8*CL_decol,1.2*CL_decol,n)
            
            X, Y = np.meshgrid(x,y)
            
            Z = np.zeros([n,n])
            
            with apoio.esconderPrint():
                for i in range(n):
                    for j in range(n):
                        if tipo == 1:
                            element = decol.Analise(S_asaf, m_a, CL_decol, X[i,j], Y[i,j], p=False)
                        elif tipo == 2:
                            element = decol.Analise(S_asaf, m_a, X[i,j], CLmax, Y[i,j], p=False)
                        elif tipo == 3:
                            element = decol.Analise(S_asaf, m_a, Y[i,j], X[i,j], CD_decol, p=False)
                        
                        Z[i,j] = element.c_paga_max()
            
            @apoio.executarNaPasta('Graficos/Desempenho')
            def grafico():
                plt.figure()
                plt.grid('on')
#                plt.axis('equal')
                plt.title('Carga Paga Maxima ($kg$)')
                if tipo == 1:
                    plt.xlabel('$\Delta C_{L_{max}}$')
                    plt.ylabel('$\Delta C_D$')
                    CS = plt.contour(X-CLmax, Y-CD_decol, Z, colors='k')
                    plt.clabel(CS, inline=1, fontsize=10)
                    plt.savefig("mapaDecolagem1.png", bbox_inches='tight', dpi=200) 
                elif tipo == 2:
                    plt.xlabel('$\Delta C_L$')
                    plt.ylabel('$\Delta C_D$')
                    CS = plt.contour(X-CL_decol, Y-CD_decol, Z, colors='k')
                    plt.clabel(CS, inline=1, fontsize=10)
                    plt.savefig("mapaDecolagem2.png", bbox_inches='tight', dpi=200)
                elif tipo == 3:
                    plt.xlabel('$\Delta C_{L_{max}}$')
                    plt.ylabel('$\Delta C_L$')
                    CS = plt.contour(X-CLmax, Y-CL_decol, Z, colors='k')
                    plt.clabel(CS, inline=1, fontsize=10)
                    plt.savefig("mapaDecolagem3.png", bbox_inches='tight', dpi=200)
             
            grafico()   
        
        if p:
            dec.vel_corrida(CPaga)
            mapa(1)
            mapa(2)
            mapa(3)  
        
        print("\nCalculando parametros de decolagem para pista de %.1f m" %(dec.c))
        print("         Area de referencia :          %f m^2" %(S_asaf))        
        print("         CL total de corrida :         %f" %(CL_decol))
        print("         CL maximo :                   %f" %(CLmax))
        print("         CDi :                         %f" %(CDi_decol))
        print("         Carga paga maxima do aviao :  %.3f kg" %(CPaga))
        print("         Peso maximo de voo :          %.3f kg" %(CPaga+m_a))
        print("         Velocidade de rotacao :       %.3f m/s" %(vd))
        return CPaga                
    
    alfa_estol,mac = alfa_estol()
    CL,CDi = getCL_CDi()
    funcao_CL,CLmax = polar_sustentacao()
    funcao_CDi = polar_arrasto()
    CPaga = decolagem()
    return CPaga,alfa_estol,mac
            
            
            
            
            
            