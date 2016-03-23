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
        return Asa_frontal.alfa_estol
    
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
                plt.savefig("PolarArrasto.png", bbox_inches='tight', dpi=200)
        
        if p: grafico()        
        return funcao_CDi
    
    def decolagem():
        CL_decol = funcao_CL(alfa_decol)
        CDi_decol = funcao_CDi(CL_decol)
        dec = decol.Analise(S_asaf, peso_vazio, 
                            CL_decol, CLmax, 2*CDi_decol, p)
        CPaga = dec.c_paga_max()
        vd = dec.vel_decol(CPaga)
        
        if p:
            dec.vel_corrida(CPaga)
        
        print("\nCalculando parametros de decolagem para pista de %.1f m" %(dec.c))
        print("         Area de referencia :          %f" %(S_asaf))        
        print("         CL total de corrida :         %f" %(CL_decol))
        print("         CL maximo :                   %f" %(CLmax))
        print("         CDi :                         %f" %(CDi_decol))
        print("         Carga paga maxima do aviao :  %.3f kg" %(CPaga))
        print("         Peso maximo de voo :          %.3f kg" %(CPaga+peso_vazio))
        print("         Velocidade de rotacao :       %.3f m/s" %(vd))
        return CPaga                
    
    alfa_estol = alfa_estol()
    CL,CDi = getCL_CDi()
    funcao_CL,CLmax = polar_sustentacao()
    funcao_CDi = polar_arrasto()
    CPaga = decolagem()
    return CPaga
            
            
            
            
            
            