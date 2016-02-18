# -*- coding: utf-8 -*-
"""
Created on Sun Feb 08 12:32:47 2015

@author: Diego Chou
"""
print("\nCarregando modulos de 'Avaliador'...")
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import time
    print("Modulos de 'Avaliador' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Avaliador'\n")
    raise

import construtor as con
import apoio
import perfil
import aerodinamica as aero
import estabilidade as est
import fuselagem
import aeronave3D
"""
Esse arquivo define os avaliadores primários, ou seja,
classes que realizam uma análise para ser iterada com 
o algoritmo genético, afim de alcançar uma configuração ótima, 
ou próxima dela.
"""

class AvaliadorLivre(con.ConstrutorLivre):
    """
    Essa classe avalia o indivíduo
    """
    #Leitura do arquivo de configuração
    config = apoio.readConfFile('avaliador.ini')
    
    def __init__(self,nome,p=False):
        super(AvaliadorLivre,self).__init__(nome)
        self.analise(p)
    
    def analise(self,p):
        # Consultar longitudinal.py
        self.i_b = np.arctan((self.h_cauda-(self.h_trem+self.hccarga/2.0+self.d_roda/2.0))/self.l_boom)*180/np.pi
        self.i_fb = np.arctan(0.5*(self.hccarga-fuselagem.Parede_fogo().h)/self.comp_bico)*180/np.pi
        
        args_aero = [self.m,self.helice,self.alfa,self.vel,self.peso_total,self.xcg,self.d_cg_raiz,
                     self.perfilw_raiz,self.perfilw_ponta,self.iw,self.Sw,self.Seffw,
                     self.bw,self.afilw,self.epsilonw,self.perfilEH,self.Seh,self.afileh,
                     self.AReh,self.lt,self.vht,self.perfilEV,self.Sev,self.afilev,
                     self.ARev,self.lt,self.iev,self.vvt,self.h_cauda,self.esp_boom, 
                     self.l_boom,self.comp_bico+self.comp_ccarga,self.hccarga,
                     self.largccarga,self.Vol_ccarga,self.i_b,self.i_fb,self.comp_bico,
                     self.angulo_base_bico,self.h_trem,self.d_roda,self.esp_roda, p]
        
        # Calculo aerodinâmico
        self.CPaga, self.vdecol, self.CL, self.aw, self.enflex_ca,\
        self.beh,self.Creh,self.ieh,self.bev,self.Crev,\
        args_est = aero.aerodinamica(*args_aero)   
        
        #Modelo 3D para calculo do momento de inercia
#        self.desenho3D(p,True)
        
        #Calculo de estabilidade
        self.qualidade_est = est.estabilidade(*args_est)
        
        #Define a pontuação avaliada para o indivíduo
        self.pontuacao = self.CPaga*self.qualidade_est
        
        if p:
            self.vista()
        self.desenho3D(p)
   
    def vista(self):
        """
        Vista lateral e superior do avião, sem motor.
        """
        # Rodando! Girando e rodando! Maôe!
        def gira(listax, listay, ponto, alfa):
            """
            Gira os pontos definidos pelas listas em função
            de um ponto.
            """
            alfa = alfa*np.pi/180
            
            # Posição relativa ao ponto de referência
            n_listax = np.array(listax)-ponto[0]
            n_listay = np.array(listay)-ponto[1]
            
            # Matriz de rotação horária
            M_rot = np.matrix([[np.cos(alfa),np.sin(alfa)],[-np.sin(alfa),np.cos(alfa)]])            
            
            #Supõe-se que listax e listay são do mesmo tamanho
            for i in range(len(listax)):
                vetor = np.array([[n_listax[i]],[n_listay[i]]])
                vetor = M_rot*vetor
                listax[i] = float(vetor[0][0])+ponto[0]
                listay[i] = float(vetor[1][0])+ponto[1]
            return listax, listay
            
        def pos_lateral():
            """
            Define as posições para que a vista lateral seja desenhada.
            """
            #Definição do trem de pouso principal
            r_roda = self.d_roda/2.0
            h_trem = self.h_trem
            pos_trem = self.pos_trem
            ponto_trem = [pos_trem,r_roda]
            trem = plt.Circle(ponto_trem,r_roda,color='black',fill=False)
            
            #Definição da bequilha
            pos_bequilha = -self.pos_bequilha
            ponto_bequilha = [pos_bequilha,r_roda]
            beq = plt.Circle(ponto_bequilha,r_roda,color='black',fill=False)         
            
            #Definição da asa
            h_fus = self.hccarga
            corda = self.raiz
            corda_antesCG = self.d_cg_raiz
            corda_depoisCG = corda-corda_antesCG 
            xraiz,yraiz = perfil.Analise(self.perfilw_raiz).getPoints() 
            xraiz = corda*np.array(xraiz)-corda_antesCG
            yraiz = corda*np.array(yraiz)+h_fus+h_trem+r_roda
            xraiz,yraiz = gira(xraiz,yraiz,[-corda_antesCG,h_fus+h_trem+r_roda],self.iw)
            xraiz,yraiz = gira(xraiz,yraiz,ponto_trem,self.alfa)
            
            #Pontos que definem a fuselagem
            parede_fogo = fuselagem.Parede_fogo()
            h_bico = parede_fogo.h
            l_fus = self.comp_ccarga
            l_bico = self.comp_bico
            l_boom = self.l_boom
            h_cauda = self.h_cauda
            
            val = [l_fus/2.0,corda_depoisCG]
            
            xfus = [-l_fus/2.0,-l_fus/2.0,max(val),l_fus/2.0,-l_fus/2.0]
            yfus = [h_trem+r_roda,h_trem+h_fus+r_roda,h_trem+h_fus+r_roda,h_trem+r_roda,h_trem+r_roda]
            xfus,yfus = gira(xfus,yfus,ponto_trem,self.alfa)
            
            xbico = [-l_fus/2.0,-l_fus/2.0-l_bico,-l_fus/2.0-l_bico,-l_fus/2.0] 
            ybico = [h_trem+r_roda,h_trem+(h_fus-h_bico)+r_roda,h_trem+h_fus+r_roda,h_trem+h_fus+r_roda]
            xbico,ybico = gira(xbico,ybico,ponto_trem,self.alfa)
            
            #Definição da Eh
            Creh = self.Creh
            xeh,yeh = perfil.Analise(self.perfilEH).getPoints()
            xeh = Creh*np.array(xeh)+l_boom+l_fus/2.0-0.25*Creh
            yeh = (-1)*Creh*np.array(yeh)+h_cauda
            xeh,yeh = gira(xeh,yeh,[l_fus/2.0+l_boom,h_cauda],self.ieh)
            xeh,yeh = gira(xeh,yeh,ponto_trem,self.alfa)
            
            #Definição da EV
            afilev = self.afilev
            bev = self.bev
            Crev = self.Crev
            Ctev = Crev*afilev
            d = l_fus/2.0+l_boom
            xev = [d-0.25*Crev,d-0.25*Ctev,d+0.75*Ctev,d+0.75*Crev,d-0.25*Crev]
            yev = [h_cauda,h_cauda+bev,h_cauda+bev,h_cauda,h_cauda]
            xev,yev = gira(xev,yev,ponto_trem,self.alfa)
            
            # Definição do boom 
            xboom = [max(val),l_fus/2.0+l_boom,l_fus/2.0]
            yboom = [h_trem+h_fus+r_roda,h_cauda,h_trem+r_roda]
            xboom,yboom = gira(xboom,yboom,ponto_trem,self.alfa)            
            
            return trem,beq,xraiz,yraiz,xfus,yfus,xbico,ybico,xboom,yboom,\
                   xeh,yeh,xev,yev
                   
        @apoio.executarNaPasta('Graficos/Vistas')
        def vista_lateral():
            plt.figure()
            fig = plt.gcf()
            plt.axis('equal') 
            plt.grid('on') 
            fig.gca().add_artist(trem)
            fig.gca().add_artist(beq)
            plt.plot(xraiz,yraiz,'k-')
            plt.plot(xfus,yfus,'k-')
            plt.plot(xbico,ybico,'k-')
            plt.plot(xboom,yboom,'k-')
            plt.plot(xeh,yeh,'k-')
            plt.plot(xev,yev,'k-')
            #plt.plot(x__,y__,'k-')
        
            #Salva a figura
            plt.title(r'Vista lateral da aeronave, $\alpha=%.1f$' %(self.alfa))
            plt.savefig("vista_lateral.png", bbox_inches='tight', dpi=200)
        
        def pos_superior():
            """
            Define as posições para que a vista superior seja desenhada.
            """
            #Definição da asa
            bw = self.bw
            raiz = self.raiz 
            pos_trans = self.pos_trans 
            ponta = self.ponta
            offset = self.offset_ponta
            mac = self.mac
            pos_mac = self.pos_mac
            raiz = self.raiz
            
            tpw = self.bw/2.0-pos_trans # semi-envergadura em que a asa é trapezoidal
            offset_mac = (pos_mac-pos_trans)*offset/(tpw)
            
            raiz_antesCG = self.d_cg_raiz
            raiz_depoisCG = raiz-raiz_antesCG 
            ponta_antesCG = raiz_antesCG-offset
            ponta_depoisCG = ponta-ponta_antesCG
            mac_antesCG = raiz_antesCG-offset_mac
            mac_depoisCG = mac-mac_antesCG
            
            xasa = [pos_trans,bw/2.0,bw/2.0,pos_trans,
                     -pos_trans,-bw/2.0,-bw/2.0,-pos_trans, pos_trans]
            yasa = [-raiz_antesCG,-ponta_antesCG,
                     ponta_depoisCG,raiz_depoisCG,
                     raiz_depoisCG,ponta_depoisCG,
                     -ponta_antesCG,-raiz_antesCG,-raiz_antesCG]
            xmac = [pos_mac,pos_mac]
            ymac = [-mac_antesCG,mac_depoisCG]
            
            #Definição da fuselagem
            parede_fogo = fuselagem.Parede_fogo()
            larg_bico = parede_fogo.larg
            l_fus = self.comp_ccarga
            larg_fus = self.largccarga
            l_bico = self.comp_bico
            l_boom = self.l_boom
            
            val = [l_fus/2.0,raiz_depoisCG]    
            
            xfus = [larg_fus/2.0,larg_fus/2.0,-larg_fus/2.0,-larg_fus/2.0,larg_fus/2.0]
            yfus = [max(val),-l_fus/2.0,-l_fus/2.0,max(val),max(val)]
            
            xbico = [larg_fus/2.0,larg_bico/2.0,-larg_bico/2.0,-larg_fus/2.0] 
            ybico = [-l_fus/2.0,-l_fus/2.0-l_bico,-l_fus/2.0-l_bico,-l_fus/2.0]
        
            #Definição da EH
            afileh = self.afileh
            beh = self.beh
            Creh = self.Creh
            Cteh = Creh*afileh
            d = l_fus/2.0+l_boom
            xeh = [0.0,beh/2.0,beh/2.0,0.0,-beh/2.0,-beh/2.0,0.0]
            yeh = [d-0.25*Creh,d-0.25*Cteh,d+0.75*Cteh,d+0.75*Creh,
                   d+0.75*Cteh,d-0.25*Cteh,d-0.25*Creh]
            
            #Definição da EV
            Crev = self.Crev
            yev,xev = perfil.Analise(self.perfilEV).getPoints()
            yev = Crev*np.array(yev)+l_boom+l_fus/2.0-0.25*Crev
            xev = Crev*np.array(xev)
    
            #Definição do boom
            xboom_d = [larg_fus/2.0 for i in range(4)]
            xboom_e = [-larg_fus/2.0 for i in range(4)]
            yboom = [max(val),l_fus/2.0+l_boom-larg_fus-0.25*Crev,
                     l_fus/2.0+l_boom-0.25*Crev,l_fus/2.0+l_boom]            
            
            return xasa,yasa,xmac,ymac,xfus,yfus,xbico,ybico,\
                   xboom_d,xboom_e,yboom,xeh,yeh,xev,yev
                   
        @apoio.executarNaPasta('Graficos/Vistas')
        def vista_superior():
            plt.figure()
            plt.axis('equal') 
            plt.grid('on')
            plt.plot(xasa,yasa,'k-')
            plt.plot(xmac,ymac,'b-',label='CMA da asa')
            plt.plot(xfus_sup,yfus_sup,'k-')
            plt.plot(xbico_sup,ybico_sup,'k-')
            plt.plot(xboom_d,yboom_sup,'k-')
            plt.plot(xboom_e,yboom_sup,'k-')
            plt.plot(xeh_sup,yeh_sup,'k-')
            plt.plot(xev_sup,yev_sup,'k-')
            
            #Salva a figura
            plt.legend(loc='best')
            plt.title('Vista superior da aeronave')
            plt.savefig("vista_superior.png", bbox_inches='tight', dpi=200)    
        
        trem,beq,xraiz,yraiz,xfus,yfus,xbico,ybico,xboom,yboom,\
        xeh,yeh,xev,yev = pos_lateral()
        
        xasa,yasa,xmac,ymac,xfus_sup,yfus_sup,\
        xbico_sup,ybico_sup,xboom_d,xboom_e,yboom_sup,\
        xeh_sup,yeh_sup,xev_sup,yev_sup = pos_superior()
        
        vista_lateral()
        vista_superior()
    
    def desenho3D(self, p = False, inercia = False):
        """
        Essa função escreve um arquivo .stl que é lido pelo script vtk3D.py
        para gerar um desenho 3D do avião.
        """
        perfilw_raiz = self.perfilw_raiz
        perfilw_ponta = self.perfilw_ponta
        perfilEH = self.perfilEH
        perfilEV = self.perfilEV
        
        iw = self.iw
        ieh = self.ieh
        iev = self.iev
        
        #ASA
        raiz = self.raiz
        ponta = self.ponta
        mac = self.mac
        d_cg_raiz = self.d_cg_raiz
        dist_trans = self.pos_trans
        semi_bw = self.bw/2.0
        washout = self.epsilonw 
        offset = self.offset_ponta
        
        #EH
        raiz_eh = self.Creh
        semi_beh = self.beh/2.0 
        afil_eh = self.afileh
        
        #EV
        raiz_ev = self.Crev
        bev = self.bev
        afil_ev = self.afilev

        #FUSELAGEM        
        parede = fuselagem.Parede_fogo()          
        h_bico = parede.h
        larg_bico = parede.larg
        comp_bico = self.comp_bico 
        h_ccarga = self.hccarga
        larg_ccarga = self.largccarga
        comp_baixo_ccarga = self.comp_ccarga
        comp_cima_ccarga = comp_baixo_ccarga/2.0 + max([comp_baixo_ccarga/2.0,raiz-d_cg_raiz])
        
        #TAILBOOM
        comp_boom = self.l_boom
        esp_boom = self.esp_boom
        h_cauda_rel_asa = self.h_cauda-(self.h_trem+self.hccarga+self.d_roda/2.0)
        
        args_asa = (perfilw_raiz, perfilw_ponta, iw, raiz, ponta, mac, d_cg_raiz, 
                    dist_trans, semi_bw, washout, offset)
        
        args_eh = (perfilEH, ieh, raiz_eh,raiz_eh*afil_eh, (1+afil_eh)*raiz_eh/2.0,
                   0.0,semi_beh, 0.0, 0.25*(1.0-afil_eh)*raiz_eh)
                   
        args_ev = (perfilEV, iev, raiz_ev,raiz_ev*afil_ev, (1+afil_ev)*raiz_ev/2.0,
                   0.0,bev, 0.0, 0.25*(1.0-afil_ev)*raiz_ev)
                   
        args_fus = (h_bico, larg_bico, comp_bico,
                    h_ccarga, larg_ccarga, comp_cima_ccarga,
                    comp_baixo_ccarga)
                    
        args_boom = (comp_boom, esp_boom, h_cauda_rel_asa)
        
        if inercia:
            aeronave3D.main('aeronave',args_asa,args_eh,args_ev,args_fus, args_boom,inercia)

        else: 
            aeronave3D.main('ind_score_%.6f' %(self.pontuacao), args_asa,
                        args_eh,args_ev,args_fus, args_boom)
            if p:
        	    import vtk3D
        	    vtk3D.main('ind_score_%.6f.stl' %(self.pontuacao))
    
    @apoio.executarNaPasta('Aeronaves')     
    def salvarParametros(self,arquivo):
        """
        Função que escreve o arquivo .ini com os parâmetros do avião.
        Foi feita para ser rodada apenas pelas classes filhas da classe AvaliadorLivre.
        """        
        if self.__class__.__name__ == "AvaliadorLivre":
            raise IOError("O metodo salvarParametros nao deve ser instanciado por AvaliadorLivre")            
            
        # Abre o arquivo .ini
        arquivo = open(str(arquivo)+'.ini', 'w')
        
        #escreve as partes do arquivo
        arquivo.write('[Analise]\n')
        arquivo.write('velocidade : %f\n' %(self.vel))
        arquivo.write('alfa : %f\n' %(self.alfa))
        arquivo.write('\n')
        arquivo.write('[Miscelanea]\n')
        arquivo.write('xcg : %f\n' %(self.xcg))
        arquivo.write('motor: %s\n' %(self.m))
        arquivo.write('helice: %s\n' %(self.helice))
        arquivo.write('\n')
        arquivo.write('[Asa]\n')
        arquivo.write('perfil_raiz: %s\n' %(self.perfilw_raiz))
        arquivo.write('perfil_ponta: %s\n' %(self.perfilw_ponta))   
        arquivo.write('area_planta: %f\n' %(self.Sw))
        arquivo.write('envergadura: %f\n' %(self.bw))
        arquivo.write('afilamento: %f\n' %(self.afilw))
        arquivo.write('incidencia: %f\n' %(self.iw))
        arquivo.write('washout: %f\n' %(self.epsilonw))
        arquivo.write('\n')
        arquivo.write('[EH]\n')
        arquivo.write('perfil: %s\n' %(self.perfilEH))
        arquivo.write('area_planta: %f\n' %(self.Seh))
        arquivo.write('afilamento: %f\n' %(self.afileh))            
        arquivo.write('alongamento: %f\n' %(self.AReh))
        arquivo.write('densidade: %f\n' %(self.k_eh))
        arquivo.write('\n')
        arquivo.write('[EV]\n')
        arquivo.write('perfil: %s\n' %(self.perfilEV))
        arquivo.write('area_planta: %f\n' %(self.Sev))
        arquivo.write('afilamento: %f\n' %(self.afilev))            
        arquivo.write('alongamento: %f\n' %(self.ARev))
        arquivo.write('densidade: %f\n' %(self.k_ev))
        arquivo.write('\n')
        arquivo.write('[Compartimento Carga]\n')
        arquivo.write('largura: %f\n' %(self.largccarga))            
        arquivo.write('altura: %f\n' %(self.hccarga))
        arquivo.write('comprimento: %f\n' %(self.comp_ccarga))
        arquivo.write('\n')
        arquivo.write('[Bico]\n')
        arquivo.write('afilamento: %f\n' %(self.ang_afil_bico))
        arquivo.write('\n')
        arquivo.write('[Tailboom]\n')
        arquivo.write('tipo: trelicado\n')
        arquivo.write('quantidade: %d\n' %(self.n_boom))
        arquivo.write('densidade: %f\n' %(self.k_boom))
        arquivo.write('espessura: %f\n' %(self.esp_boom))
        arquivo.write('\n')
        arquivo.write('[Trem de pouso]\n')
        arquivo.write('angulo: %f\n' %(self.ang_trem))
        arquivo.write('altura: %f\n' %(self.h_trem))
        arquivo.write('diametro_roda: %f\n' %(self.d_roda))     
        arquivo.write('espessura_roda: %f\n' %(self.esp_roda))
        arquivo.write('\n')
        arquivo.write('[Cauda]\n')
        arquivo.write('lt: %f\n' %(self.lt))
        arquivo.write('h_cauda: %f\n' %(self.h_cauda))
        
        #fecha o arquivo
        arquivo.close()

################################################################################

class Avaliador2015(con.Construtor2015,AvaliadorLivre):
    """
    Essa classe avalia o indivíduo
    """
    #Leitura do arquivo de configuração
    config = apoio.readConfFile('avaliador.ini')
    
    def __init__(self, vel, alfa, xcg,                  # Variáveis gerais
                 iw, Sw, bw, afilw,                     # Variáveis da asa 
                 AReh, ARev, vvt_vht,                   # Variáveis do estabilizador horizontal e vertical
                 hccarga, largccarga,                   # Variáveis do compartimento de carga 
                 p=False, out=False):             
        con.Construtor2015.__init__(self,vel, alfa, xcg,            
                                    iw, Sw, bw, afilw,               
                                    AReh, ARev, vvt_vht,
                                    hccarga, largccarga)                              
        self.analise(p)

        if out:   
            string = 'output_avaliador2015_%s_%s' %(time.strftime("%Y %m %d"), time.strftime("%H.%M.%S"))
            self.salvarParametros(string)
