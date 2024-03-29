# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 21:05:25 2016

@author: Diego Chou
"""

print("\nCarregando modulos de 'Construtor'...")

try:
    import numpy as np
    print("Modulos de 'Construtor' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Construtor'\n")
    raise

import apoio
import write_input_file as wif
import perfil
import geral

"""
 Todos os pesos estão em kg. 
 Todas as distâncias em metros.
 Todos as funções retornam angulos em graus.
"""

class Construtor2016(object):
    """
    Classe responsável por gerar a aeronave, com os critérios estabelecidos
    pelo regulamento do ano vigente. O regulamento de 2016 estabeleceu 
    hangaragem em um cone de 2.5m de diâmetro e 0.75m de altura.
    
    As coordenadas estão dispostas no padrão do software Athena Vortex Lattice.
    Por simplificaçao inicial de projeto, a asa é retangular, não possui 
    enflexamento nem diedro.    
    
    
    INPUTS:
    ------
    * dz_asas: Altura entre as asas
    * asaf: Asa frontal
    * asat: Asa traseira
    * c: Corda
    * ang: Angulo de incidencia
    * epsilon: Angulo de torcao
    """
    
    d_cone = 2.9 # Diâmetro do cone
    h_cone = 0.75 # Altura do cone
    m_motor = 0.731 # Peso de motor e hélice
    z_motor = 0.2 
    k_asa = [0.41716, 0.5] # Constantes para o calculo da densidade da asa 
    ro_estabilizador = 0.9921 #Densidade do EV
    
    # Parametros gerométricos
    x_min = 0.330-0.5*d_cone # Ponto mais a frente da aeronave
    z_min = 0.07 # Altura do ponto mais baixo da aeronave    
    
    #Carga    
    m_carga = 12 # Carga projetada em 2015
    ro_carga = 8.8e3 #Densidade do bronze
    
    #Parametros de discretizaçao das asas
    Nchord = 8
    Cspace = 1.0
    Nspan = 16
    Sspace = -2.0
    
    def __init__(self,name,dz_asas, vel, x_motor,
                 x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                 x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                 c_ev,perfil_ev,p):    
        # Parametros gerais iteráveis
        self.name = name
        self.dz_asas = dz_asas 
        self.vel = vel
        self.p = p
        self.g = geral.gravidade
        
        # Ar
        self.ar = geral.Ar()
        self.mi_ar = self.ar.mi
        self.ro_ar = self.ar.rho()
        self.vel_som = geral.vel_som
        
        #Motor
        if x_motor > x_ba_asaf+0.063:
            raise ValueError("O motor colide com a asa frontal")
        elif x_motor < self.x_min+0.063:
            raise ValueError("O motor foi colocado muito a frente")
        self.pos_motor = np.array([x_motor, 0.0, self.z_motor])
        
        # Parametros de asa frontal 
        self.x_ba_asaf = x_ba_asaf
        self.c_asaf = c_asaf
        self.ang_asaf = ang_asaf
        self.epsilon_asaf = epsilon_asaf
        self.perfilr_asaf = str(perfilr_asaf)
        self.perfilp_asaf = str(perfilp_asaf)
        
        #Parametros de asa traseira
        self.x_bf_asat = x_bf_asat
        self.c_asat = c_asat
        self.ang_asat = ang_asat
        self.epsilon_asat = epsilon_asat
        self.perfilr_asat = str(perfilr_asat)
        self.perfilp_asat = str(perfilp_asat)       
        
        # Estabilizador Vertical
        self.c_ev = c_ev
        self.perfil_ev = perfil_ev
        
        #Calculos para a asa frontal
        self.pos_ba_cr_asaf, self.pos_bf_cr_asaf, self.ang_cr_asaf,\
        self.pos_ba_ct_asaf, self.pos_bf_ct_asaf, self.ang_ct_asaf,\
        self.S_asaf, self.bw_asaf, self.AR_asaf = self.asa_frontal()
        
        #Calculos para a asa traseira 
        self.pos_ba_cr_asat, self.pos_bf_cr_asat, self.ang_cr_asat,\
        self.pos_ba_ct_asat, self.pos_bf_ct_asat, self.ang_ct_asat,\
        self.S_asat, self.bw_asat, self.AR_asat = self.asa_traseira()
        
        #Calculos para o EV
        self.pos_ba_cr_ev, self.pos_bf_cr_ev,\
        self.pos_ba_ct_ev, self.pos_bf_ct_ev,\
        self.S_ev, self.b_ev, self.AR_ev = self.ev()
        
        #Posição do CG levando em conta pesos do motor e asas
        self.pos_cg, self.m_vazio = self.cg()
        
        #Tensor de inercia do aviao
        self.J = self.inercia()
        
        # Escreve arquivo '.avl'
        header, asas = self.input_val(self.Nchord,self.Cspace,self.Nspan,self.Sspace)
        wif.Input_geometry(self.name,header,asas)
        
        #Escreve arquivo '.mass'
        inertia = [self.m_vazio+self.m_carga,
                   self.pos_cg[0],self.pos_cg[1],self.pos_cg[2],
                   self.J[0,0],self.J[1,1],self.J[2,2],
                   self.J[0,1],self.J[0,2],self.J[1,2]]
        wif.Input_mass(name,self.g,self.ro_ar,inertia)
        
    def densidade_asas(self,AR):
        """
        Retorna a densidade por área de asa
        """
        return self.k_asa[0]*AR**self.k_asa[1]
        
    def cone(self,a,z):
        """
        Função que retorna as coordenadas do cone
        """
        z0 = self.h_cone
        c = self.d_cone/(2*self.h_cone)
        if z > z0:
            raise ValueError("Altura máxima ultrapassada")
        R = c*(z0-z)
        b1 = np.sqrt(R**2-a**2)
        b2 = -b1
        return b1,b2

    def cone2(self,x,y):
        z0 = self.h_cone
        c = self.d_cone/(2.*z0)
        z = -np.sqrt(x**2+y**2)/c+z0
        if z < 0.0:
            raise ValueError("Posicao fora dos limites do cone")
        return z
    
    def asa_frontal(self):
        """
        Geometria da asa frontal
        """
        if self.x_ba_asaf < self.x_min:
            raise ValueError("A asa frontal foi colocada muito a frente")
        corda = self.c_asaf
        iw = self.ang_asaf
        epsilon = self.epsilon_asaf
        
        def pos_cr():
            """
            Retorna o ponto do bordo de ataque e a incidencia
            da corda raiz
            """
            y = 0.0
            x_ba = self.x_ba_asaf
            x_bf = x_ba+corda*np.cos(iw*np.pi/180)
            z_bf = self.z_min            
            z_ba = z_bf+corda*np.sin(iw*np.pi/180)
            ang_cr = iw
            
            return np.array([x_ba,y,z_ba]), np.array([x_bf,y,z_bf]), ang_cr
            
        def pos_ct():
            """
            Retorna o ponto do bordo de ataque e a incidencia
            da corda da ponta
            """
            z0 = self.h_cone
            c = self.d_cone/(2.*z0)
            z_ba = pos_ba_cr[2]
            ang_ct = iw+epsilon            
            
            def calc(string):
                """
                Calculo das posições (x,y)
                
                STRING:
                ------
                * 'ba': Em relação ao bordo de ataque
                * 'bf': Em relação ao bordo de fuga
                """
                R = c*(z0-z_ba)
                if string == 'ba':
                    y = np.sqrt(R**2-pos_ba_cr[0]**2) #float(max(Y))
                elif string == 'bf':
                    y = np.sqrt(R**2-pos_bf_cr[0]**2)
                return y
                
            y = calc('ba')
            
            # Se o valor de y para o ponto do bordo de fuga estiver fora do cone
            # a função calc é refeita em função do bordo de fuga.
            if y > max(self.cone(pos_bf_cr[0],z_ba)):
                y = calc('bf')    
                
            return np.array([pos_ba_cr[0],y,pos_ba_cr[2]]),\
                   np.array([pos_bf_cr[0],y,pos_bf_cr[2]]), ang_ct
            
        def Sw():
            S1 = np.cross(pos_bf_cr-pos_ba_cr,pos_ba_ct-pos_ba_cr)
            S2 = np.cross(pos_ba_ct-pos_bf_ct,pos_bf_cr-pos_bf_ct)
            return np.linalg.norm(S1+S2)
            
        def AR():
            bw = 2*pos_ba_ct[1]
            return bw,bw**2/Sw 
        
        pos_ba_cr, pos_bf_cr, ang_cr = pos_cr()
        pos_ba_ct, pos_bf_ct, ang_ct = pos_ct()
            
        Sw = Sw()
        bw,AR = AR()  
        
        return pos_ba_cr, pos_bf_cr, ang_cr, pos_ba_ct, pos_bf_ct, ang_ct, Sw, bw, AR
    
    def asa_traseira(self):
        """
        Geometria da asa traseira
        """
        corda = self.c_asat
        iw = self.ang_asat
        epsilon = self.epsilon_asat        
        
        def pos_cr():
            """
            Retorna o ponto do bordo de ataque e a incidencia
            da corda raiz
            """
            y = 0.0
            z_ba = self.pos_ba_cr_asaf[2]+self.dz_asas
            z_bf = z_ba - corda*np.sin(iw*np.pi/180)

            if min([z_ba,z_bf]) < self.z_min:
                raise ValueError("A asa traseira foi colocada muito baixa")
            
            x_bf = self.x_bf_asat
            x_ba = x_bf-corda*np.cos(iw*np.pi/180)
            ang_cr = iw
            
            x_max = 0.5*self.d_cone*(1 - z_ba/self.h_cone)
            
            if self.x_bf_asat > x_max:
                raise ValueError("A asa traseira foi colocada muito atras")
            
            return np.array([x_ba,y,z_ba]), np.array([x_bf,y,z_bf]), ang_cr
            
        def pos_ct():
            """
            Retorna o ponto do bordo de ataque e a incidencia
            da corda da ponta
            """
            z0 = self.h_cone
            c = self.d_cone/(2.*z0)
            z_ba = pos_ba_cr[2]
            ang_ct = iw+epsilon            
            
            def calc(string):
                """
                Calculo das posições (x,y)
                
                STRING:
                ------
                * 'ba': Em relação ao bordo de ataque
                * 'bf': Em relação ao bordo de fuga
                """
                R = c*(z0-z_ba)
                if string == 'ba':
                    y = np.sqrt(R**2-pos_ba_cr[0]**2) #float(max(Y))

                elif string == 'bf':
                    y = np.sqrt(R**2-pos_bf_cr[0]**2)
                
                return y
                
            y = calc('bf')
            
            # Se o valor de y para o ponto do bordo de fuga estiver fora do cone
            # a função calc é refeita em função do bordo de fuga.
            if y > max(self.cone(pos_ba_cr[0],z_ba)):
                y = calc('ba')    
            
            return np.array([pos_ba_cr[0],y,pos_ba_cr[2]]),\
                   np.array([pos_bf_cr[0],y,pos_bf_cr[2]]), ang_ct
        
        def Sw():
            S1 = np.cross(pos_bf_cr-pos_ba_cr,pos_ba_ct-pos_ba_cr)
            S2 = np.cross(pos_ba_ct-pos_bf_ct,pos_bf_cr-pos_bf_ct)
            return np.linalg.norm(S1+S2)
            
        def AR():
            bw = 2*pos_ba_ct[1]
            return bw, bw**2/Sw        
        
        pos_ba_cr, pos_bf_cr, ang_cr = pos_cr()
        pos_ba_ct, pos_bf_ct, ang_ct = pos_ct()
        
        Sw = Sw()
        bw, AR = AR() 
        
        return pos_ba_cr, pos_bf_cr, ang_cr, pos_ba_ct, pos_bf_ct, ang_ct, Sw, bw, AR

    def ev(self):
        """
        Geometria do estabilizador vertical
        """
        corda = self.c_ev        
        
        def pos_cr():
            """
            Retorna a posiçao da corda raiz
            """
            y = self.bw_asat/6.0
            z_ba = z_bf = self.pos_ba_cr_asat[2]
            x_ba = self.pos_ba_cr_asat[0]
            x_bf = x_ba+self.c_asat
            return np.array([x_ba,y,z_ba]), np.array([x_bf,y,z_bf])
            
        def pos_ct():
            """
            Retorna a posiçao da corda da ponta
            """
            y = self.bw_asat/6.0
            x_ba = self.pos_ba_cr_asat[0]
            x_bf = x_ba+corda            
            z_ba = z_bf =  min(self.cone2(x_ba,y),self.cone2(x_bf,y))
                
            return np.array([x_ba,y,z_ba]), np.array([x_bf,y,z_bf])
            
        def Sw():
            S1 = np.cross(pos_bf_cr-pos_ba_cr,pos_ba_ct-pos_ba_cr)
            S2 = np.cross(pos_ba_ct-pos_bf_ct,pos_bf_cr-pos_bf_ct)
            return 0.5*np.linalg.norm(S1+S2)
            
        def AR():
            bw = pos_ba_ct[2]-pos_ba_cr[2]
            return bw, bw**2/Sw        
        
        pos_ba_cr, pos_bf_cr = pos_cr()
        pos_ba_ct, pos_bf_ct = pos_ct()
        
        Sw = Sw()
        bw, AR = AR()
    
        return pos_ba_cr, pos_bf_cr, pos_ba_ct, pos_bf_ct, Sw, bw, AR

    def formato(self):
        """
        Função que retorna os parâmetros das asas da forma necessária para
        as funções cg() e inercia() 
        """
        def espelhar(asa):
            asa_esp = np.array(asa)
            for i in asa_esp:
                i[1] *= -1
            return asa_esp      
            
        # Densidades superficiais
        ro_asaf = self.densidade_asas(self.AR_asaf)
        ro_asat = self.densidade_asas(self.AR_asat)        
        
        # Definição das asas
        asaf = np.array([self.pos_ba_cr_asaf,
                         self.pos_ba_ct_asaf,
                         self.pos_bf_ct_asaf,
                         self.pos_bf_cr_asaf,
                         self.pos_ba_cr_asaf])
        
        asat = np.array([self.pos_ba_cr_asat,
                         self.pos_ba_ct_asat,
                         self.pos_bf_ct_asat,
                         self.pos_bf_cr_asat,
                         self.pos_ba_cr_asat])
                         
        ev = np.array([self.pos_ba_cr_ev,
                       self.pos_ba_ct_ev,
                       self.pos_bf_ct_ev,
                       self.pos_bf_cr_ev,
                       self.pos_ba_cr_ev])
            
        asaf_esp = espelhar(asaf)
        asat_esp = espelhar(asat)
        ev_esp = espelhar(ev)
        
        # Array de triangulação
        tri = np.array([[0,3,1],[1,3,2]]) #Calculado pela função tri() em apoio.py
                         
        tri_asaf = []
        tri_asat = []
        tri_ev = []
        triangulo = lambda asa: [asa[tri[i][0]],asa[tri[i][1]],asa[tri[i][2]]]         
        for i in range(len(tri)):
                tri_asaf.append(triangulo(asaf))
                tri_asaf.append(triangulo(asaf_esp))
                tri_asat.append(triangulo(asat))
                tri_asat.append(triangulo(asat_esp))
                tri_ev.append(triangulo(ev))
                tri_ev.append(triangulo(ev_esp))
                         
        return ro_asaf, ro_asat, np.array(tri_asaf), np.array(tri_asat), np.array(tri_ev)
        
    def cg(self):
        """
        Função que retorna a posição estimada do CG do avião.
        """
        ro_asaf, ro_asat, asaf, asat, ev = self.formato()
        ro_ev = self.ro_estabilizador        
        
        m_asaf = ro_asaf*self.S_asaf
        m_asat = ro_asat*self.S_asat
        m_ev =  2.0*ro_ev*self.S_ev
        m_asas = m_asaf+m_asat+m_ev
        
        def bar(asa):
            n = len(asa)
            bar = 0
            for i in asa:
                p1,p2,p3 = i[0],i[1],i[2]
                bar += apoio.baricentro(p1,p2,p3)
            return bar/float(n)
        
        cg = (m_asaf*bar(asaf)+m_asat*bar(asat)+m_ev*bar(ev))/(m_asaf+m_asat+m_ev)
        cg = (m_asas*cg+self.m_motor*self.pos_motor)/(m_asas+self.m_motor)  
        p_vazio = m_asas+self.m_motor
        return cg, p_vazio
    
    def inercia(self, carga = False):
        """
        Calcula o tensor de inercia da aeronave
        """
        ro_asaf, ro_asat, asaf, asat, ev = self.formato() 
        ro_ev = self.ro_estabilizador        
        
        def J(asa,ro):
            """
            Calcula o tensor de inercia da asa em relacao ao CG
            """
            J = np.zeros((3,3)) #Tensor de inercia  
            for i in asa:
                p1,p2,p3 = i[0],i[1],i[2]
                J += apoio.inercia_tri(p1,p2,p3,self.pos_cg,ro)
            return J
            
        J_asaf = J(asaf,ro_asaf)
        J_asat = J(asat,ro_asat)
        J_ev = J(ev,ro_ev)
        J_motor = apoio.inercia_ponto(self.pos_motor,self.pos_cg,self.m_motor)
            
        if carga:
            # A carga é modelada como um cubo cujo baricentro está 
            # localizado no CG da aeronave
            L = (self.m_carga/self.ro_carga)**(1./3)
            J_carga = apoio.inercia_cubo(self.m_carga,L)
            return J_asaf+J_asat+J_motor+J_carga
        else:
            return J_asaf+J_asat+J_ev+J_motor
    
    def input_val(self,Nchord,Cspace,Nspace,Sspace):
        """
        Funçao que retorna os valores necessários para criar o arquivo 
        de input do AVL.
        """
        # Definiçao do header
        mach = self.vel/self.vel_som
        iYsym = iZsym = 0
        Zsym = 0.0
        Sref = self.S_asaf
        Cref = self.c_asaf
        Bref = self.bw_asaf
        Xref,Yref,Zref = self.pos_cg
        header = [mach,iYsym,iZsym,Zsym,Sref,Cref,Bref,Xref,Yref,Zref]
           
        # Definicao das asas
        DISC = [Nchord ,Cspace, Nspace, Sspace]
        YDUPLICATE = ANGLE = 0.0
            
        AFILE_raiz_asaf = "\"Perfis\%s.dat\"" %(self.perfilr_asaf)
        AFILE_ponta_asaf = "\"Perfis\%s.dat\"" %(self.perfilp_asaf)
        AFILE_raiz_asat = "\"Perfis\%s.dat\"" %(self.perfilr_asat)
        AFILE_ponta_asat = "\"Perfis\%s.dat\"" %(self.perfilp_asat)
        AFILE_ev = "\"Perfis\%s.dat\"" %(self.perfil_ev)
        
        def CLAF_asaf(string):
            if string == 'raiz':
                aerofolio = self.perfilr_asaf 
            elif string == 'ponta':
                aerofolio = self.perfilp_asaf 
            corda = self.c_asaf
            Re = self.ro_ar*corda*self.vel/self.mi_ar
            dcl = perfil.Analise(aerofolio).ajustelinear(Re,mach,self.p)['dCl'][0]
            return dcl*180/(2*np.pi**2)
            
        def CLAF_asat(string):
            if string == 'raiz':
                aerofolio = self.perfilr_asat 
            elif string == 'ponta':
                aerofolio = self.perfilp_asat 
            corda = self.c_asat
            Re = self.ro_ar*corda*self.vel/self.mi_ar
            dcl = perfil.Analise(aerofolio).ajustelinear(Re,mach,self.p)['dCl'][0]
            return dcl*180/(2*np.pi**2)
        
        if self.perfilr_asaf == 'x':
            CLAF_raiz_asaf = 1.0
        else:
            CLAF_raiz_asaf = CLAF_asaf('raiz')
        
        if self.perfilr_asat == 'x':
            CLAF_raiz_asat = 1.0
        else:
            CLAF_raiz_asat = CLAF_asat('raiz')
            
        if self.perfilp_asaf == 'x':
            CLAF_ponta_asaf = 1.0
        else:
            CLAF_ponta_asaf = CLAF_asaf('ponta')
            
        if self.perfilp_asat == 'x':
            CLAF_ponta_asat = 1.0
        else:
            CLAF_ponta_asat = CLAF_asat('ponta')
            
        CLAF_ev = 1.0

        def COORD(string):
            if string == 'raiz_asaf':
                Xle,Yle,Zle = self.pos_ba_cr_asaf
                Chord = self.c_asaf
                Ainc = self.ang_cr_asaf
            elif string == 'raiz_asat':
                Xle,Yle,Zle = self.pos_ba_cr_asat
                Chord = self.c_asat
                Ainc = self.ang_cr_asat
            elif string == 'ponta_asaf':
                Xle,Yle,Zle = self.pos_ba_ct_asaf
                Chord = self.c_asaf
                Ainc = self.ang_ct_asaf
            elif string == 'ponta_asat':
                Xle,Yle,Zle = self.pos_ba_ct_asat
                Chord = self.c_asat
                Ainc = self.ang_ct_asat
            elif string == 'raiz_ev':
                Xle,Yle,Zle = self.pos_ba_cr_ev
                Chord = self.c_asat
                Ainc = 0.0
            elif string == 'ponta_ev':
                Xle,Yle,Zle = self.pos_ba_ct_ev
                Chord = self.c_ev
                Ainc = 0.0
            return [Xle,Yle,Zle,Chord,Ainc,0,0]
            
        COORD_raiz_asaf = COORD('raiz_asaf')
        COORD_raiz_asat = COORD('raiz_asat')
        COORD_raiz_ev = COORD('raiz_ev')
        COORD_ponta_asaf = COORD('ponta_asaf')
        COORD_ponta_asat = COORD('ponta_asat')
        COORD_ponta_ev = COORD('ponta_ev')
        
        SECTION_frontal = [[COORD_raiz_asaf,AFILE_raiz_asaf,CLAF_raiz_asaf], 
                           [COORD_ponta_asaf,AFILE_ponta_asaf,CLAF_ponta_asaf]]        
        asa_frontal = [DISC, YDUPLICATE, ANGLE, SECTION_frontal]
    
        SECTION_traseira = [[COORD_raiz_asat,AFILE_raiz_asat,CLAF_raiz_asat], 
                            [COORD_ponta_asat,AFILE_ponta_asat,CLAF_ponta_asat]]  
        asa_traseira = [DISC, YDUPLICATE, ANGLE, SECTION_traseira]
        
        SECTION_ev = [[COORD_raiz_ev,AFILE_ev,CLAF_ev], 
                      [COORD_ponta_ev,AFILE_ev,CLAF_ev]]     
        est_vert = [DISC, YDUPLICATE, ANGLE, SECTION_ev]
            
        asas = {'Asa Frontal':asa_frontal,'Asa Traseira':asa_traseira, 'EV':est_vert}        
        
        return header, asas
