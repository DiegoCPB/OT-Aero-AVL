# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 21:05:25 2016

@author: Diego Chou
"""

print("\nCarregando modulos de 'Construtor'...")

try:
    import numpy as np
    import sympy as sym
    print("Modulos de 'Construtor' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Construtor'\n")
    raise

import apoio
import write_input_file as wif
import perfil

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
    
    INPUTS:
    ------
    * ang_clear: angulo de folga para que a ponta da asa não toque o chão
    * dz_asas: Altura entre as asas
    * asaf: Asa frontal
    * asat: Asa traseira
    * cr: Corda raiz
    * ct: Corda ponta
    * ang: Angulo de incidencia
    * enflex: Enflexamento do bordo de ataque
    * epsilon: Angulo de torcao
    """
    
    d_cone = 2.5 # Diâmetro do cone
    h_cone = 0.75 # Altura do cone
    x_frente = 0.330-0.5*d_cone # Ponto mais a frente da aeronave
    z_min = 0.04 # Altura do ponto mais baixo da aeronave 
    m_motor = 0.731 # Peso de motor e hélice
    pos_motor = np.array([x_frente+0.063, 0.0, 0.198])
    k_asa = [0.41716, 0.5] # Constantes para o calculo da densidade da asa 
    m_carga = 12 # Carga projetada em 2015
    ro_carga = 8.8e3 #Densidade do bronze
    vel_som = 340.0
    mi_ar = 1.962e-5
    ro_ar = 1.086

    #Parametros de discretizaçao das asas
    Nchord = 8
    Cspace = 1.0
    Nspan = 12
    Sspace = -2.0
    
    def __init__(self,name,ang_clear,dz_asas, vel, perfil_raiz, perfil_ponta,
                 cr_asaf,ct_asaf,ang_asaf,enflex_asaf,epsilon_asaf,
                 cr_asat,ang_asat):    
        # Parametros gerais iteráveis
        self.ang_clear = ang_clear
        self.dz_asas = dz_asas 
        self.vel = vel
        self.perfil_raiz = str(perfil_raiz)
        self.perfil_ponta = str(perfil_ponta)
        
        # Parametros de asa frontal        
        self.cr_asaf = cr_asaf
        self.ct_asaf = ct_asaf
        self.ang_asaf = ang_asaf
        self.enflex_asaf = enflex_asaf
        self.epsilon_asaf = epsilon_asaf
        
        #Parametros de asa traseira
        self.cr_asat = cr_asat
        self.ct_asat = ct_asaf 
        self.ang_asat = ang_asat
        self.epsilon_asat = ang_asaf+epsilon_asaf-ang_asat
        
        #Calculos para a asa frontal
        self.pos_ba_cr_asaf, self.pos_bf_cr_asaf, self.ang_cr_asaf,\
        self.pos_ba_ct_asaf, self.pos_bf_ct_asaf, self.ang_ct_asaf,\
        self.S_asaf, self.bw_asaf, self.AR_asaf = self.asa_frontal()
        
        #Calculos para a asa traseira 
        self.pos_ba_cr_asat, self.pos_bf_cr_asat, self.ang_cr_asat,\
        self.pos_ba_ct_asat, self.pos_bf_ct_asat, self.ang_ct_asat,\
        self.S_asat, self.bw_asat, self.AR_asat = self.asa_traseira()
        
        #Posição do CG levando em conta pesos do motor e asas
        self.pos_cg, self.p_vazio = self.cg()
        
        #Tensor de inercia do aviao
        self.J = self.inercia()
        
        # Escreve arquivo '.avl'
        header, asas = self.input_val(self.Nchord,self.Cspace,self.Nspan,self.Sspace)
        wif.Input_geometry(name,header,asas)
        
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
        else:
            R = c*(z0-z)
            b1 = np.sqrt(R**2-a**2)
            b2 = -b1
            return b1,b2
    
    def asa_frontal(self):
        """
        Geometria da asa frontal
        """
        cr = self.cr_asaf
        ct = self.ct_asaf
        iw = self.ang_asaf
        enflex = self.enflex_asaf
        epsilon = self.epsilon_asaf
        
        def pos_cr():
            """
            Retorna o ponto do bordo de ataque e a incidencia
            da corda raiz
            """
            y = 0.0
            x_ba = self.x_frente
            x_bf = x_ba+cr
            z_bf = self.z_min            
            z_ba = z_bf+cr*np.sin(iw*np.pi/180)
            
            ang_cr = iw
            
            return np.array([x_ba,y,z_ba]), np.array([x_bf,y,z_bf]), ang_cr
            
        def pos_ct():
            """
            Retorna o ponto do bordo de ataque e a incidencia
            da corda da ponta
            """
            z0 = self.h_cone
            c = self.d_cone/(2.*self.h_cone)
            z_ba = 2.5*np.tan(self.ang_clear*np.pi/180)
            ang_ct = iw+epsilon            
            
            def calc(string):
                """
                Calculo das posições (x,y)
                
                STRING:
                ------
                * 'ba': Enflexamento em relação ao bordo de ataque
                * 'bf': Enflexamento em relação ao bordo de fuga
                """
                z_bf = z_ba - np.tan(ang_ct*np.pi/180)*ct
                Y = sym.Symbol('Y')
                R = c*(z0-z_ba)
                if string == 'ba':
                    x = lambda y: y*np.tan(enflex*np.pi/180)+self.x_frente
                    Y = sym.solvers.solve(x(Y)**2+Y**2-R**2, Y)
                    y = float(max(Y))
                    x_ba = x(y)
                    x_bf = x_ba+ct
                    
                elif string == 'bf':
                    x = lambda y: y*np.tan(enflex*np.pi/180)+self.x_frente+ct
                    Y = sym.solvers.solve(x(Y)**2+Y**2-R**2, Y)           
                    y = float(max(Y))
                    x_bf = x(y)
                    x_ba = x_bf-ct
                
                return x_ba, x_bf, y, z_bf
                
            x_ba, x_bf, y , z_bf = calc('ba')
            
            # Se o valor de y para o ponto do bordo de fuga estiver fora do cone
            # a função calc é refeita em função do bordo de fuga.
            if y > max(self.cone(x_bf,z_ba)):
                x_ba, x_bf, y, z_bf = calc('bf')    
            
            return np.array([x_ba,y,z_ba]), np.array([x_bf,y,z_bf]), ang_ct
            
        def Sw():
            S1 = np.cross(pos_bf_cr-pos_ba_cr,pos_ba_ct-pos_ba_cr)
            S2 = np.cross(pos_ba_ct-pos_bf_ct,pos_bf_cr-pos_bf_ct)
            return np.linalg.norm(S1+S2)
            
        def AR():
            bw = 2*((pos_ba_ct[1]-pos_ba_cr[1])**2+(pos_ba_ct[2]-pos_ba_cr[2])**2)
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
        cr = self.cr_asat
        iw = self.ang_asat
        
        def pos_cr():
            """
            Retorna o ponto do bordo de ataque e a incidencia
            da corda raiz
            """
            r_cone = 0.5*self.d_cone
            h_cone = self.h_cone
            y = 0.0
            z_ba = self.pos_ba_cr_asaf[2]+self.dz_asas
            z_bf = z_ba - cr*np.sin(iw*np.pi/180)
            x_bf = -z_ba*r_cone/h_cone + r_cone
            x_ba = x_bf-cr
            ang_cr = iw
            
            return np.array([x_ba,y,z_ba]), np.array([x_bf,y,z_bf]), ang_cr
            
        def pos_ct():
            """
            Retorna o ponto do bordo de ataque e a incidencia
            da corda da ponta. Coincidente com a asa frontal.
            """
            return self.pos_ba_ct_asaf, self.pos_bf_ct_asaf, self.ang_ct_asaf
        
        def Sw():
            S1 = np.cross(pos_bf_cr-pos_ba_cr,pos_ba_ct-pos_ba_cr)
            S2 = np.cross(pos_ba_ct-pos_bf_ct,pos_bf_cr-pos_bf_ct)
            return np.linalg.norm(S1+S2)
            
        def AR():
            bw = 2*((pos_ba_ct[1]-pos_ba_cr[1])**2+(pos_ba_ct[2]-pos_ba_cr[2])**2)
            return bw, bw**2/Sw        
        
        pos_ba_cr, pos_bf_cr, ang_cr = pos_cr()
        pos_ba_ct, pos_bf_ct, ang_ct = pos_ct()
        
        Sw = Sw()
        bw, AR = AR() 
        
        return pos_ba_cr, pos_bf_cr, ang_cr, pos_ba_ct, pos_bf_ct, ang_ct, Sw, bw, AR

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
            
        asaf_esp = espelhar(asaf)
        asat_esp = espelhar(asat)
        
        # Array de triangulação
        tri = np.array([[0,3,1],[1,3,2]]) #Calculado pela função tri() em apoio.py
                         
        tri_asaf = []
        tri_asat = []
        triangulo = lambda asa: [asa[tri[i][0]],asa[tri[i][1]],asa[tri[i][2]]]         
        for i in range(len(tri)):
                tri_asaf.append(triangulo(asaf))
                tri_asaf.append(triangulo(asaf_esp))
                tri_asat.append(triangulo(asat))
                tri_asat.append(triangulo(asat_esp))
                         
        return ro_asaf, ro_asat, np.array(tri_asaf), np.array(tri_asat)
        
    def cg(self):
        """
        Função que retorna a posição estimada do CG do avião.
        """
        ro_asaf, ro_asat, asaf, asat = self.formato()
        
        m_asas = ro_asaf*self.S_asaf+ro_asat*self.S_asat
        
        def bar(asa):
            n = len(asa)
            bar = 0
            for i in asa:
                p1,p2,p3 = i[0],i[1],i[2]
                bar += apoio.baricentro(p1,p2,p3)
            return bar/float(n)
        
        cg = (ro_asaf*bar(asaf)+ro_asat*bar(asat))/(ro_asaf+ro_asat)
        cg = (m_asas*cg+self.m_motor*self.pos_motor)/(m_asas+self.m_motor)  
        p_vazio = m_asas+self.m_motor
        return cg, p_vazio
    
    def inercia(self, carga = False):
        """
        Calcula o tensor de inercia da aeronave
        """
        ro_asaf, ro_asat, asaf, asat = self.formato() 
        
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
        J_motor = apoio.inercia_ponto(self.pos_motor,self.pos_cg,self.m_motor)
            
        if carga:
            # A carga é modelada como um cubo cujo baricentro está 
            # localizado no CG da aeronave
            L = (self.m_carga/self.ro_carga)**(1./3)
            J_carga = apoio.inercia_cubo(self.m_carga,L)
            return J_asaf+J_asat+J_motor+J_carga
        else:
            return J_asaf+J_asat+J_motor
    
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
        Cref = self.cr_asaf
        Bref = self.bw_asaf
        Xref,Yref,Zref = self.pos_cg
        header = [mach,iYsym,iZsym,Zsym,Sref,Cref,Bref,Xref,Yref,Zref]
           
        # Definicao das asas
        DISC = [Nchord ,Cspace, Nspace, Sspace]
        YDUPLICATE = ANGLE = 0.0
        SCALE = TRANSLATE = [0.0,0.0,0.0]
            
        AFILE_raiz = "\"Perfis\%s.dat\"" %(self.perfil_raiz)
        AFILE_ponta = "\"Perfis\%s.dat\"" %(self.perfil_ponta)
        
        def CLAF_raiz(asa):
            if asa == 'asaf':
                corda = self.cr_asaf
            elif asa == 'asat':
                corda = self.cr_asat 
            Re = self.ro_ar*corda*self.vel/self.mi_ar
            dcl = perfil.Analise(self.perfil_raiz).ajustelinear(Re,mach)['dCl'][0]
            return dcl*180/(2*np.pi**2)
            
        def CLAF_ponta():
            Re = self.ro_ar*self.ct_asaf*self.vel/self.mi_ar
            dcl = perfil.Analise(self.perfil_ponta).ajustelinear(Re,mach)['dCl'][0]
            return dcl*180/(2*np.pi**2)
            
        CLAF_raiz_asaf = CLAF_raiz('asaf')
        CLAF_raiz_asat = CLAF_raiz('asat')
        CLAF_ponta = CLAF_ponta()
        
        def COORD(string):
            if string == 'raiz_asaf':
                Xle,Yle,Zle = self.pos_ba_cr_asaf
                Chord = self.cr_asaf
                Ainc = self.ang_cr_asaf
            elif string == 'raiz_asat':
                Xle,Yle,Zle = self.pos_ba_cr_asat
                Chord = self.cr_asat
                Ainc = self.ang_cr_asat
            elif string == 'ponta_asaf':
                Xle,Yle,Zle = self.pos_ba_ct_asaf
                Chord = self.ct_asaf
                Ainc = self.ang_ct_asaf
            elif string == 'ponta_asat':
                Xle,Yle,Zle = self.pos_ba_ct_asat
                Chord = self.ct_asat
                Ainc = self.ang_ct_asat
            return [Xle,Yle,Zle,Chord,Ainc,0,0]
            
        COORD_raiz_asaf = COORD('raiz_asaf')
        COORD_raiz_asat = COORD('raiz_asat')
        COORD_ponta_asaf = COORD('ponta_asaf')
        COORD_ponta_asat = COORD('ponta_asat')
        
        SECTION_frontal = [[COORD_raiz_asaf,AFILE_raiz,CLAF_raiz_asaf], 
                           [COORD_ponta_asaf,AFILE_ponta,CLAF_ponta]]        
        asa_frontal = [DISC, YDUPLICATE, SCALE, TRANSLATE, ANGLE, SECTION_frontal]
    
        SECTION_traseira = [[COORD_raiz_asat,AFILE_raiz,CLAF_raiz_asat], 
                            [COORD_ponta_asat,AFILE_ponta,CLAF_ponta]]  
        asa_traseira = [DISC, YDUPLICATE, SCALE, TRANSLATE, ANGLE, SECTION_traseira]
        
        asas = {'Asa_frontal':asa_frontal,'Asa_traseira':asa_traseira}        
        
        return header, asas
            
    def plot(self):
        """
        Função de debug pra verificar em um gráfico 3d se os pontos foram
        gerados nas posições corretas.
        """
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        
        asaf, asat = self.formato()[2:]        
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d([-1.25,1.25])
        ax.set_ylim3d([-1.25,1.25])
        ax.set_zlim3d([-1.25,1.25])
        
        cg = self.pos_cg
        ax.plot([cg[0]],[cg[1]],[cg[2]],'ro', label='CG')        
        
        for i in range(len(asaf)):
            M = asaf[i]
            x = np.append(M[:,0],M[0,0])
            y = np.append(M[:,1],M[0,1])
            z = np.append(M[:,2],M[0,2])
            ax.plot(x,y,z)   
                    
        for i in range(len(asat)):
            M = asat[i]
            x = np.append(M[:,0],M[0,0])
            y = np.append(M[:,1],M[0,1])
            z = np.append(M[:,2],M[0,2])
            ax.plot(x,y,z)
        
        plt.legend(loc='best')
        plt.show()
        
if __name__ == '__main__':
    name = 'Test plane'
    ang_clear = 3.0
    dz_asas = 0.35
    vel = 15
    perfil_raiz = 'S1223 MOD2015'
    perfil_ponta = 'MIN ponta2016'
    cr_asaf = 0.5
    ct_asaf = 0.25
    ang_asaf = 5.0
    enflex_asaf = 40.0
    epsilon_asaf = -3.0
    cr_asat = 0.3
    ang_asat = 2.0
    
    aviao = Construtor2016(name,ang_clear,dz_asas, vel, perfil_raiz, perfil_ponta,
                           cr_asaf,ct_asaf,ang_asaf,enflex_asaf,epsilon_asaf,
                           cr_asat,ang_asat)
    aviao.plot()
