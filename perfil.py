# -*- coding: utf-8 -*-
"""
Created on Fri Oct 03 22:45:49 2014

@author: Diego Chou
"""

print("\nCarregando modulos de 'Perfil'...")

try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import InterpolatedUnivariateSpline
    print("Modulos de 'Perfil' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Perfil'\n")
    raise
    
import apoio

# Lista de variáveis:

# Nome do arquivo                               name                INPUT
# Número de Reynolds                            Re                  INPUT
# Lista os arquivos na pasta 'Perfis'           lerDir              FUNÇÃO
# Abre o arquivo para leitura                   openFile            FUNÇÃO
# Cl Cd e Cm por nº de Reynolds                 getClCdCm           OUTPUT 
# Cl Cd e Cm por nº de Reynolds e angulo        getClCdCm_alfa      OUTPUT                 
# CL Máximo                                     getClmax            OUTPUT
# CD mínimo                                     getCdmin            OUTPUT  

class Read(object): #Lê os arquivos na pasta Perfis 
    
    def __init__(self, name):
        self.name = name
        self.listaRe, self.listaMach = apoio.lerDirPerfil_Re_Mach(self.name,local=__file__)
    
    # Todos os valores necessários para a análise dos perfis
    # devem ser conseguidos pela classe Analise(Read) 
    
    @apoio.executarNaPasta('Perfis')
    def getPoints(self):
        """
        Função que retorna os pontos do arquivo do perfil.
        O arquivo deve estar no formato aceito pelo XFLR5
        """
        name = self.name
        filename = '%s.dat' %(name)
        
        try:
            f = open(filename, 'r')
        except IOError as er:
            print(er)
            print('')
            raise
            
        flines = f.readlines()[1:]
        listaX = []
        listaY = []
        
        for i in range(len(flines)):
#            print flines[i]
            words = flines[i].split()
#            print words[1]
            try:
                x = float(words[0])
                y = float(words[1])
                listaX.append(x)
                listaY.append(y)
            except:
                pass    
            
        #normalização do perfil
        corda = max(listaX)
        listaX = np.array(listaX)/corda
        listaY = np.array(listaY)/corda
        
        return list(listaX),list(listaY)
    
    @apoio.executarNaPasta('Perfis')
    def openFile(self, Re, mach):
        name = self.name
        filename = '%s_T1_Re%.3f_M%.2f_N9.0.txt' %(name, Re*1e-6, mach)
        try:
            f = open(filename, 'r')
        except IOError as er:
            print(er)
            print('')
            raise
            
        flines = f.readlines()
        return flines
    
    def getClCdCm(self, Re, mach):
        flines = self.openFile(Re, mach)
#        dicionario = {}
        lista = []
        
        for i in range(11,len(flines)):
#            print flines[i]
            words = flines[i].split()
#            print words[1]
            try:
                Cl = float(words[1])
                Cd = float(words[2])
                Cm = float(words[4])
                angulo = float(words[0])
#                dicionario.update({alpha:[Cl, Cd, Cm]})
                lista.append([angulo, [Cl, Cd, Cm]])
            except:
                pass
          
        return lista#, dicionario
        
    def getClCdCm_alfa(self, alfa, Re, mach):       
        lista = self.getClCdCm(Re, mach)
        
        for i in range(len(lista)):
            Cd = lista[i][1][1]
            Cl = lista[i][1][0]
            Cm = lista[i][1][2]
            angulo = lista[i][0]
            if angulo == alfa:
                lista = [angulo,[Cl, Cd, Cm]]
                break
            
        return lista
        
    def getXcp(self, Re, mach):  
        flines = self.openFile(Re, mach)
        lista = []
        
        for i in range(11,len(flines)):
            words = flines[i].split()
            try:
                xCp = float(words[9])
                angulo = float(words[0])
                lista.append([angulo, xCp])
            except:
                pass
            
        return lista
        
    def getXcp_alfa(self, alfa, Re, mach): 
        lista = self.getXcp(Re, mach)
        
        for i in range(len(lista)):
            xCp = lista[i][1]
            angulo = lista[i][0]
            if angulo == alfa:
                lista = [angulo,xCp]
                break
        
        return lista
        
                    
    def getClmax(self, Re, mach):
        flines = self.openFile(Re, mach)
        
        CLmax = 0
        angulo_clmax = 0
        for i in range(11,len(flines)):
#            print flines[i]
            words = flines[i].split()

            try:
                CL = float(words[1])
#                print CL
                angulo = float(words[0])
                if(CL>CLmax):
                    CLmax = CL
                    angulo_clmax = angulo 
            except:
                pass
                    
        return CLmax, angulo_clmax
    
    def getCdmin(self, Re, mach):
        flines = self.openFile(Re, mach)
        
        CDmin = 10
        angulo_cdmin = 0
        for i in range(11,len(flines)):
#            print flines[i]
            words = flines[i].split()
            
            try:
                CD = float(words[2])
#                print CD
                angulo = float(words[0])
                if(CD<CDmin):
                    CDmin = CD
                    angulo_cdmin = angulo 
            except:
                pass
                    
        return CDmin, angulo_cdmin
        

"""        
 Lista de variáveis:

 Principais valores de perfil para um Re           valores_perfil      OUTPUT
 Tabela com os principais valores do perfil        tab_perfil          OUTPUT  
 ajuste linear dos valores de Cl e Cm              ajustelinear        OUTPUT
 Interpolaçao para todos os valores de Re e alfa   interpReAlfa_Cd     OUTPUT
 Centro aerodinâmico do perfil (% da corda)        ac                  OUTPUT  
 Gráfico de Cl e Cd por angulo de ataque           plotClCd            OUTPUT  
"""

# Aqui devem ser definidas as analises dos perfis 
class Analise(Read): 
    
    def __init__(self, name, p = False):
        super(Analise,self).__init__(name)
        self.p = p
        
    def area(self, corda):
        """
        Calcula a área do perfil em função da corda.
        """
        
        def triangulacao_perfil(): 
            """
            Essa função triangula o perfil para que a área possa ser calculada
            pelo produto vetorial dos vetores dos triângulos gerados.
            """
            tri = np.array([[0.0, 0.0, 0.0]])
            N = N_total//5
            Ni = N//2 
            
            if Ni == N/2.0: # Caso com número par de pontos
                val = 1
            else:
                val = 0
           
            # Triangulacao do 1º perfil
            # Vale a regra da mão direita
            for i in range(0,Ni-1):
                tri = np.append(tri,[[i, i+1, N-i-2]], axis=0)
            
            if y[0] != y[-1]:
                tri = np.append(tri,[[0, N-2, N-1]], axis=0)
                
            for i in range(1,Ni-val):
                tri = np.append(tri,[[i, N-i-2, N-i-1]], axis=0)
            return tri[1:]        
        
        x,y = self.getPoints()
        x = np.array(x)
        y = np.array(y)
        N_total = len(x)    
        area = 0.0        

        tri = triangulacao_perfil()
        
        for i in range(len(tri)):
            # Calculando o vetor normal
            x0 = x[tri[i][0]]
            x1 = x[tri[i][1]]
            x2 = x[tri[i][2]]
            
            y0 = y[tri[i][0]]
            y1 = y[tri[i][1]]
            y2 = y[tri[i][2]]   
            
            AB = np.array([x1-x0, y1-y0,0.0])
            AC = np.array([x2-x0, y2-y0,0.0])
            p_vetorial = np.cross(AB, AC)
            n = np.linalg.norm(p_vetorial)
            area += n
        
        return float(area)*corda**2
        
    def spline(self):
        """
        Gera 2 splines a partir dos pontos do perfil. 
        Um para o extradorso e outro para o intradorso.
        Se o bordo de fuga tiver espessura, ela é zerada.
        """        
        name = self.name
        x,y = self.getPoints()
        
        ba = min(x) # abscissa do bordo de ataque
        index = x.index(ba) # posição do ponto na lista

        x1 = x[:index]
        y1 = y[:index]
        
        x2 = x[index:]
        y2 = y[index:]
        
        # Resolve o problema de x1 ou x2 poder estar em ordem decrescente
        if x1[0]>x1[-1]:
            x1 = x1[::-1] #Inverte a lista
            y1 = y1[::-1]
            
        if x2[0]>x2[-1]:
            x2 = x2[::-1]
            y2 = y2[::-1]
            
        if y1[1]>y2[-2]:
            if abs(y1[-1]) != 0.0:
                y1[-1] = 0.0
                y2[-1] = 0.0
            
            s_extra = InterpolatedUnivariateSpline(x1,y1)
            s_intra = InterpolatedUnivariateSpline(x2,y2)
        else:
            if abs(y2[-1]) != 0.0:
                y1[-1] = 0.0
                y2[-1] = 0.0    
            
            s_extra = InterpolatedUnivariateSpline(x2,y2)
            s_intra = InterpolatedUnivariateSpline(x1,y1)
        
        @apoio.executarNaPasta('Graficos/Aerodinamica')
        def grafico():
            beta = np.linspace(0,np.pi,200)
            xs = 0.5*(1.0-np.cos(beta))
            spline1 = s_extra(xs)
            spline2 = s_intra(xs) 
            
            plt.figure()
            plt.grid('on')
            plt.axis('equal')
            plt.title('%s' %(name))
            plt.xlim([-0.1,1.1])
            plt.plot(x,y,'.r',label = 'Pontos')
            plt.plot(xs,spline1,'-r')
            plt.plot(xs,spline2,'-b',label = 'Spline')
            plt.legend(loc='best')
            plt.savefig('perfil %s.png' %(name), bbox_inches='tight', dpi=200)
            
        if self.p:
            grafico()
            
        return s_extra,s_intra
        
    def linha_media(self):
        """
        Gera a função da linha média do perfil.
        """
        extra, intra = self.spline()
        beta = np.linspace(0,np.pi,200)
        xs = 0.5*(1.0-np.cos(beta))
        media = 0.5*(extra(xs)+intra(xs))        
        s_media = InterpolatedUnivariateSpline(xs,media)
        
        
        @apoio.executarNaPasta('Graficos/Aerodinamica')
        def grafico():
            name = self.name
            plt.figure()
            plt.grid('on')
            plt.axis('equal')
            plt.title('%s' %(name))
            plt.xlim([-0.1,1.1])
            plt.plot(xs,extra(xs),'-r')
            plt.plot(xs,intra(xs),'-r',label = 'Spline')
            plt.plot(xs,media,'--k',label = 'Linha Media')
            plt.legend(loc='best')
            plt.savefig('linha_media_%s.png' %(name), bbox_inches='tight', dpi=200)
            
        if self.p:
            grafico()        
        
        return s_media
        
    def coords_padrao(self,n):
        """
        A partir dos pontos lidos no arquivo .dat do perfil, gera um conjunto
        de coordenadas padrao, com pontos em posicoes especificas dos intra e
        extradorsos.
        
        Input:
        -----
        * n : numero de pontos em cada dorso do perfil (intra e extra)
        """
        beta = np.linspace(0,np.pi,n)
        s1,s2 = self.spline()
        xs = 0.5*(1.0-np.cos(beta))
        
        extradorso = s1(xs[::-1])
        intradorso = s2(xs[1:])
        x = np.append(xs[::-1],xs[1:])
        y = np.append(extradorso,intradorso)
        
        return list(x),list(y)
        
    def espessura_relativa(self, xp):
        """
        Calculo de e/c (razão de espessura) do perfil para uma posição relativa xp (% da corda).
        """
        x = self.getPoints()[0]
        ba = min(x) # abscissa do bordo de ataque
        bf = max(x) # abscissa do bordo de fuga
        c = bf-ba
        s1,s2 = self.spline()
        
        # xp é relativo (% da corda) e xa é absoluto
        xa = xp*c 
        
        # Calcula a espessura absoluta no ponto xs       
        e = abs(s1(xa)-s2(xa))
        
        return e/c # Retorna a espessura relativa à corda no ponto
        
    def t_c(self):
        """
        Calculo de t/c (thickness ratio) do perfil a partir dos pontos do arquivo.
        """
        x = self.getPoints()[0]
        ba = min(x) # abscissa do bordo de ataque
        bf = max(x) # abscissa do bordo de fuga
        c = bf-ba
        s1,s2 = self.spline()
        
        # Gera pontos que melhor contornam o perfil,
        # com maior densidade de pontos próximo ao bordo de ataque e de fuga 
        # onde as curvaturas são geralmente maiores
        beta = np.linspace(0,np.pi,500)
        xs = c*(0.5*(1.0-np.cos(beta))) 
        
        # Calcula os valores        
        spline1 = s1(xs)
        spline2 = s2(xs)
        
        val = abs(spline1-spline2)
        
        t = max(val)
        
        i = int(np.where(val==t)[0])
        xt_c = (xs[i]-ba)/c
        return t/c, xt_c
        
    def delta_y(self):
        """
        Coeficiente de afilamento do bordo de ataque do perfil em % da corda.
        Importante para o calculo do Ks e delta_alfa_s da asa.
        """
        x = self.getPoints()[0]
        ba = min(x) # abscissa do bordo de ataque
        bf = max(x) # abscissa do bordo de fuga
        c = bf-ba
        s1,s2 = self.spline()
        
        def y(x):
            y = max([s1(x),s2(x)]) # s1 e s2 podem estar trocados dependendo do perfil
            return y
        
        # Definição do Raymer
        x1 = 0.0015*c
        y1 = y(x1)
        
        x2 = 0.06*c
        y2 = y(x2)
         
        delta_y = 100*(y2-y1)/c # Em % da corda
        return delta_y
        
    def plotClCd(self, Re, mach):
        name = self.name
        pontos = self.getClCdCm(Re, mach)
        x = []
        y1 = []
        y2 = []              
        
        for i in range(len(pontos)):
            x.append(pontos[i][0])
            y1.append(pontos[i][1][0])
            y2.append(pontos[i][1][1])
        
        @apoio.executarNaPasta('Graficos/Aerodinamica')
        def grafico():
            # Gera o gráfico do Cl por angulo de ataque
            fig, ax1 = plt.subplots()
            plt.title(r"%s : $Re = %g$, $Mach = %.2f$" %(name, Re, mach))
            ax1.grid('on')
            ax1.plot(x, y1, 'b-')
            ax1.set_xlabel(r'$\alpha$ $(graus)$')
            ax1.set_ylabel(r'$C_l$', color='b')
            
            # Gera o gráfico do Cd por angulo de ataque
            ax2 = plt.twinx()
            ax2.plot(x, y2, 'r-')
            ax2.set_ylabel('$C_d$', color='r')
            plt.savefig("ClCdPorAlpha%s.png" %(name), bbox_inches='tight', dpi=200)
            
        
        grafico()
        
    def funcaoClCdCm(self, Re, mach):
        """
        Função que retorna os coeficientes do perfil para qualquer
        angulo dentro do intervalo de pontos,
        """
        if Re not in self.listaRe:
            raise ValueError("Re = %d nao listado para %s" %(Re,self.name))
        if mach not in self.listaMach:
            raise ValueError("Mach = %d nao listado para %s" %(mach,self.name))
            
        lista = self.getClCdCm(Re,mach)
        alfa = []
        Cl = []
        Cd = []
        Cm = []
        
        for i in range(len(lista)):
            alfa.append(lista[i][0])
            Cl.append(lista[i][1][0])
            Cd.append(lista[i][1][1])
            Cm.append(lista[i][1][2])
            
        sCl = InterpolatedUnivariateSpline(alfa,Cl,k=2)
        sCd = InterpolatedUnivariateSpline(alfa,Cd,k=2)
        sCm = InterpolatedUnivariateSpline(alfa,Cm,k=2)
            
        return sCl, sCd, sCm
        
    def funcaoXcp(self,Re,mach):
        """
        Função que retorna a posição do Cp do perfil em porcentagem
        para qualquer angulo dentro do intervalo de pontos,
        """
        if Re not in self.listaRe:
            raise ValueError("Re = %d nao listado para %s" %(Re,self.name))
        if mach not in self.listaMach:
            raise ValueError("Mach = %d nao listado para %s" %(mach,self.name))
            
        lista = self.getXcp(Re,mach)
        alfa = []
        xCp = []

        for i in range(len(lista)):
            alfa.append(lista[i][0])
            xCp.append(lista[i][1])
            
        sXcp = InterpolatedUnivariateSpline(alfa,xCp,k=2)
            
        return sXcp

    def interpXcp(self,alfa,Re,mach):
        """
        Interpola os valores tabelados para todos os ângulos, Reynolds e Machs
        dentro do intervalo existente.
        """
        
        ReAntes, ReDepois = apoio.intervalo(self.listaRe, Re)
        machAntes, machDepois = apoio.intervalo(self.listaMach, mach)
        
        #Interpolação linear do número de Mach
        val_ReAntes_machAntes = self.funcaoXcp(ReAntes,machAntes)(alfa)
        val_ReAntes_machDepois = self.funcaoXcp(ReAntes,machDepois)(alfa)
        
        val_ReAntes = apoio.interplinear(machAntes,val_ReAntes_machAntes,
                                         machDepois,val_ReAntes_machDepois,mach)
                                       
        val_ReDepois_machAntes = self.funcaoXcp(ReDepois,machAntes)(alfa)
        val_ReDepois_machDepois = self.funcaoXcp(ReDepois,machDepois)(alfa)
        
        val_ReDepois = apoio.interplinear(machAntes,val_ReDepois_machAntes,
                                         machDepois,val_ReDepois_machDepois,mach)
        
        #Interpolação linear do número de Reynolds
        val = apoio.interplinear(ReAntes,val_ReAntes,
                                 ReDepois,val_ReDepois,Re)
                                 
        return val

    def interpReMachAlfa(self, alfa, Re, mach, string):
        """
        Interpola os valores tabelados para todos os ângulos, Reynolds e Machs
        dentro do intervalo existente.
        """
        # Garante que as variáveis de entrada são do tipo 'float' do Python
        alfa = float(alfa)
        Re = float(Re)
        mach = float(mach)  
        
        if string == 'Cl':
            var = 0      
        elif string == 'Cd':
            var = 1
        elif string == 'Cm':
            var = 2
        else:
            raise ValueError("string = 'Cl','Cd' ou 'Cm'")
        
        ReAntes, ReDepois = apoio.intervalo(self.listaRe, Re)
        machAntes, machDepois = apoio.intervalo(self.listaMach, mach)
        
        #Interpolação linear do número de Mach
        val_ReAntes_machAntes = self.funcaoClCdCm(ReAntes,machAntes)[var](alfa)
        val_ReAntes_machDepois = self.funcaoClCdCm(ReAntes,machDepois)[var](alfa)
        
        val_ReAntes = apoio.interplinear(machAntes,val_ReAntes_machAntes,
                                         machDepois,val_ReAntes_machDepois,mach)
                                       
        val_ReDepois_machAntes = self.funcaoClCdCm(ReDepois,machAntes)[var](alfa)
        val_ReDepois_machDepois = self.funcaoClCdCm(ReDepois,machDepois)[var](alfa)
        
        val_ReDepois = apoio.interplinear(machAntes,val_ReDepois_machAntes,
                                         machDepois,val_ReDepois_machDepois,mach)
        
        #Interpolação linear do número de Reynolds
        val = apoio.interplinear(ReAntes,val_ReAntes,
                                 ReDepois,val_ReDepois,Re)
                                 
        return val
        
    def interpSplineReMach(self,Re, mach, string):
        """
        Função que retorna o spline da polar definida pela string, para qualquer
        números de Reynolds e Mach dentro do intervalo existente
        """
        # Garante que as variáveis de entrada são do tipo 'float' do Python
        Re = float(Re)
        mach = float(mach)  
        
        if string == 'Cl':
            var = 0      
        elif string == 'Cd':
            var = 1
        elif string == 'Cm':
            var = 2
        else:
            raise ValueError("string = 'Cl','Cd' ou 'Cm'")
            
        if (Re in self.listaRe) and (mach in self.listaMach):
            sCoeff = self.funcaoClCdCm(Re,mach)[var]
        else:
            ReAntes, ReDepois = apoio.intervalo(self.listaRe, Re)
            machAntes, machDepois = apoio.intervalo(self.listaMach, mach)
            
            #Interpolação linear do número de Mach
            s_Ra_Ma = self.funcaoClCdCm(ReAntes,machAntes)[var]
            s_Ra_Md = self.funcaoClCdCm(ReAntes,machDepois)[var]
            inicio_Ra = [s_Ra_Ma.get_knots()[0],s_Ra_Md.get_knots()[0]]
            fim_Ra = [s_Ra_Ma.get_knots()[-1],s_Ra_Md.get_knots()[-1]]
        
            s_Rd_Ma = self.funcaoClCdCm(ReDepois,machAntes)[var]
            s_Rd_Md = self.funcaoClCdCm(ReDepois,machDepois)[var]
            inicio_Rd = [s_Ra_Ma.get_knots()[0],s_Rd_Md.get_knots()[0]]
            fim_Rd = [s_Ra_Ma.get_knots()[-1],s_Rd_Md.get_knots()[-1]]
            
            x_R = np.linspace(max(inicio_Ra+inicio_Rd),min(fim_Ra+fim_Rd))
            y_Ra = ((machDepois-mach)*s_Ra_Ma(x_R)+(mach-machAntes)*s_Ra_Md(x_R))/(machDepois-machAntes)
            y_Rd = ((machDepois-mach)*s_Rd_Ma(x_R)+(mach-machAntes)*s_Rd_Md(x_R))/(machDepois-machAntes)
            
            #lista com os valores do coeficiente definido pela string para os
            # números de Reynolds e Mach de entrada.
            y_R = ((ReDepois-Re)*y_Ra+(Re-ReAntes)*y_Rd)/(ReDepois-ReAntes)
            sCoeff = InterpolatedUnivariateSpline(x_R,y_R,k=1)
                                 
        return sCoeff
        
    def interpReMach_Clmax(self, Re, mach):
        # Garante que as variáveis de entrada são do tipo 'float' do Python
        Re = float(Re)
        mach = float(mach)
        
        Re_antes, Re_depois = apoio.intervalo(self.listaRe, Re)
        mach_antes, mach_depois = apoio.intervalo(self.listaMach, mach)
        
        C_Ra_Ma, alfa_Ra_Ma = self.getClmax(Re_antes, mach_antes)
        C_Rd_Ma, alfa_Rd_Ma = self.getClmax(Re_depois, mach_antes)
        C_R_Ma = apoio.interplinear(Re_antes,C_Ra_Ma,Re_depois,C_Rd_Ma,Re)
        alfa_R_Ma = apoio.interplinear(Re_antes,alfa_Ra_Ma,Re_depois,alfa_Rd_Ma,Re)
        
        C_Ra_Md, alfa_Ra_Md = self.getClmax(Re_antes, mach_depois)
        C_Rd_Md, alfa_Rd_Md = self.getClmax(Re_depois, mach_depois)
        C_R_Md = apoio.interplinear(Re_antes,C_Ra_Md,Re_depois,C_Rd_Md,Re)
        alfa_R_Md = apoio.interplinear(Re_antes,alfa_Ra_Md,Re_depois,alfa_Rd_Md,Re)
        
        C_R_M = apoio.interplinear(mach_antes,C_R_Ma,mach_depois,C_R_Md,mach)
        alfa_R_M = apoio.interplinear(mach_antes,alfa_R_Ma,mach_depois,alfa_R_Md,mach)
        return C_R_M, alfa_R_M
        
    def interpReMach_Cdmin(self, Re, mach):
        # Garante que as variáveis de entrada são do tipo 'float' do Python
        Re = float(Re)
        mach = float(mach)
        
        Re_antes, Re_depois = apoio.intervalo(self.listaRe, Re)
        mach_antes, mach_depois = apoio.intervalo(self.listaMach, mach)
        
        C_Ra_Ma, alfa_Ra_Ma = self.getCdmin(Re_antes, mach_antes)
        C_Rd_Ma, alfa_Rd_Ma = self.getCdmin(Re_depois, mach_antes)
        C_R_Ma = apoio.interplinear(Re_antes,C_Ra_Ma,Re_depois,C_Rd_Ma,Re)
        alfa_R_Ma = apoio.interplinear(Re_antes,alfa_Ra_Ma,Re_depois,alfa_Rd_Ma,Re)
        
        C_Ra_Md, alfa_Ra_Md = self.getCdmin(Re_antes, mach_depois)
        C_Rd_Md, alfa_Rd_Md = self.getCdmin(Re_depois, mach_depois)
        C_R_Md = apoio.interplinear(Re_antes,C_Ra_Md,Re_depois,C_Rd_Md,Re)
        alfa_R_Md = apoio.interplinear(Re_antes,alfa_Ra_Md,Re_depois,alfa_Rd_Md,Re)
        
        C_R_M = apoio.interplinear(mach_antes,C_R_Ma,mach_depois,C_R_Md,mach)
        alfa_R_M = apoio.interplinear(mach_antes,alfa_R_Ma,mach_depois,alfa_R_Md,mach)
        return C_R_M, alfa_R_M
    
    def ajustelinear(self, Re, mach, p=False):
        name = self.name
        
#        if Re not in self.listaRe:
#            Re = apoio.intervalo(self.listaRe,Re)[0]
#            
#        if mach not in self.listaMach:
#            mach = apoio.intervalo(self.listaMach,mach)[0]

#        lista = self.getClCdCm(Re,mach)
        
        alpha_CLmax = self.interpReMach_Clmax(Re,mach)[1]
        sCL = self.interpSplineReMach(Re,mach,'Cl')
        sCM = self.interpSplineReMach(Re,mach,'Cm')
        
        alphaCL = np.linspace(sCL.get_knots()[0],alpha_CLmax)
        CL = sCL(alphaCL)        
        alphaCM = np.linspace(0.0,alpha_CLmax-5.0)
        CM= sCM(alphaCM)
        
        a1, b1 = np.polyfit(alphaCL, CL, deg=1)
        a3, b3 = np.polyfit(alphaCM, CM, deg=1)
        
        @apoio.executarNaPasta('Graficos/Aerodinamica')
        def grafico():
            def sinal(b):
                if b >= 0.0:
                    sinal = '+'
                else:
                    sinal = '-'
                return sinal
            
            y1 = a1*alphaCL + b1
            y3 = a3*alphaCM + b3
            plt.figure()
            ax1 = plt.subplot2grid((2,1), (0,0), colspan=2)
            plt.title('Perfil %s; $Re = %g$, $Mach=%.2f$' %(name,Re,mach))
            ax1.grid('on')
            ax1.plot(alphaCL, y1, 'r-', label=r'$C_l = %f*\alpha %s %f$' %(a1,sinal(b1),abs(b1))) 
            ax1.plot(alphaCL, CL, 'o')
            ax1.set_ylabel('$C_l$')
            ax1.legend(loc='best')
            
            ax3 = plt.subplot2grid((2,1), (1,0), rowspan=1)
            ax3.grid('on')
            ax3.plot(alphaCM, y3, 'r-', label=r'$C_m = %f*\alpha %s %f$' %(a3,sinal(b3),abs(b3)))
            ax3.plot(alphaCM, CM, 'o')
            ax3.set_xlabel(r'$\alpha$ ($graus$)')
            ax3.set_ylabel('$C_m$')
            ax3.legend(loc='best')
            
            plt.savefig("ajuste%s_%g.png" %(name, Re), bbox_inches='tight', dpi=200)
        
        if self.p or p:
            grafico()
        
        return {'dCl':(a1,b1), 'dCm':(a3,b3)}    
        
    @staticmethod
    def ac(a0, m0):
        """
        Posição do centro aerodinamico do perfil relativa a corda. 
        Os valores de Cm para o perfil são tirados a 25% da corda no XFLR5
        
        m0 : dcm/dalfa
        a0 : dcl/dalfa
        """
        xac = -m0/a0  
        ac = 0.25+xac
        return ac
        
    def ac_medio(self, mach):
        """
        O centro aerodinamico do perfil varia muito pouco com o
            número de Reynolds.
        Essa função tira a media do centro aerodinamico para todos os 
            Reynolds disponíveis para o numero de Mach do perfil em questão
        """
        lista = self.listaRe
        ac = []
        
        for Re in lista:
            a = self.ajustelinear(Re, mach)
            ac.append(self.ac(a['dCl'][0],a['dCm'][0]))
    
        ac_medio = sum(ac)/len(ac)
        return ac_medio
        
        
    #Valores do perfil para um único número de Reynolds
    def valor_perfil(self, Re, mach):
        arquivo = self.name
        print("\n-> %s (Re = %g, Mach = %.3f)" %(arquivo, Re, mach))
        
        dicionario = {'Alfa_clmax':0.0, 
                      'Clmax':0.0,
                      'dcl/dalfa':0.0,
                      'Alfa_cdmin':0.0, 
                      'Cdmin':0.0,
                      'dcm/dalfa':0.0,
                      'AC':0.0} 
                      
        p = self.interpReMach_Clmax(Re,mach)
        dicionario['Alfa_clmax'] = p[1]
        dicionario['Clmax'] = p[0]
        
        d = self.interpReMach_Cdmin(Re,mach)
        dicionario['Alfa_cdmin'] = d[1]
        dicionario['Cdmin'] = d[0]
        
        a = self.ajustelinear(Re,mach,False)
        dicionario['dcl/dalfa'] = a['dCl'][0]
        dicionario['dcm/dalfa'] = a['dCm'][0]
        
        c = self.ac(a['dCl'][0],a['dCm'][0])
        dicionario['AC'] = c
        
        for keys,values in dicionario.items():
            print('\t%s : %f' %(keys,values))
            
        return dicionario
    
    # Valores do perfis para vários números de Reynolds
    def tab_perfil(self, mach):
        lista = self.listaRe
        arquivo = self.name        
        print("\n-> %s, Mach = %.3f" %(arquivo,mach))
        dicionario = {'1. Alfa_clmax':[], 
                      '2. Clmax':[],
                      '3. Alfa_cdmin':[], 
                      '4. Cdmin':[],
                      '5. dcl/dalfa':[],                 
                      '6. dcm/dalfa':[],
                      '7. XAC':[],
                      'Re':[]} 
        
        for Re in lista:
            dicionario['Re'].append(Re)
            p = self.getClmax(Re, mach)
            dicionario['1. Alfa_clmax'].append(p[1])
            dicionario['2. Clmax'].append(p[0])
            d = self.getCdmin(Re, mach)
            dicionario['3. Alfa_cdmin'].append(d[1])
            dicionario['4. Cdmin'].append(d[0])
            a = self.ajustelinear(Re,mach)
            dicionario['5. dcl/dalfa'].append(a['dCl'][0])
            dicionario['6. dcm/dalfa'].append(a['dCm'][0])
            c = self.ac(a['dCl'][0],a['dCm'][0])
            dicionario['7. XAC'].append(c)
            
        df = pd.DataFrame(dicionario).set_index('Re')
        return df
        
    def K2_Nicolai(self, Re, mach):
        """
        Este método refere-se ao cálculo do coeficiente K'' 
        de arrasto induzido viscoso citado no paper Nicolai
        """
        name = self.name
        angulos = np.array(range(0,11))
        Cl = self.interpSplineReMach(Re,mach,'Cl')(angulos)
        Cd = self.interpSplineReMach(Re,mach,'Cd')(angulos)
        
        def Clmin():
            alfa = self.interpReMach_Cdmin(Re, mach)[1]
            Clmin = self.interpReMachAlfa(alfa,Re,mach,'Cl')
            return Clmin
            
        x = (Cl-Clmin())**2
        y = Cd
        
        X = np.vstack([x, np.ones(len(x))])
        a1, b1 = np.linalg.lstsq(X.T, y)[0]
        
        @apoio.executarNaPasta('Graficos/Aerodinamica')
        def grafico():
            yr = a1*x + b1
            plt.figure()
            
            plt.title(r"Fator K'' para perfil %s; $Re = %g$, $Mach = %.2f$" %(name,Re,mach))
            plt.grid('on')
            plt.plot(x, yr, 'r-', label=r"$K'' = \frac{\Delta y}{\Delta x} = %f$" %(a1)) 
            plt.plot(x, y, 'o')
            plt.xlabel('$(C_l-C_{l_{min}})^2$')
            plt.ylabel('$C_d$')
            plt.legend(loc='best')
            
            plt.savefig("K''%s_%g.png" %(name, Re), bbox_inches='tight', dpi=200)            
            
        if self.p:
            grafico()            
            
        return a1
        
    def principal(self, mach):
        tabelaAsa = self.tab_perfil(mach)
        print(tabelaAsa)

if __name__ == "__main__":    
    perfil = 'S1223 MOD2015'   
    analise = Analise(perfil,p=True)
    
    analise.linha_media()