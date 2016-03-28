# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 21:10:02 2014

@author: Diego Chou Pazo Blanco
"""

# --------------------------------- LICENCE  -------------------------------- #                 
#                                                                             #              
#    Copyrighted 2015 by Minerva Aerodesign UFRJ - aerodesign@poli.ufrj.br    #              
#                                                                             #            
#    This program is free software: you can redistribute it and/or modify     #            
#    it under the terms of the GNU General Public License as published by     #          
#    the Free Software Foundation, either version 3 of the License, or        #             
#    (at your option) any later version.                                      #             
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #                                     
#    but WITHOUT ANY WARRANTY without even the implied warranty of            #               
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
# --------------------------------------------------------------------------- #                                                          
#                                                                                                                                 
# É importante verificar se o perfil escolhido cumpre os requisitos de formato necessários.                          
# Os pontos devem ser definidos a partir do bordo de fuga e retornando a ele apos seguir o contorno do perfil.
                                                      
print("\nCarregando modulos de 'aeronave3D'...")

try:
    import numpy as np
    import matplotlib.pyplot as plt
    print("Modulos de 'perfil2asa' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'aeronave3D'\n")
    raise

from apoio import executarNaPasta
import perfil as airfoil

class Asa3D(object):
    """
    Essa classe cuida da triangulação da asa.
    """
    
    def __init__(self, nome_raiz, nome_ponta, alpha, washout, offset, 
                 corda_raiz, distancia1, corda_ponta, distancia2, corda_mac, espelhar):
        self.perfil_raiz = airfoil.Analise(nome_raiz)
        self.perfil_ponta = airfoil.Analise(nome_ponta)
        self.espelhar = espelhar
        self.a = alpha
        self.wo = washout
        self.off = offset
        self.c1 = corda_raiz
        self.d1 = distancia1
        self.c2 = corda_ponta
        self.d2 = distancia2
        self.c_mac = corda_mac
    
    def coords_perfis(self, n=50):
        """
        Retorna as coordenadas dos perfis da raiz e da ponta de asa.
        Se os perfis sao iguais, apenas retorna os pontos lido no .dat
        no formato apropriado. Senao, aplica splines no perfis para que
        ambos possuam coordenadas nas mesmas posicoes da corda.
        
        Input:
        -----
        * n : numero de pontos em cada dorso do perfil (intra e extra)
        
        Output:
        ------
        * coords : (coordenadas raiz, coordenadas ponta)
        """        
        if self.perfil_raiz == self.perfil_ponta:
            c = self.perfil_raiz.getPoints()
            coords = []            
            for i in range(len(c[0])):
                ponto = np.array([[c[0][i]],[c[1][i]]])
                coords.append(ponto)
            
            return coords,coords
            
        else:
            cr = self.perfil_raiz.coords_padrao(n)
            cp = self.perfil_ponta.coords_padrao(n)
            coords_r = []
            coords_p = []
            for i in range(len(cr[0])):
                ponto = np.array([[cr[0][i]],[cr[1][i]]])
                coords_r.append(ponto)
                ponto = np.array([[cp[0][i]],[cp[1][i]]])
                coords_p.append(ponto)
            
            return coords_r,coords_p
    
    @staticmethod        
    def escala(corda, perfil):
        dim = corda
        
        pos0 = perfil
    
        for i in range(len(pos0)):
            pos0[i] = dim*pos0[i] 
            
        return pos0
    
    @staticmethod    
    def translacao(x, y, perfil):
        trans2d = np.array([[x],[y]])
    
        pos0 = perfil
        pos1 = []
    
        for i in range(len(pos0)):
            a = trans2d + pos0[i] 
            pos1.append(a)
        
        return pos1
    
    @staticmethod
    def rotacao2d(alpha, perfil):
        angulo = alpha*np.pi/180
        giro2d = np.matrix([[np.cos(angulo), np.sin(angulo)],[-np.sin(angulo), np.cos(angulo)]])
        
        pos0 = perfil
    
        for i in range(len(pos0)):
            pos0[i] = giro2d * pos0[i]
            
        return pos0
    
    @staticmethod
    def grafico2d(perfil, corda, alpha, titulo): 
        x = np.zeros(len(perfil))
        y = np.zeros(len(perfil))
        
        for i in range(len(perfil)):
            x[i] = perfil[i][0]
            y[i] = perfil[i][1]
            
        plt.figure()
        plt.axis('equal')
        plt.grid('on')
        plt.title(titulo[0:-1] + ": %.1f graus" %(alpha))
        plt.plot(x, y, '.-r')
        plt.xlim([-0.25,0.25])#([min(x)-0.1*corda, max(x)+0.1*corda])
    
    def tornar_3D(self, perfil_central, perfil_trans, perfil_ponta):
        if len(perfil_central) != len(perfil_ponta):
            raise ValueError("Os perfis nao tem a mesma quantidade de pontos")

        dim = len(perfil_central) 
        
        # Perfil central         
        x0 = np.zeros(dim)
        z0 = np.zeros(dim)
        y0 = np.zeros(dim)
        
        # Perfil de transição
        x1 = np.zeros(dim)
        z1 = np.zeros(dim)
        y1 = self.d1*np.ones(dim)
        
        # Perfil da ponta
        x2 = np.zeros(dim)
        z2 = np.zeros(dim)
        y2 = self.d2*np.ones(dim)
    
        for i in range(dim):
            x0[i] = perfil_central[i][0]
            z0[i] = perfil_central[i][1] 
            
            x1[i] = perfil_trans[i][0]
            z1[i] = perfil_trans[i][1]
            
            x2[i] = perfil_ponta[i][0]
            z2[i] = perfil_ponta[i][1]
        
        # Cria uma lista com todas as posições dos 5 perfis ou 3 perfis
        if self.espelhar:
            x = np.append(np.append(np.append(np.append(x2,x1),x0),x1),x2)
            z = np.append(np.append(np.append(np.append(z2,z1),z0),z1),z2)
            y = np.append(np.append(np.append(np.append(-y2,-y1),y0),y1),y2)
        else: 
            x = np.append(np.append(x2,x1),x0)
            z = np.append(np.append(z2,z1),z0)
            y = np.append(np.append(-y2,-y1),y0)

        N = len(x)

        return x, y, z, N 
        
    # A triangulaçao dos perfis da ponta é específica para o caso em que estes 
    # são o 1º e o 5º na ordem dada pela funçao "tornar_3D()"    
    def triangulacao(self,z, N_total, tri): 
        if self.espelhar:
            # São 5 perfis no total
            N = N_total//5
            # São 4 seçoes de casca
            lista = np.linspace(0,3*N,4)
        else:
            # São 3 perfis no total
            N = N_total//3
            # São 2 seçoes de casca
            lista = np.linspace(0,N,2)
            
        Ni = N//2 
        
        if Ni == N/2.0: # Caso com número par de pontos
            val = 1
        else:
            val = 0
       
        # Triangulacao do 1º perfil
        # Vale a regra da mão direita
        for i in range(0,Ni-1):
            tri = np.append(tri,[[i, i+1, N-i-2]], axis=0)
        
        if z[0] != z[-1]:
            tri = np.append(tri,[[0, N-2, N-1]], axis=0)
            
        for i in range(1,Ni-val):
            tri = np.append(tri,[[i, N-i-2, N-i-1]], axis=0) 
        
        # Triangulação da superfície da casca da asa
        # Vale a regra da mão direita 
        # A triangulação da "casca" da asa é específica para o ordenamento dos perfis
        # dado pela funçao "tornar_3D()" 
        for j in lista:
            if z[0] != z[-1]:
                tri = np.append(tri,[[j, j+N-1, j+N],[j+N, j+N-1, j-1]], axis=0)   
            
            for i in range(N-1):
                tri = np.append(tri, [[j+i, j+N+i, j+i+1]], axis=0)
                tri = np.append(tri, [[j+i+1, j+N+i, j+N+i+1]], axis=0)
        
        if self.espelhar:
            # Triangulação do 5º perfil
            # É importante lembrar que agora vale a regra da mão esquerda
            for i in range(0,Ni-1):
                tri = np.append(tri,[[4*N+i, 5*N-i-2, 4*N+i+1]], axis=0)
            
            if z[0] != z[-1]:
                tri = np.append(tri,[[4*N, 5*N-1, 5*N-2]], axis=0)
                
            for i in range(1,Ni-val):
                tri = np.append(tri,[[4*N+i, 5*N-i-1, 5*N-i-2]], axis=0)
        else:
            # Triangulação do 3º perfil
            # É importante lembrar que agora vale a regra da mão esquerda
            for i in range(0,Ni-1):
                tri = np.append(tri,[[2*N+i, 3*N-i-2, 2*N+i+1]], axis=0)
            
            if z[0] != z[-1]:
                tri = np.append(tri,[[2*N, 3*N-1, 3*N-2]], axis=0)
                
            for i in range(1,Ni-val):
                tri = np.append(tri,[[2*N+i, 3*N-i-1, 3*N-i-2]], axis=0)           
            
        return tri
        
    def gerar(self):
        coords_raiz, coords_ponta = self.coords_perfis()   
        
        perfil_central = self.rotacao2d(self.a, self.translacao(0, 0, self.escala(self.c1, coords_raiz)))
        perfil_trans = list(perfil_central)
        perfil_ponta = self.rotacao2d(self.a + self.wo, self.translacao(self.off, 0, self.escala(self.c2, coords_ponta)))
        
        # Para Debug
#        self.grafico2d(perfil_central, self.c1, self.a, 'perfil raiz ')    
#        self.grafico2d(perfil_trans, self.c1, self.a, 'perfil trans ')
#        self.grafico2d(perfil_ponta, self.c2, self.a + self.wo, 'perfil ponta ')
        
        tri = np.array([[0.0, 0.0, 0.0]])
        x, y, z, N = self.tornar_3D(perfil_central, perfil_trans, perfil_ponta)
                    
        tri = self.triangulacao(z, N, tri)
        tri = tri[1:]
    
        return x,y,z,tri
         

class Fuselagem3D(object):
    """
    Essa classe cuida da triangulação da fuselagem.
    """    
    
    def __init__(self, h_bico,larg_bico, comp_bico, 
                 h_ccarga, larg_ccarga, 
                 comp_cima_ccarga, comp_baixo_ccarga):
        
        # Definição do bico do avião
        self.h_bico = h_bico
        self.larg_bico = larg_bico
        self.comp_bico = comp_bico
        
        # Definição do compartimento de carga do avião
        self.h_ccarga = h_ccarga
        self.larg_ccarga = larg_ccarga
        self.comp_cima_ccarga = comp_cima_ccarga
        self.comp_baixo_ccarga = comp_baixo_ccarga
        
    def tornar_3D(self):
        """
        Gera a lista de coordenadas da fuselagem.
        
        A origem dos eixos de encontra no centro da parte de cima do 
        compartimento de carga.
        
        O avião é orientado com a fuselagem apontando para o sentido 
        negativo de x e a envergadura na direção y.
        """
        h_bico = self.h_bico
        larg_bico = self.larg_bico
        comp_bico = self.comp_bico
        h_ccarga = self.h_ccarga
        larg_ccarga = self.larg_ccarga
        comp_cima_ccarga = self.comp_cima_ccarga
        comp_baixo_ccarga = self.comp_baixo_ccarga        
        
        # A fuselagem é composta por 14 pontos
        dim = 14
        
        # O argumento dim+2 se deve ao fato de o primeiro e o
        # último pontos das laterais da fuselagem serem o mesmo
        # ponto dos dois lados da mesma, o que adiciona 2 pontos.
        N = dim+2
        x = np.zeros(N)
        y = np.zeros(N)
        z = np.zeros(N)
        
        # Coordenadas x
        x[0] = x[7] = comp_cima_ccarga-comp_baixo_ccarga/2.0
        x[1] = x[6] = comp_baixo_ccarga/2.0
        x[2] = x[5] = -comp_baixo_ccarga/2.0
        x[3] = x[4] = -(comp_bico+comp_baixo_ccarga/2.0)
        x[8:] = x[:8]        
        
        # Coordenadas y
        y[0] = y[1] = y[2] = y[5] = y[6] = y[7] = larg_ccarga/2.0
        y[3] = y[4] = larg_bico/2.0
        y[8:] = -y[:8]
        
        # Coordenadas z
        z[0] = z[4] = z[5] = z[6] = z[7] = 0.0
        z[1] = z[2] = -h_ccarga
        z[3] = -h_bico
        z[8:] = z[:8]
        
        return x,y,z,N
        
    def triangulacao_lado(self, N_total, tri): 
        # São dois lados
        N = N_total//2
        Ni = N//2 
        
        if Ni == N/2.0: # Caso com número par de pontos
            val = 1
        else:
            val = 0
       
        # Triangulacao do 1º lado
        # Vale a regra da mão direita
        for i in range(0,Ni-1):
            tri = np.append(tri,[[i, i+1, N-i-2]], axis=0)
            
        for i in range(1,Ni-val):
            tri = np.append(tri,[[i, N-i-2, N-i-1]], axis=0) 
        
        # Triangulação do 2º lado
        # É importante lembrar que agora vale a regra da mão esquerda
        for i in range(0,Ni-1):
            tri = np.append(tri,[[N+i, 2*N-i-2, N+i+1]], axis=0)
            
        for i in range(1,Ni-val):
            tri = np.append(tri,[[N+i, 2*N-i-1, 2*N-i-2]], axis=0)
                
        return tri
    
             
    def triangulacao_casca(self, N_total, tri):
        """
        O conceito de casca da fuselagem segue o mesmo princípio da casca da asa
        A triangulação da "casca" da fuselagem é específica para o ordenamento
        dado pela funçao "tornar_3D()"   
        """
        # São dois lados        
        N = N_total//2
        
        # Triangulação da superfície da casca da asa
        # Vale a regra da mão direita           
        for i in range(N-1):
            tri = np.append(tri, [[i, N+i, i+1]], axis=0)
            tri = np.append(tri, [[i+1, N+i, N+i+1]], axis=0)
            
        return tri        
    
    def gerar(self):
        tri = np.array([[0.0, 0.0, 0.0]])
        x, y, z, N = self.tornar_3D()
                    
        tri = self.triangulacao_lado(N, tri)
        tri = self.triangulacao_casca(N, tri)
        tri = tri[1:]
        return x,y,z,tri
        
class Tailboom3D(Fuselagem3D):
    """
    Cuida da triangulação do tailboom treliçado
    """     
     
    def __init__(self, comp_baixo_ccarga, comp_cima_ccarga, h_ccarga,
                 comp_boom, esp_boom, h_cauda):        
        self.h_ccarga = h_ccarga
        self.comp_cima_ccarga = comp_cima_ccarga
        self.comp_baixo_ccarga = comp_baixo_ccarga
        self.comp_boom = comp_boom
        self.esp_boom = esp_boom
        self.h_cauda = h_cauda
        
    def tornar_3D(self):
        """
        Gera a lista de coordenadas do tailboom.
        
        A origem dos eixos de encontra no centro da parte de cima do 
        compartimento de carga.
        
        O avião é orientado com a fuselagem apontando para o sentido 
        negativo de x e a envergadura na direção y.
        """
        h_ccarga = self.h_ccarga
        comp_cima_ccarga = self.comp_cima_ccarga
        comp_baixo_ccarga = self.comp_baixo_ccarga 
        comp_boom = self.comp_boom
        esp_boom = self.esp_boom
        h_cauda = self.h_cauda
        
        # Cada tailboom é composto por 6 pontos
        dim = 6
        
        # O argumento dim+2 se deve ao fato de o primeiro e o
        # último pontos das laterais da fuselagem serem o mesmo
        # ponto dos dois lados da mesma, o que adiciona 2 pontos.
        N = dim+2
        x = np.zeros(N)
        y = np.zeros(N)
        z = np.zeros(N)
        
        # Coordenadas x
        x[0] = x[3] = comp_boom+comp_baixo_ccarga/2.0
        x[1] = comp_baixo_ccarga/2.0
        x[2] = comp_cima_ccarga-comp_baixo_ccarga/2.0
        x[4:] = x[:4]        
        
        # Coordenadas y
        y[0] = y[1] = y[2] = y[3] = esp_boom/2.0
        y[4:] = -y[:4]
        
        # Coordenadas z
        z[0] = z[3] = h_cauda
        z[2] = 0.0 
        z[1] = -h_ccarga
        z[4:] = z[:4]
        
        return x,y,z,N
             
def escrever_ascii_stl(arquivo, *objetos):   
    """
    Função que escreve o arquivo .stl

    objetos = (obj_1,obj_2,...,obj_n)
    obj_i = (nome, x, y, z, tri)
    
    nome : nome do objeto a ser gerado
    x,y,z : coordenadas dos pontos
    tri : array definindo a ordem da triangulação
    """     
    # Abre o arquivo .stl
    arquivo = open(str(arquivo)+'.stl', 'w')
    
    for i in range(len(objetos)):    
        objeto, x, y, z, tri = objetos[i]
        arquivo.write('solid %s\n' %(str(objeto)))    
        
        for i in range(len(tri)):
            # Calculando o vetor normal
            x0 = x[tri[i][0]]
            x1 = x[tri[i][1]]
            x2 = x[tri[i][2]]
            
            y0 = y[tri[i][0]]
            y1 = y[tri[i][1]]
            y2 = y[tri[i][2]]
            
            z0 = z[tri[i][0]]
            z1 = z[tri[i][1]]
            z2 = z[tri[i][2]]
            
            AB = np.array([x1-x0, y1-y0, z1-z0])
            AC = np.array([x2-x0, y2-y0, z2-z0])
            p_vetorial = np.cross(AB, AC)
            n = p_vetorial/np.linalg.norm(p_vetorial)
                    
            if np.isnan(np.linalg.norm(n)) == False:
                # Escrevendo faces
                arquivo.write('  facet normal %e %e %e\n' %(n[0],n[1],n[2]))
                arquivo.write('    outer loop\n')
                arquivo.write('      vertex %e %e %e\n' %(x0, y0, z0))
                arquivo.write('      vertex %e %e %e\n' %(x1, y1, z1))
                arquivo.write('      vertex %e %e %e\n' %(x2, y2, z2))
                arquivo.write('    endloop\n')
                arquivo.write('  endfacet\n')
                
        arquivo.write('endsolid %s\n' %(str(objeto)))
        
    arquivo.close()


def main(nome_arquivo, args_asa, args_eh, args_ev, args_fus, args_boom, inercia = False):

    perfil_asa_raiz, perfil_asa_ponta, iw, c_raiz_w, c_ponta_w,\
    c_mac_w, dist_cg_raiz_w,dist_trans_w, semi_bw, washout_w,\
    offset_w = args_asa
    
    perfil_eh, ieh, c_raiz_eh, c_ponta_eh, c_mac_eh,\
    dist_trans_eh, semi_beh, washout_eh, offset_eh = args_eh

    perfil_ev, iev, c_raiz_ev, c_ponta_ev, c_mac_ev,\
    dist_trans_ev, bev, washout_ev, offset_ev = args_ev    
    
    h_bico, larg_bico, comp_bico,h_ccarga,\
    larg_ccarga, comp_cima_ccarga,comp_baixo_ccarga = args_fus
    
    comp_boom, esp_boom, h_cauda = args_boom
   
    # ASA
    x_asa, y_asa, z_asa, tri_asa = Asa3D(perfil_asa_raiz, perfil_asa_ponta, iw, washout_w, offset_w, 
                                         c_raiz_w, dist_trans_w, c_ponta_w, semi_bw,
                                         c_mac_w,True).gerar()    
    x_asa = x_asa - dist_cg_raiz_w + 0.25*c_mac_w
    obj_asa = ('asa',x_asa,y_asa,z_asa,tri_asa)
    
    #EH
    x_eh, y_eh, z_eh, tri_eh = Asa3D(perfil_eh, perfil_eh, -ieh, washout_eh, offset_eh, 
                                     c_raiz_eh, dist_trans_eh, c_ponta_eh, semi_beh,
                                     c_mac_eh,True).gerar()  
    x_eh = x_eh + comp_baixo_ccarga/2.0 + comp_boom
    z_eh = -z_eh + h_cauda
    y_eh = -y_eh
    obj_eh = ('eh',x_eh,y_eh,z_eh,tri_eh)
    
    #EV
    x_ev, z_ev, y_ev, tri_ev = Asa3D(perfil_ev, perfil_eh, iev, washout_ev, offset_ev, 
                                     c_raiz_ev, dist_trans_ev, c_ponta_ev, bev,
                                     c_mac_ev,False).gerar()     
    x_ev = x_ev + comp_baixo_ccarga/2.0 + comp_boom
    z_ev = -z_ev + h_cauda
    obj_ev = ('ev',x_ev,y_ev,z_ev,tri_ev)
    
    # FUSELAGEM
    fus =  Fuselagem3D(h_bico,larg_bico, comp_bico, 
                       h_ccarga, larg_ccarga, 
                       comp_cima_ccarga, comp_baixo_ccarga)                     
    x_fus, y_fus, z_fus, tri_fus = fus.gerar()  
    obj_fus = ('fuselagem',x_fus,y_fus,z_fus,tri_fus)
    
    #TAILBOOM
    tailboom = Tailboom3D(comp_baixo_ccarga, comp_cima_ccarga, h_ccarga,
                          comp_boom, esp_boom, h_cauda)
    x_boom, y_boom, z_boom, tri_boom = tailboom.gerar()
    obj_boom1 = ('tailboom1',x_boom, y_boom+(larg_ccarga-esp_boom)/2.0, z_boom, tri_boom)
    obj_boom2 = ('tailboom2',x_boom, y_boom-(larg_ccarga-esp_boom)/2.0, z_boom, tri_boom)
    
    # ARQUIVO .STL
    if inercia:
        executarNaPasta('Graficos/Inercia')(escrever_ascii_stl)(nome_arquivo+'_asa', obj_asa)
        executarNaPasta('Graficos/Inercia')(escrever_ascii_stl)(nome_arquivo+'_eh',obj_eh)
        executarNaPasta('Graficos/Inercia')(escrever_ascii_stl)(nome_arquivo+'_ev', obj_ev)
        executarNaPasta('Graficos/Inercia')(escrever_ascii_stl)(nome_arquivo+'_fus', obj_fus)
        executarNaPasta('Graficos/Inercia')(escrever_ascii_stl)(nome_arquivo+'_boom1', obj_boom1)
        executarNaPasta('Graficos/Inercia')(escrever_ascii_stl)(nome_arquivo+'_boom2', obj_boom2)
    else:
        executarNaPasta('Graficos/STL3D')(escrever_ascii_stl)(nome_arquivo, obj_asa, obj_eh, obj_ev, obj_fus, obj_boom1, obj_boom2)


#nome_raiz,       nome_ponta,       alpha, washout,    offset,    corda_raiz, distancia1,    corda_ponta, distancia2, corda_mac, espelhar   
#perfil_asa_raiz, perfil_asa_ponta, iw,    washout_w,  offset_w,  c_raiz_w,   dist_trans_w,  c_ponta_w,   semi_bw,    c_mac_w,   True
#perfil_eh,       perfil_eh,       -ieh,   washout_eh, offset_eh, c_raiz_eh,  dist_trans_eh, c_ponta_eh,  semi_beh,   c_mac_eh,  True

if __name__ == '__main__':
    args_asa = ['S1223 MOD2015', 'S1223 MOD2015', 2.86, 0.448107662504, 0.291765181867,
                '?', 0.37447265625, 0.99834, -3.0,
                0.0356087340335,-0.201,0.12]
    args_eh = ['NACA 0011', 1.17, 0.288571, 0.288571,'?',
               0.0, 0.5650875, 2.8, 0.0,0.425,0.1855]

    args_ev = ['NACA 0011', 0.0, 0.288571, 0.281557,'?',
               0.0, 0.122112, 0.0, 0.0,0.425,0.1855]

    def f(nome_arquivo, args_asa, args_eh, args_ev):

        perfil_asa_raiz, perfil_asa_ponta, iw, c_raiz_w, c_ponta_w,\
        c_mac_w, dist_trans_w, semi_bw, washout_w,\
        offset_w,x_w,z_w = args_asa
        
        perfil_eh, ieh, c_raiz_eh, c_ponta_eh, c_mac_eh,\
        dist_trans_eh, semi_beh, washout_eh, offset_eh,xeh,zeh = args_eh

        perfil_ev, iev, c_raiz_ev, c_ponta_ev, c_mac_ev,\
        dist_trans_ev, bev, washout_ev, offset_ev,xev,zev = args_ev    
       
        # ASA
        x_asa, y_asa, z_asa, tri_asa = Asa3D(perfil_asa_raiz, perfil_asa_ponta, iw, washout_w, offset_w, 
                                             c_raiz_w, dist_trans_w, c_ponta_w, semi_bw,
                                             c_mac_w,True).gerar()    
        x_asa = x_asa + x_w
        z_asa = z_asa + z_w
        obj_asa = ('asa',x_asa,y_asa,z_asa,tri_asa)
        
        #EH
        x_eh, y_eh, z_eh, tri_eh = Asa3D(perfil_eh, perfil_eh, -ieh, washout_eh, offset_eh, 
                                         c_raiz_eh, dist_trans_eh, c_ponta_eh, semi_beh,
                                         c_mac_eh,True).gerar()  
        x_eh = x_eh + xeh
        z_eh = z_eh + zeh
        obj_eh = ('eh',x_eh,y_eh,z_eh,tri_eh)
        
        #EV
        x_ev, z_ev, y_ev, tri_ev = Asa3D(perfil_ev, perfil_eh, iev, washout_ev, offset_ev, 
                                         c_raiz_ev, dist_trans_ev, c_ponta_ev, bev,
                                         c_mac_ev,False).gerar()     
        x_ev = x_ev + xev
        z_ev = -z_ev + zev
        obj_ev = ('ev',x_ev,y_ev,z_ev,tri_ev)
        
        
        executarNaPasta('Graficos/STL3D')(escrever_ascii_stl)(nome_arquivo, obj_asa, obj_eh, obj_ev)
    
    f('Aeronave2016', args_asa, args_eh, args_ev)
