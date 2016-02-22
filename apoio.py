# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 19:19:22 2016

@author: Diego Chou
"""

print("\nCarregando modulos de 'Apoio'...")

try:
    import numpy as np
    from os import chdir, makedirs, listdir, getcwd
    from os.path import dirname, realpath, isfile, join
    from bisect import bisect_right, bisect_left
    print("Modulos de 'Apoio' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Apoio'\n")
    raise

def baricentro(p1,p2,p3):
    """
    Calculo do baricentro de um triangulo.
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    return (p1+p2+p3)/3.0

def inercia_ponto(p,p_ref,M):
    """
    Cálculo de inércia de uma massa pontual em função de um ponto de referência
    
    INPUT:
    -----
    * p = (x,y,z)
    * p_ref - ponto de referência
    * M - massa pontual
    """
    X,Y,Z = p_ref - p
    J = np.zeros((3,3))
    
    # Eixos paralelos
    J[0,0] += (Y**2 + Z**2)*M
    J[1,1] += (X**2 + Z**2)*M
    J[2,2] += (X**2 + Y**2)*M
    
    return J

def inercia_tri(p1,p2,p3,p_ref,ro):
    """
    Cálculo de inércia de um triângulo em função de um ponto de referência
    
    INPUT:
    -----
    * p_i = (x,y,z)
    * vertices - p1, p2, p3
    * ro - densidade superficial
    * p_ref - ponto de referência
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    p_ref = np.array(p_ref)    
    
    bar = baricentro(p1,p2,p3)
    X,Y,Z = p_ref - bar
    
    I = np.identity(3)  
    S = (1.0/24)*(np.ones((3,3))+I)
    V = np.matrix([p1,p2,p3])
    a = np.linalg.norm(np.cross(p2-p1,p3-p1))
    C = a*V.T*S*V
    
    # Inercia em torno do baricentro
    J = float(C.trace())*I-C
    
    # Eixos paralelos
    J[0,0] += (Y**2 + Z**2)*0.5*a
    J[1,1] += (X**2 + Z**2)*0.5*a
    J[2,2] += (X**2 + Y**2)*0.5*a
    
    # Transformação de área para massa
    J *= ro
    
    return J
    
def inercia_cubo(M,L):
    """
    Cálculo de inércia de um cubo em seus eixos principais, 
    relativo ao seu baricentro.
    
    INPUT:
    -----
    * M - massa do cubo
    * L - lado do cubo
    """
    return np.identity(3)*M*L**2/6.0
    
def tri(N):
    """
    Função de triangulação de um cortorno fechado com
    N pontos. O primeiro e o último pontos devem ser coincidentes.
    """
    Ni = N//2
    tri = np.array([[0,0,0]])
    
    for i in range(0,Ni-1):
        tri = np.append(tri,[[i, N-i-2, i+1]], axis=0)

    if Ni == N/2.0: # Caso com número par de pontos
        for i in range(1,Ni-1):
            tri = np.append(tri,[[i, N-i-1, N-i-2]], axis=0) 
            
    else: # Caso com número ímpar de pontos         
        for i in range(1,Ni): 
            tri = np.append(tri,[[i, N-i-1, N-i-2]], axis=0)
    
    return tri[1:] 

def executarNaPasta(string):
   """
   Muda o diretório de execução para o endereço especifidado pela string
   """
   def nested(func):
       def wrapper(*args, **kwargs):
            _cwd = getcwd()
            try:
                chdir(string)
            except OSError:
                # O makedirs cria diretórios recursivamente
                # mkdir('dir1/dir2')    -> ERRO
                # makedirs('dir1/dir2') -> OK
                makedirs(string)
                chdir(string)
            result = func(*args, **kwargs)
            chdir(_cwd)
            return result 
       return wrapper
   return nested
   
def interplinear(x1,y1,x2,y2,x):   
    if x1<=x<=x2 or x2<=x<=x1:
        try:
            y = ((x-x1)*y2+(x2-x)*y1)/(x2-x1)
        except ZeroDivisionError:
            y = y1
    else:
        print('\nx1 = %f' %(x1))
        print('x =  %f' %(x))
        print('x2 = %f' %(x2))
        raise ValueError("O ponto está fora do intervalo")
    return y 
  
def intervalo(lista, valor):
    lista = sorted(lista)
#    print lista
    if valor < lista[0] or valor > lista[-1]:
        print('\nINPUT:                    %f' %(valor))
        print('Intervalo de valores:    %s' %([lista[0],lista[-1]]))
        raise ValueError('Input fora do intervalo de valores')
        
    def achar_maior_antes():
        """
        Acha o maior valor da lista menor ou igual a x
        """
        i = bisect_right(lista,valor)
        if lista[i-1] == lista[-1]:
            return lista[i-2]
        else:
            return lista[i-1]
        
    def achar_menor_depois():
        """
        Acha o menor valor da lista maior ou igual a x
        """
        i = bisect_left(lista,valor)
        if lista[i] == lista[0]:
            return lista[i+1]
        elif (valor in lista) and (valor != lista[-1]):
            return lista[i+1]
        else:            
            return lista[i]
    
    valor_antes = achar_maior_antes()
    valor_depois = achar_menor_depois()
            
    return valor_antes, valor_depois
    
def lerDir(pasta, name, local = __file__):
    diretorio = dirname(realpath(local))
    diretorio = diretorio + pasta
    lista = []
    for i in listdir(diretorio):
        if isfile(join(diretorio,i)) and name+'_' in i:
            lista.append(i)
    return sorted(lista)
      
def lerDirPerfil_Re_Mach(name,local=__file__):
    lista = lerDir('/Perfis', name, local) 
    listaRe = []
    listaMach = []
    for i in range(len(lista)):
        listaRe.append(int(float(lista[i][len(name)+6:len(name)+11])*10**6))
        listaMach.append(float(lista[i][len(name)+13:len(name)+17]))
    listaRe = sorted(set(listaRe))
    listaMach = sorted(set(listaMach))   
    
    return listaRe, listaMach