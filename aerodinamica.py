# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 23:55:59 2015

@author: Diego Chou
"""

print("\nCarregando modulos de 'Aerodinamica'...")
try:
    import subprocess as sp
    print("Modulos de 'Aerodinamica' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Aerodinamica'\n")
    raise

from apoio import issueCmd

def aerodinamica(name,vel,trim,p=False):                                                  
    """
    Analise aerodinâmica do avião gerado.
    Essa funçao executa o programa avl.exe
    """
    ps = sp.Popen(['avl.exe'],stdin=sp.PIPE,stdout=None,stderr=None)
            
    issueCmd(ps,'load Runs/%s.avl' %(name))
    issueCmd(ps,'mass Runs/%s.mass' %(name))
    issueCmd(ps,'mset 1')
    issueCmd(ps,'oper')
    issueCmd(ps,'a')
    issueCmd(ps,'a')
    issueCmd(ps,'%f' %(trim))
    issueCmd(ps,'m')
    issueCmd(ps,'v')
    issueCmd(ps,'%f' %(vel))
    issueCmd(ps,'')
    issueCmd(ps,'x')
