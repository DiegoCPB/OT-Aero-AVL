# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 23:55:59 2015

@author: Diego Chou
"""

print("\nCarregando modulos de 'Aerodinamica'...")
try:
    print("Modulos de 'Aerodinamica' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Aerodinamica'\n")
    raise


def aerodinamica(ps,name,vel,trim,p=False):                                                  
    """
    Analise aerodinâmica do avião gerado.
    Essa funçao executa o programa avl.exe
    """
