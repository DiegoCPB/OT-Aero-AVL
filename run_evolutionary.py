# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 14:17:09 2015

@author: diego
"""

import evolutionary as evo
from apoio import tempoDeExecucao

@tempoDeExecucao
def main():
    ev = evo.Evolutionary()
    ev.principal()

main()
