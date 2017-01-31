# -*- coding: utf-8 -*-
"""
Created on Mon Sep 01 11:07:04 2014

@author: Mateus Pereira, Diego Chou
"""

print("\nCarregando modulos de 'Evolutionary'...")

try:
    import random
    from deap import algorithms, base, creator, tools
    import matplotlib.pyplot as plt
    import numpy as np
    import multiprocessing as mtp
    import sys
    import time
    import os
    import glob
#    import networkx
    print("Modulos de 'Evolutionary' foram carregados com sucesso!")
except ImportError:
    print("ERRO ao importar para 'Evolutionary'\n")
    raise

from apoio import salvarPrint, executarNaPasta
import avaliador

# Lista de variaveis:
#    
# Probabilidade de ocorrer mutacao em um individuo      pbmut       INPUT
# Probabilidade de ocorrer crossover                    pbcros      INPUT
# Numero de geracoes de individuos                      ngen        INPUT
# Limites inferiores das variaveis                      low         INPUT
# Limites superiores das variaveis                      up          INPUT
# Numero de individuos em cada geracao                  mu          INPUT
# Numero de novos individuos gerados a cada geracao     lamb        INPUT
# Numero de individuos armazenados no Hall Of Fame      size_hof    INPUT

# Velocidade de análise
vel = 20

# Ângulo de ataque de corrida do avião
alfa = 0.0

# Variável para a contagem de individuos
numero = 1

# Os perfis das asas são constantes durante a iteraçao
perfilr_asaf = perfilp_asaf = 'S1223 MOD2015'
perfilr_asat = perfilp_asat = 'NACA 0011'
perfil_ev = 'NACA 0011'


def pontuacao(individual):
    """
    Calcula a pontuacao de cada individuo
    IMPORTANTE: ESSA FUNCAO DEVE RETORNAR UMA TUPLA
    ----------  POR CAUSA DA ESTRUTURA INTERNA DO DEAP
    """
    global numero
    print("\n\n\nAVALIACAO DO INDIVIDUO %i" %(numero))

    dz_asas,x_motor,\
    x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,\
    x_bf_asat,c_asat,ang_asat,epsilon_asat,\
    c_ev = individual
    
    name = 'ind%d' %(numero)
    
    try:
        #with esconderPrint():
        A = avaliador.Avaliador2016(name,dz_asas,alfa,vel,x_motor, 
                                    x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                                    x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                                    c_ev,perfil_ev, p = False)
        score = A.pontuacao
        print('\nAeronave gerada com sucesso.')
    except Exception, e:
        print('\n%s' %(e))
        print("Erro na avaliacao. Pontuacao zerada.")
        score = 0

    if numero%10 == 0:
        files = glob.glob('Runs/*')
        time.sleep(2)
        for f in files:
            os.remove(f)
            
    numero += 1
    return score, #Retorna tupla


class Evolutionary(object):
    ngen = 50
    pbmut = 0.5
    pbcros = 0.5
    mu = 25
    lamb = 12
    #     dz_asas x_motor x_ba_asaf c_asaf ang_asaf epsilon_asaf x_bf_asat c_asat ang_asat epsilon_asat c_ev
    low = [0.05,  -1.0,   -1.0,     0.25,  0.0,     -5.0,        0.5,      0.25,  -5.0,    -5.0,        0.1]
    up =  [0.5,    0.0,    0.5,     0.5,   5.0,      0.0,        1.0,      0.5,    5.0,     5.0,        0.3]
    size_hof = 10
    
    bool_limites = False
    
    def __init__(self):
        # A ordem dos limites eh definida pela posicao das variaveis na lista
        # O hof armazena os melhores individuos dentre todas as iteracoes          
        self.hof = tools.HallOfFame(maxsize=self.size_hof)
        
        

    @salvarPrint('Output/melhor_individuo_%s_%s.txt' %(time.strftime("%Y %m %d"), 
                                                       time.strftime("%H.%M")))        
    def calculo_final(self,individual):
        """
        Refaz o cálculo completo do melhor avião
        """
        dz_asas,x_motor,\
        x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,\
        x_bf_asat,c_asat,ang_asat,epsilon_asat,\
        c_ev = individual
        
        name = 'best'
    
        print("\n\n\nREAVALIACAO DO MELHOR INDIVIDUO")        
        
        A = avaliador.Avaliador2016(name,dz_asas,alfa,vel,x_motor, 
                                    x_ba_asaf,c_asaf,ang_asaf,epsilon_asaf,perfilr_asaf, perfilp_asaf,
                                    x_bf_asat,c_asat,ang_asat,epsilon_asat,perfilr_asat, perfilp_asat,
                                    c_ev,perfil_ev, p = True)
                                    
        files = glob.glob('Runs/ind*')
        for f in files:
            os.remove(f)
                                             
        
    def mutate_VELHO(self, ind):
        """
        Gera mutacao com uma probabilidade indpbmut
        de ocorrer mutacao em cada variavel do indviduo
        """
        for i, xl, xu in zip(xrange(len(ind)), self.low, self.up):
            if random.random() < 0.1: #inserir indpbmut
                ind[i] = random.uniform(xl, xu)
        return ind,
        
        
    def mutate(self, ind):
        """
        Gera mutacao em uma variavel do individuo
        """
        k = random.randint(0, len(ind)-1)
#        print ind
#        print sum(ind), ind[k]
        xl = self.low[k]
        xu = self.up[k]
        ind[k] = random.uniform(xl, xu)
#        print sum(ind), ind[k]
#        print 30*'-'
        return ind,
        
        
    def individual(self):
        """
        Cria um novo individuo com valores aleatorios dentro 
        do intervalo de variancia
        """        
        ind = []
        for xl, xu in zip(self.low, self.up):
            ind.append(random.uniform(xl, xu))
        return ind
        
        
    def limites(self, pop):
        """
        Modifica os intervalos de variancia se a media do hof 
        estiver no 25% final ou inicial do intervalo
        """
        if self.bool_limites == True:
            print('\n\n%s' %(30*'-'))
            print('Intervalos variaveis\n')
            for i in range(len(self.low)):
                xl = self.low[i]
                xu = self.up[i]
                p_inf = xl + 0.25*(xu-xl)
                p_sup = xu - 0.25*(xu-xl)
    
                media = 0
                for j in range(len(self.hof)):
                    media += self.hof[j][i]
                media /= float(j)
                
                if media < p_inf:
                    self.low[i] = xl - 0.75*(xu-xl)
                    self.up[i] = p_inf
                    
                elif media > p_sup:
                    self.low[i] = p_sup
                    self.up[i] = xu + 0.75*(xu-xl)
                    
                print(media)
            
            print('low : %s' %(self.low))
            print('up :  %s' %(self.up))
            print(' %s' %(30*'-'))
            return
        else:
            print('\n\n%s' %(30*'-'))
            print('Intervalos fixos')
            print('%s' %(30*'-'))
            return
            
        
    def loop(self, log=False):
        """
        Realiza o loop genetico
        """
        creator.create('FitnessMax', base.Fitness, weights=(1.0,))
        creator.create('Individual', list, fitness=creator.FitnessMax)
  
        toolbox = base.Toolbox()
        
        toolbox.register('individual', tools.initIterate, creator.Individual, self.individual)
        toolbox.register('population', tools.initRepeat, list, toolbox.individual, n=self.mu)
        toolbox.register('mutate', self.mutate)
        
        toolbox.register('evaluate', pontuacao)
        
        toolbox.register('mate', tools.cxTwoPoint)
        toolbox.register('select', tools.selBest)

        if sys.platform == 'linux2':
# O paralelismo so funciona em linux
            toolbox.register('map', mtp.Pool().map)
            
#        history = tools.History()
#        toolbox.decorate('mate', history.decorator)
#        toolbox.decorate('mutate', history.decorator)

        st = tools.Statistics(lambda ind: ind.fitness.values)
        st.register("Media", np.mean)
        st.register("DP", np.std)   
        st.register("Minimo", np.min)
        st.register("Maximo", np.max)
        
        # Os limites variam somente se essa opção constar no arquivo de configuração
        st.register("Limites", self.limites)

        pop = toolbox.population()
#        history.update(pop)
        pop, logbook = algorithms.eaMuPlusLambda(pop, toolbox, self.mu, self.lamb, self.pbcros, 
                                                 self.pbmut, self.ngen, stats=st, halloffame=self.hof,
                                                 verbose=False)
        
#        plt.figure()
##        graph = networkx.DiGraph(history.genealogy_tree)
##        graph = graph.reverse()     # Make the grah top-down
##        colors = [toolbox.evaluate(history.genealogy_history[i])[0] for i in graph]
##        networkx.draw(graph, node_color=colors)
#        gen_best = history.getGenealogy(self.hof[0])
#        graph = networkx.DiGraph(gen_best).reverse()
#        networkx.draw(graph)
#        plt.show()        
        
        if log: 
            return pop, logbook

        else:   
            return pop

    
    @executarNaPasta('Graficos/Evolutionary')
    def grafico(self, logbook):
        """
        Plota o grafico Pontuacao x Geracao
        """
        log = str(logbook).split()
        av = []
        ma = []
        for i in range(7, len(log), 7): 
#        de 7 em 7 porque ha 5 registros em tools.Statistics + gen e nevals

#            for j in range(len(log[i+2])):
#                if log[i+2][j] == '.': log[i+2] = log[i+2][:j] + log[i+2][j+1:]
##            os 2 loops sao para tirar o ponto do numero
#            for j in range(len(log[i+5])):
#                if log[i+5][j] == '.': log[i+5] = log[i+5][:j] + log[i+5][j+1:]

            for j in range(len(log[i+2])):
                if log[i+2][j] == ',': 
                    log[i+2] = log[i+2][:j] + '.' + log[i+2][j+1:]
#            os 2 loops sao para trocar a virgula do decimal do numero para um ponto
            for j in range(len(log[i+5])):
                if log[i+5][j] == ',': 
                    log[i+5] = log[i+5][:j] + '.' + log[i+5][j+1:]

            av.append(float(log[i+2]))
            ma.append(float(log[i+5]))

        gen = range(self.ngen+1)
        
        plt.figure()
        plt.grid('on')
        plt.title('Pontuacao por geracao')
        plt.plot(gen, av, 'r-', label='$P_{media}$')
        plt.plot(gen, ma, '--', label='$P_{max}$')
        plt.xlabel('$Gen$')
        plt.ylabel('$P$')
        plt.legend(loc='best')
        plt.savefig("evolutionary.png", bbox_inches='tight', dpi=200)
#        plt.show()


    def principal(self): 
        population, logbook = self.loop(True)
        self.grafico(logbook)
    
        @salvarPrint('Output/evolutionary_%s_%s.txt' %(time.strftime("%Y %m %d"), 
                                                       time.strftime("%H.%M")))
        def print_final():
            print('')
            print(logbook)        
    
            print('\nLimites Inferiores :')
            print(self.low)
    
            print('\nLimites Superiores :')
            print(self.up)
            
            hall = self.hof
            print('\nMelhores individuos :')
            for i in range(len(hall)):
                print(hall[i]) 
                print('')
                     
        
        print_final()       
        
        self.calculo_final(self.hof[0])
