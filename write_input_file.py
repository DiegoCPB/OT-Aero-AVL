# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 19:16:21 2016

@author: Diego Chou
"""
import apoio

class Input_geometry(object):
    """
    Classe que escreve o arquivo de input geométrico no formato
    do software AVL v3.35
    """
    def __init__(self,name,header,asas):
        # Header inputs        
        self.name = name
        self.header = header 
        
        # Surface data inputs
        self.asas = asas

        # Escreve o arquivo de geometria
        self.create_file()

    @apoio.executarNaPasta('Runs')
    def create_file(self):
        self.arquivo = open(self.name+'.avl','w')
        self.write_header()
        self.write_surfaces()
        self.arquivo.close()
                
    def write_header(self):
        """
        Parâmetros do header
        
        INPUT
        -----
        * header = [mach,iYsym,iZsym,Zsym,Sref,Cref,Bref,Xref,Yref,Zref]
        """
        name = self.name
        
        mach, iYsym, iZsym, Zsym,\
        Sref, Cref, Bref,\
        Xref, Yref, Zref = self.header
        
        self.arquivo.write('%s\n' %(name))
        self.arquivo.write('%.1f\n' %(mach))
        self.arquivo.write('%d %d %.1f\n' %(iYsym, iZsym, Zsym))
        self.arquivo.write('%.1f %.1f %.1f\n' %(Sref, Cref, Bref))
        self.arquivo.write('%.1f %.1f %.1f\n' %(Xref, Yref, Zref))
    
    def write_surfaces(self):
        """
        Parâmetros de superfícies.
    
        INPUT
        -----
        * asas :         {surface_1, surface_2, ...}
         * surface_i :   [DISC, YDUPLICATE, SCALE, TRANSLATE, ANGLE, SECTION]
          * DISC :       [Nchord ,Cspace, Nspace, Sspace]
          * YDUPLICATE : Ydupl
          * SCALE :      [Xscale, Yscale, Zscale]
          * TRANSLATE :  [dX, dY, dZ]
          * ANGLE :      dAinc
          * SECTION :    [section_1, section_2, ...]
           * section_i : [COORDS,AFILE,CLAF]
            * COORDS :   [Xle, Yle, Zle, Chord, Ainc]
            * AFILE :    airfoil filename
            * CLAF :     dCL/da scaling factor
        """
        for SURFACE in self.asas:
            self.arquivo.write('#'+79*'='+'\n')
            self.arquivo.write('SURFACE\n')    
            self.arquivo.write('%s\n' %(SURFACE))
            SECT = self.asas.get(SURFACE)
            
            for i in range(len(SECT)):
                VAL = SECT[i]
                if i == 0:
                    self.arquivo.write('%d %.1f %d %.1f\n\n' %(VAL[0],VAL[1],
                                                                   VAL[2],VAL[3]))
                                                         
                elif i == 1:
                    self.arquivo.write('YDUPLICATE\n')
                    self.arquivo.write('%.2f\n\n' %(VAL))
                    
                elif i == 2:
                    self.arquivo.write('SCALE\n')
                    self.arquivo.write('%.2f %.2f %.2f\n\n' %(VAL[0],VAL[1],
                                                         VAL[2]))
                    
                elif i == 3:
                    self.arquivo.write('TRANSLATE\n')
                    self.arquivo.write('%.2f %.2f %.2f\n\n' %(VAL[0],VAL[1],
                                                         VAL[2]))
                                                         
                elif i == 4:
                    self.arquivo.write('ANGLE\n')
                    self.arquivo.write('%.2f\n' %(VAL))
                    
                elif i == 5:
                    for j in range(len(VAL)):
                        self.arquivo.write('#'+79*'-'+'\n')
                        self.arquivo.write('SECTION\n')
                        PAR=VAL[j]
                        for k in range(len(PAR)):
                            NUM = PAR[k]
                            if k == 0:
                                self.arquivo.write('%.2f %.2f %.2f %.2f %.2f %d %d \n\n' %(NUM[0],NUM[1],
                                                                                           NUM[2],NUM[3],
                                                                                           NUM[4],NUM[5],NUM[6]))
                            elif k == 1:
                                self.arquivo.write('AFILE\n')
                                self.arquivo.write('%s\n\n' %(NUM))
                                
                            elif k == 2:
                                self.arquivo.write('CLAF\n')
                                self.arquivo.write('%.3f\n' %(NUM))
                            
class Input_mass(object):
    """
    Classe que escreve o arquivo de input da inércia do avião no formato
    do software AVL v3.35
    """
if __name__ == '__main__':
    header = [0.0,1,0,0.0,0.8,0.35,2.2,0.3,0.0,0.5]
    asas = {'Asa_frontal': [[8, 1.0, 12, -2.0], 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, [[[-0.91999999999999993, 0.0, 0.083577871373829077, 0.5, 5.0], 'Perfis/S1223 MOD2015.dat', 0.99259071560308809], [[-0.067942265947765423, 1.0154428656544325, 0.131019448207603, 0.25, 2.0], 'Perfis/S1223 MOD2015.dat', 1.0978393244087321]]], 
            'Asa_traseira': [[8, 1.0, 12, -2.0], 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, [[[0.22737021437695154, 0.0, 0.43357787137382908, 0.3, 2.0], 'Perfis/S1223 MOD2015.dat', 1.0897423608814247], [[-0.067942265947765423, 1.0154428656544325, 0.131019448207603, 0.25, 2.0], 'Perfis/S1223 MOD2015.dat', 1.0978393244087321]]]}
    f = Input_geometry('teste',header,asas)
                                        
                    
                    
                                         
                    
                    
                    
                
