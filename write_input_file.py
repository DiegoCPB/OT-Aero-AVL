# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 19:16:21 2016

@author: Diego Chou
"""

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

        self.create_file()

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
        * header = [mach,iYsym,iZsym,Zsym,Sref,Cref,Bref,Xref,Yref,Zref,CDp]
        """
        name = self.name
        
        mach, iYsym, iZsym, Zsym,\
        Sref, Cref, Bref,\
        Xref, Yref, Zref, CDp = self.header
        
        self.arquivo.write('%s\n' %(name))
        self.arquivo.write('%.1f\n' %(mach))
        self.arquivo.write('%d    %d    %.1f\n' %(iYsym, iZsym, Zsym))
        self.arquivo.write('%.1f    %.1f    %.1f\n' %(Sref, Cref, Bref))
        self.arquivo.write('%.1f    %.1f    %.1f\n' %(Xref, Yref, Zref))
        self.arquivo.write('%f\n\n' %(CDp))
        self.arquivo.write('#'+79*'='+'\n\n')
    
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
            * COORDS :   [Xle, Yle, Zle, Chord, Ainc, Nspanwise, Sspace]
            * AFILE :    airfoil filename
            * CLAF :     dCL/da scaling factor
        """
        for SURFACE in self.asas:
            self.arquivo.write('SURFACE\n')    
            self.arquivo.write('%s\n' %(SURFACE))
            VAL = self.asas.get(SURFACE)
            
            for i in range(len(VAL)):
                if i == 0:
                    self.arquivo.write('%f %f %f %f\n\n' %(VAL[0],VAL[1],
                                                         VAL[2],VAL[3]))
                                                         
                elif i == 1:
                    self.arquivo.write('YDUPLICATE\n')
                    self.arquivo.write('%f\n' %(VAL))
                    
                elif i == 2:
                    self.arquivo.write('SCALE\n')
                    self.arquivo.write('%f %f %f\n\n' %(VAL[0],VAL[1],
                                                         VAL[2]))
                    
                elif i == 3:
                    self.arquivo.write('TRANSLATE\n')
                    self.arquivo.write('%f %f %f\n\n' %(VAL[0],VAL[1],
                                                         VAL[2]))
                                                         
                elif i == 4:
                    self.arquivo.write('ANGLE\n')
                    self.arquivo.write('%f\n' %(VAL))
                    self.arquivo.write('#'+79*'-'+'\n')
                    
                elif i == 5:
                    for j in len(VAL):
                        self.arquivo.write('#'+79*'-'+'\n')
                        self.arquivo.write('SECTION\n')
                        PAR=VAL[j]
                        
                        if j == 0:
                            self.arquivo.write('%f %f %f %f %f %f %f\n\n' %(PAR[0],PAR[1],PAR[2],PAR[3], 
                                                                            PAR[4],PAR[5],PAR[6]))
                        elif j == 1:
                            self.arquivo.write('AFILE\n')
                            self.arquivo.write('%s\n' %(PAR))
                            
                        elif j == 2:
                            self.arquivo.write('CLAF\n')
                            self.arquivo.write('%s\n' %(PAR))
                            
            self.arquivo.write('#'+79*'='+'\n\n')
                            
class Input_mass(object):
    """
    Classe que escreve o arquivo de input da inércia do avião no formato
    do software AVL v3.35
    """
if __name__ == '__main__':
    header = [0.0,1,0,0.0,0.8,0.35,2.2,0.3,0.0,0.5,0.0]
    asas = {}
    f = Input_geometry('teste',header,asas)
                                        
                    
                    
                                         
                    
                    
                    
                
