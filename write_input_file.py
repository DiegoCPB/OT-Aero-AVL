# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 19:16:21 2016

@author: Diego Chou
"""

class Input_geometry(object):
    """
    Classe que escreve o arquivo de input geométrico no formato
    do softwae AVL v3.35
    """
    def __init__(self,name,header,asas,):
        # Header inputs        
        self.name = name
        self.header = header 
        
        # Surface and body data inputs
        self.asas = asas #inclui todas as superfícies sustentadoras

        self.create_file()

    def create_file(self):
        self.arquivo = open(self.name+'.avl','w')

        self.write_header()
#        self.write_surfaces()
        
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
        self.arquivo.write('%f\n' %(CDp))
        self.arquivo.write('#'+79*'='+'\n\n')
    
    def write_surfaces(self):
        """
        Parâmetros de superfícies.
    
        INPUT
        -----
        * asas = {surface_1, surface_2, ...}
            * surface_i = [DISC,CONPONENT,YDUPLICATE,SCALE,TRANSLATE,ANGLE,SECTION]
                * SECTION = [section_1, section_2, ...]
                    * section_i = {...}
        """
        for SURFACE in self.asas:
            self.arquivo.write('SURFACE\n')    
            self.arquivo.write('%s\n' %(SURFACE))
            
            for i in len(self.asas.get(SURFACE)):
                VAL = self.asas.get(SURFACE)
                if i == 0:
                    self.arquivo.write('%f %f %f %f\n\n' %(VAL[0],VAL[1],
                                                         VAL[2],VAL[3]))
                                                         
                if i == 1:
                    self.arquivo.write('COMPONENT\n')
                    self.arquivo.write('%f\n' %(VAL))
                   
                if i == 2:
                    self.arquivo.write('YDUPLICATE\n')
                    self.arquivo.write('%f\n' %(VAL))
                    
                if i == 3:
                    self.arquivo.write('SCALE\n')
                    self.arquivo.write('%f %f %f\n\n' %(VAL[0],VAL[1],
                                                         VAL[2]))
                    
                if i == 4:
                    self.arquivo.write('TRANSLATE\n')
                    self.arquivo.write('%f %f %f\n\n' %(VAL[0],VAL[1],
                                                         VAL[2]))
                                                         
                if i==5:
                    self.arquivo.write('ANGLE\n')
                    self.arquivo.write('%f\n' %(VAL))
                    
                if i == 6:
                    for j in len(VAL):
                        self.arquivo.write('SECTION\n')
                        lista=VAL[j]
                        self.arquivo.write('%f %f %f %f %f %f %f\n\n' %(lista[0],lista[1],lista[2],lista[3], 
                                                                        lista[4],lista[5],lista[6]))
                    
if __name__ == '__main__':
    header = [0.0,1,0,0.0,0.8,0.35,2.2,0.3,0.0,0.5,0.0]
    asas = {}
    f = Input_geometry('teste',header,asas)
                                                         
                    
                    
                                         
                    
                    
                    
                
