from guiparts import PlotPanel
import math            
import settings        
import time 
import wx
import matplotlib

#===================================================================================================
class EigenPlotPanel(PlotPanel):
    """Plots several lines in distinct colors."""
    
    XTOL = 0.2
    YTOL = 0.2
    PICK_TOLERANCE = 5
    DRAW_BUTTON = 1
    DEL_BUTTON = 3
    AXIS_MIN_MAX = [-2,2,-2,2]
    #------------------------------------------------------------------------------------------------                          
    def __init__( self, parent,**kwargs ):
        
        self.solver = parent.huckel_solver
        self.molecule = parent.sketch_pad.molecule
        self.pointer = parent.level_pointer
        self._bonding = []
        self.atom_bond_stack2 = []
        self.current_atom_type = 'C'
        self.set_focus_to_huck = None
        # initiate plotter
        PlotPanel.__init__( self, parent, **kwargs )

        self.back_from_dialog = False
        self.dbl_click_timer = None
        self.line = None
        self.lineXs = []
        self.lineYs = []
        self.mouseDown = False
        self._redrawFlag = False
        #self.SetColor( (255,255,255) )
        self.solver_update_needed = False

    #------------------------------------------------------------------------------------------------          
    def _onIdle( self, evt ):
        
        PlotPanel._onIdle(self,evt)
        if self._redrawFlag and self.GetParent().visual_mode.IsChecked():
            self._redrawFlag = False
            self.draw()
            
        
    #------------------------------------------------------------------------------------------------          
    def refreshFromHuckel(self):
        self._redrawFlag = False
        if self.GetParent().visual_mode.IsChecked():        
            self.draw()
        
    #------------------------------------------------------------------------------------------------              
    def initializePlot(self):
    
        self.subplot = self.figure.add_subplot( 111 )
        self.figure.suptitle('Orbitals')
        self.setPlotProperties()
        self.drawLegend()                    
        self.setPlotProperties()
        
    def setPlotProperties(self):
        #self.subplot.axis('auto')
        #self.subplot.set_aspect('auto',adjustable="datalim")
        self.subplot.set_xticklabels([""])
        self.subplot.set_xticks([])
        self.subplot.set_yticks([])
        self.subplot.set_yticklabels([""])
        
    
    def addNewAtom(self,x,y,size=10):
        if size > 0:
            fmt = 'or'
        else:
            fmt = 'ob'
        self.subplot.plot([x],[y],fmt,markersize=abs(size))#max(2,abs(size)))
    #------------------------------------------------------------------------------------------------                      
    def addNewBond(self,xcoords,ycoords):
        
        #if not self.molecule.getBond(connection[0],connection[1]):
        
        xs,ys = xcoords[0],ycoords[0]
        xf,yf = xcoords[1],ycoords[1]
        self.subplot.plot([xs,xf],[ys,yf],'r')
       

    #------------------------------------------------------------------------------------------------              
    def draw( self ):
        if not hasattr( self, 'subplot' ):
            self.initializePlot()            
        else:
            while self.subplot.lines:
                self.subplot.lines.pop()
            while self.subplot.texts:
                self.subplot.texts.pop()

        if len(self.molecule.atom_stack)>0 and self.GetParent().visual_mode.IsChecked():
            
            eigen_vecs = self.solver.eigen_vecs
            pointer = self.GetParent().level_pointer
            if pointer > len(eigen_vecs):
                pointer = len(eigen_vecs)-1

            magnitudes = [int(50*x) for x in eigen_vecs[pointer]]

            for bond in self.molecule.bond_stack:

                self.addNewBond(bond.XCoords(),bond.YCoords())

            self.xs = []
            self.ys = []
                
            for ii,atom in enumerate(self.molecule.atom_stack):
                x,y = atom.x,atom.y
                self.xs.append(x),self.ys.append(y)
                self.addNewAtom(x,y,magnitudes[ii])
            self.resize()
            
            self.relable_atoms()


        self.canvas.draw()
        
    def drawLegend(self):

        fmts = (('.w','Sign of eigenvector coefficient:'),('ro','Ci > 0'),('bo','Ci < 0'),)                    
        [self.subplot.plot([-999],[0],fmt[0],label=fmt[1],markersize=5) for fmt in fmts]
            
        font = matplotlib.font_manager.FontProperties(size=11)
        legend = self.subplot.legend(loc=8,ncol=3,prop=font,columnspacing=1,numpoints=1,markerscale=400)
        legend.draw_frame(False)
        
#------------------------------------------------------------------------------------------------                      
    def resize(self):
        
        #zoom out so atom occupies roughly 1/3 panel
        scale = 0.4
        min_x,max_x = min(self.xs),max(self.xs)
        min_y,max_y = min(self.ys),max(self.ys)
        
        dx = abs(max_x-min_x)*scale
        dy = abs(max_y-min_y)*scale
        ax = [min_x-dx,max_x+dx,min_y-dy,max_y+dy]
        self.subplot.axis('equal')
        self.subplot.axis(ax)
            

    def relable_atoms(self):

        self.subplot.texts = []
        
        delta = 0.1
        for ii in range(len(self.xs)):
            x,y = self.xs[ii],self.ys[ii]
            if x < 0 and y<0:
                ha,va='right','top'
                dx,dy = -delta,-delta
            elif x>0 and y<0:
                ha,va='left','top'
                dx,dy = delta,-delta
            elif x<0 and y>0:
                ha,va='right','bottom'
                dx,dy = -delta,delta
            else:
                ha,va='left','bottom'
                dx,dy= delta,delta
            self.subplot.text(x+dx,y+dy,str(ii+1),ha=ha,va=va,size='large',zorder=20,color='black')              
