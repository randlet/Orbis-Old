import time
from guiparts import PlotPanel
from molecule import Molecule,Atom,Bond
import matplotlib
import math            
import numpy
import settings        
import time 
import wx
import sys
from validators import NumberValidator
from editatombond import EditBondTypes,EditAtomTypes
#===================================================================================================
class AtomTypeDialog(wx.Dialog):
    
    def __init__(self,atom,atom_num):
        wx.Dialog.__init__(self,None,-1, 'Select Atom Properties for atom %d' %(atom_num+1))
        
        ok_button = wx.Button(self,wx.ID_OK,"OK")
        ok_button.SetDefault()
        cancel_button = wx.Button(self,wx.ID_CANCEL, "Cancel")

        a_txt = wx.StaticText(self, -1,u"Select Atom Type:", style=wx.ALIGN_CENTRE)
        
        self.atom_type = wx.ComboBox(self,-1,atom.sym,choices=Atom.ATOM_TYPES.keys(),style=wx.CB_DROPDOWN|wx.CB_READONLY,pos=(140,7))
        self.atom_type.SetToolTip(wx.ToolTip(Atom.ATOM_TYPES[self.atom_type.GetValue()]["description"]))
        self.Bind(wx.EVT_COMBOBOX,self.onAtomType,self.atom_type)
        
        hx_txt = wx.StaticText(self, -1,u"Set hx for atom %d:" % (atom_num+1), style=wx.ALIGN_CENTRE,pos = (10,40))
        self.hx = wx.TextCtrl(self, -1,str(atom.hx),validator=NumberValidator())
        self.hx.SetValue(str(atom.hx))
        

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        sizer = wx.GridSizer(2,2,10,10)
        sizer.Add(a_txt,flag=wx.ALIGN_RIGHT)
        sizer.Add(self.atom_type,flag=wx.ALIGN_LEFT)
        sizer.Add(hx_txt,flag=wx.ALIGN_RIGHT)
        sizer.Add(self.hx,flag=wx.ALIGN_LEFT)

        button_sizer = wx.BoxSizer()

        button_sizer.Add(ok_button,flag=wx.ALIGN_RIGHT)
        button_sizer.Add(cancel_button,flag=wx.ALIGN_RIGHT)
        
        main_sizer.Add(sizer,flag=wx.ALL,border=5)
        main_sizer.Add(button_sizer,-1,wx.ALL|wx.ALIGN_RIGHT,5)
        self.SetSizer(main_sizer)
        self.Fit()
        self.hx.SetFocus()
        self.hx.SelectAll()
        
    def onAtomType(self,event):
        atom_type = Atom.ATOM_TYPES[self.atom_type.GetValue()]
        self.hx.SetValue(str(atom_type[settings.h_delim]))
        self.atom_type.SetToolTip(wx.ToolTip(atom_type["description"]))
#===================================================================================================        
class BondTypeDialog(wx.Dialog):
    
    def __init__(self,bond,input_str):
        wx.Dialog.__init__(self,None,-1, 'Select Bond Properties')
        ok_button = wx.Button(self,wx.ID_OK,"OK")
        ok_button.SetDefault()
        cancel_button = wx.Button(self,wx.ID_CANCEL, "Cancel")
        
        self.k_xy_text = wx.StaticText(self, -1,input_str)
        self.k_xy = wx.TextCtrl(self, -1,str(bond.k_xy),validator=NumberValidator())
        self.k_xy.SetValue(str(bond.k_xy))


        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        sizer = wx.BoxSizer()
        sizer.Add(self.k_xy_text,flag=wx.ALIGN_RIGHT)
        sizer.Add(self.k_xy,flag=wx.ALIGN_LEFT)
        button_sizer = wx.BoxSizer()

        button_sizer.Add(ok_button,flag=wx.ALIGN_RIGHT)
        button_sizer.Add(cancel_button,flag=wx.ALIGN_RIGHT)
       
        main_sizer.Add(sizer,flag=wx.ALL,border=5)
        main_sizer.Add(button_sizer,-1,wx.ALL|wx.ALIGN_RIGHT,5)
        self.SetSizer(main_sizer)
        self.Fit()
        self.k_xy.SetFocus()
        self.k_xy.SelectAll()
#===================================================================================================
class MoleculePlotPanel(PlotPanel):
    """Plots several lines in distinct colors."""
    
    XTOL = 0.1
    YTOL = 0.1
    PICK_TOLERANCE = 5
    DRAW_BUTTON = 1
    DEL_BUTTON = 3
    AXIS_MIN_MAX = [-2,2,-2,2]
    TOL_RATIO = XTOL/(AXIS_MIN_MAX[1]-AXIS_MIN_MAX[0])
    
    #------------------------------------------------------------------------------------------------                          
    def __init__( self, parent,**kwargs ):

        
        self.tooltip = wx.ToolTip(tip='tip with a long %s line and a newline\n' % (' '*100))
        self.tooltip.Enable(False)
        self.tooltip.SetDelay(500)
        
        self.solver = parent.huckel_solver

        self.ignored_atom_types = []
        self.ignored_bond_types = []
        
        self.molecule = Molecule()
        self._bonding = []
        self.atom_bond_stack2 = []
        self.current_atom_type = "C"
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
        self.moved = False
        self._doRefreshFromHuckel = False
        #self.SetColor( (255,255,255) )
        self.solver_update_needed = False
        
        

    #------------------------------------------------------------------------------------------------                      

    def reset(self):
        #pop off all atoms and bonds and then tell the huckel solver it needs to be updated
        while self.subplot.lines:
            self.subplot.lines.pop()

            item2 = self.atom_bond_stack2.pop()
            if item2 in self.molecule.atom_stack:
                self.molecule.removeAtom(self.molecule.atom_stack.index(item2))
            elif item2 in self.molecule.bond_stack:
                bond = self.molecule.bond_stack[self.molecule.bond_stack.index(item2)]
                a = self.molecule.atom_stack.index(bond.a)
                b = self.molecule.atom_stack.index(bond.b)
                self.molecule.removeBond(a,b)

        self.relable_atoms()

        #self.solver_update_needed = True
        self._redrawFlag  = True
    #------------------------------------------------------------------------------------------------          
    def mouseDown(self,event):
        self.moved = False
        if event.inaxes and not self.back_from_dialog:
            self.mouseDown = True
            self.mouseDownX = event.xdata
            self.mouseDownY = event.ydata
            self.mouse_down_on_atom = self.onAtom(event.xdata,event.ydata)
        else:
            self.back_from_dialog = False
            self.mouse_down_on_atom = -1
            self.set_focus_to_huck = None
            self.line = None
    #------------------------------------------------------------------------------------------------              
    def onAtom(self,x,y):
        onAtom = -1
        for ii,atom in enumerate(self.molecule.atom_stack):
            xa,ya = atom.x,atom.y
            if abs(xa-x)< self.XTOL and abs(ya-y)< self.YTOL:
                onAtom = ii
                break
        return onAtom

    #------------------------------------------------------------------------------------------------      
    def newAtomRequired(self,x,y):
        unconnected_bond =  self.line != None
        click_on_bond_or_atom = self.set_focus_to_huck != None
        down_on_atom = self.mouse_down_on_atom >= 0
        return  not (click_on_bond_or_atom or unconnected_bond or down_on_atom)
    #------------------------------------------------------------------------------------------------                      
    def addNewAtom(self,x,y,hx=0.0,sym=""):
        if sym == "":
            sym = self.current_atom_type
        #make sure all atom types are present
        if sym not in Atom.ATOM_TYPES.keys():

            if sym not in self.ignored_atom_types and wx.MessageBox("This atom type (%s) doesn't currently exist in the Orbis data set. Create a new atom type? (Select No to use default (C) or Yes to create new atom type)" %(sym),"Create new atom type",style=wx.YES_NO) == wx.YES:
                dlg = EditAtomTypes()
                create_new = dlg.OnNew(None)

                if create_new == wx.ID_OK:
                    sym = dlg.atom_type.GetValue()
                    Atom.WriteAtomLib(dlg.atom_types)
                else:
                    self.ignored_atom_types.append(sym)
                    sym = "C"
                dlg.Destroy()
            else:
                self.ignored_atom_types.append(sym)
                sym = "C"
        hx = Atom.ATOM_TYPES[sym]['h']
        color = settings.COLOUR_LIST[sorted(Atom.ATOM_TYPES.keys()).index(sym)].lower()

        artist = self.subplot.plot([x],[y],'o',color=color,picker=self.PICK_TOLERANCE,markersize=10,zorder=20)[-1]

        self.molecule.addAtom(artist,x,y,hx,sym)
        self.atom_bond_stack2.append(self.molecule.atom_stack[-1])

    #------------------------------------------------------------------------------------------------                      
    def connectionMade(self):
        
        downAtom = self.onAtom(self.mouseDownX,self.mouseDownY)
        upAtom =  self.onAtom(self.mouseUpX,self.mouseUpY)

        if downAtom >= 0 and upAtom >= 0 and downAtom != upAtom:
            return (downAtom,upAtom)
        else:
            return False
    #------------------------------------------------------------------------------------------------                      
    def addNewBond(self,connection,k_xy=1.0):
        
        if not self.molecule.getBond(connection[0],connection[1]):
            a,b = self.molecule.atom_stack[connection[0]],self.molecule.atom_stack[connection[1]]
            xs,ys = a.x,a.y
            xf,yf = b.x,b.y
                        
            bond_types = Bond.BOND_TYPES
            b1 = "%s-%s" % (a.sym,b.sym)
            b2 = "%s-%s" % (b.sym,a.sym)

            if b1 in bond_types.keys():
                bond_type = b1
            elif b2 in bond_types:
                bond_type = b2
            else:
                if b1 not in self.ignored_bond_types and b2 not in self.ignored_bond_types and  wx.MessageBox("This bond type (%s,%s) doesn't currently exist. Create a new bond type? (Select No to use default or Yes to create new bond type)" %(a.sym,b.sym),"Create new bond type",style=wx.YES_NO) == wx.YES:
                    dlg = EditBondTypes()
                    create_new = dlg.sketchNew(a.sym,b.sym)
                    if create_new == wx.ID_OK:
                        bond_type = dlg.bond_type
                        Bond.WriteBondLib(dlg.bond_types)
                    else:
                        self.ignored_bond_types.extend([b1,b2])
                        bond_type = "C-C"
                    dlg.Destroy()
                else:
                    self.ignored_bond_types.extend([b1,b2])                    
                    bond_type = "C-C"
                
            artist = self.subplot.plot([xs,xf],[ys,yf],'r',picker=self.PICK_TOLERANCE,zorder=10)[-1]
            
            self.molecule.addBond(artist,connection[0],connection[1],bond_types[bond_type][settings.k_delim])
            self.atom_bond_stack2.append(self.molecule.bond_stack[-1])
            self._redrawFlag = True

        
    #------------------------------------------------------------------------------------------------      
    def mouseUp(self,event):

        self._doRefreshFromHuckel = False
        if event.button == self.DRAW_BUTTON:
            
            self.mouseUpX,self.mouseUpY = event.xdata , event.ydata
            self.mouseDown = False
            
            #just delete the temp line since we will draw a new one if a connection was actually made
            if self.line:
                del(self.subplot.lines[-1])
                self.line = None
                
            #check if a new atom or bond is needed after click
            if event.inaxes:
                self._redrawFlag = True

                if self.newAtomRequired(self.mouseUpX,self.mouseUpY) and not self.moved:


                    hx = Atom.ATOM_TYPES[self.current_atom_type]["h"]
                    self.addNewAtom(self.mouseUpX,self.mouseUpY,hx,self.current_atom_type)

                    #self.solver_update_needed = True
                    
                    self.solver.setData(self.molecule.huckelMatrix())
                else:
                    new_connection = self.connectionMade()
                    if new_connection:
                        self.addNewBond(new_connection)
                        self.solver.setData(self.molecule.huckelMatrix())
                        #self.solver_update_needed = True


            #now change the focus to huckel matrix if user picked an atom or a bond
            if self.set_focus_to_huck :
                huck = self.GetParent().huckel_matrix
                huck.SetGridCursor(self.set_focus_to_huck[0],self.set_focus_to_huck[1])
                #huck.SetFocus()
                huck.SelectBlock(self.set_focus_to_huck[0],self.set_focus_to_huck[1],self.set_focus_to_huck[0],self.set_focus_to_huck[1])
                self.set_focus_to_huck = None

    #------------------------------------------------------------------------------------------------              

    def double_click(self):
        dbl_click = False
        if not self.dbl_click_timer:
            self.dbl_click_timer = time.time()
        elif time.time() - self.dbl_click_timer < settings.DBL_CLICK_TIME:
            dbl_click = True
            self.dbl_click_timer = None
        else:
            self.dbl_click_timer = time.time()
        return dbl_click

    #------------------------------------------------------------------------------------------------              
    def deleteAtomOrBond(self,artist):

        #index of the artist that was picked
        idx = self.subplot.lines.index(artist)

        #delete the artist that was clicked
        self.subplot.lines.pop(idx)
        
        #pop the atom or bond 
        item2 = self.atom_bond_stack2.pop(idx)
        if item2 in self.molecule.atom_stack:
            self.molecule.removeAtom(self.molecule.atom_stack.index(item2))
            pointer = self.GetParent().level_pointer
            if pointer >len(self.molecule.atom_stack)-1:
                self.GetParent().setLevelPointer(0)
            

        elif item2 in self.molecule.bond_stack:
            bond = self.molecule.bond_stack[self.molecule.bond_stack.index(item2)]
            a = self.molecule.atom_stack.index(bond.a)
            b = self.molecule.atom_stack.index(bond.b)
            self.molecule.removeBond(a,b)
            
        self.solver.setData(self.molecule.huckelMatrix())            
#        self.solver_update_needed = True
        self._redrawFlag = True
        
    
    #------------------------------------------------------------------------------------------------              
    def handleAtomPick(self,atom):
        a = self.molecule.atom_stack.index(atom)
        huck_mat = self.GetParent().huckel_matrix
        self.set_focus_to_huck = (a,a)
        
        if self.double_click():
            atom_type = AtomTypeDialog(atom,a)
            result = atom_type.ShowModal()
            if result == wx.ID_OK:
#                atom.setData(atom.x,atom.y,-1.*abs(float(atom_type.hx.GetValue())),atom_type.atom_type.GetValue())
                atom.setData(atom.x,atom.y,float(atom_type.hx.GetValue()),atom_type.atom_type.GetValue())
                bonds = [bond for bond in self.molecule.bond_stack if bond.a == atom or bond.b == atom]
                
                for bond in bonds:
                    update = True
                    b1 = Bond.FMT % (bond.a.sym,bond.b.sym)
                    b2 = Bond.FMT % (bond.b.sym,bond.a.sym)
                    if b1 not in Bond.BOND_TYPES.keys() and b2 not in Bond.BOND_TYPES.keys():
                        if wx.MessageBox("This bond type (%s) doesn't currently exist. Create a new bond type? (Select No to use default or Yes to create new bond type)" %(b1),"Create new bond type",style=wx.YES_NO) == wx.YES:
                            dlg = EditBondTypes()
                            create_new = dlg.sketchNew(bond.a.sym,bond.b.sym)
                            if create_new == wx.ID_OK:
                                bond_type = dlg.bond_type
                                Bond.WriteBondLib(dlg.bond_types)
                            dlg.Destroy()       
                        else:
                            update = False
                    if update:
                        if b1 in Bond.BOND_TYPES.keys():
                            bond.k_xy = Bond.BOND_TYPES[b1][settings.k_delim]
                        else:
                            bond.k_xy = Bond.BOND_TYPES[b2][settings.k_delim]
                self._redrawFlag = True      
                self.solver.setData(self.molecule.huckelMatrix())
#                self.solver_update_needed = True
                
            atom_type.Destroy()
            self.back_from_dialog = True

    #------------------------------------------------------------------------------------------------              
    def handleBondPick(self,bond):
        a,b = self.molecule.atom_stack.index(bond.a),self.molecule.atom_stack.index(bond.b)

        #rather than setting focus now we will do it later in mouseUp since focus will be restolen by the sketch pad if we do it here 
        self.set_focus_to_huck = (a,b)

        if self.double_click():
            bond_type = BondTypeDialog(bond,"Set k_xy for bond between atoms %d && %d (%s-%s):" % (a+1,b+1,bond.a.sym,bond.b.sym))
            result = bond_type.ShowModal()
            if result == wx.ID_OK:
                
#                bond.k_xy = -1.*abs(float(bond_type.k_xy.GetValue()))
                bond.k_xy = float(bond_type.k_xy.GetValue())
                if abs(bond.k_xy) < settings.eps:
                    self.deleteAtomOrBond(bond.artist)
                #self.solver_update_needed = True
                self.solver.setData(self.molecule.huckelMatrix())
                #self._redrawFlag = True
                
            bond_type.Destroy()
            self.back_from_dialog = True
        
            
    #------------------------------------------------------------------------------------------------              
    def pickEvent(self,event):
        mouse_event = event.mouseevent
        
        #index of the artist that was picked
        idx = self.subplot.lines.index(event.artist)
        
        if mouse_event.button == self.DEL_BUTTON:
            self.deleteAtomOrBond(event.artist)
        
        if mouse_event.button == self.DRAW_BUTTON:
            
            item = self.atom_bond_stack2[idx]
            
            if item in self.molecule.atom_stack:
                self.handleAtomPick(item)

            else:
                if self.set_focus_to_huck == None:
                    #above condition is to ensure that we select the atom when multiple atoms/bonds are chosen
                    self.handleBondPick(item)
                    
    #------------------------------------------------------------------------------------------------      
    def mouseMove(self,event):
        

        if event.xdata != None and event.ydata != None:
            
            if self.mouseDown:
                r2 = (self.mouseDownX-event.xdata)*(self.mouseDownX-event.xdata)+(self.mouseDownY-event.ydata)*(self.mouseDownY-event.ydata)
                if r2 > self.XTOL*self.YTOL:
                    self.moved = True
            
            if len(self.atom_bond_stack2) == 0:
                self.tooltip.SetTip('Left-Click to create an atom')
                self.tooltip.Enable(True)
                
            elif self.onAtom(event.xdata,event.ydata) >=0 and self.line == None:
                tip = 'Left-Click and hold on atom then drag to another atom to create a bond\n' \
                    'Ctrl+Left-Click and hold on atom to move\n'\
                    'Double-Click on atom or bond to set hx or k_xy\n'\
                    'Right-Click to delete atom or bond\n'\
                    'Scroll-Wheel to zoom\n'\
                    'Ctrl+Scroll-Wheel to rotate plot'
                self.tooltip.SetTip(tip)
                self.tooltip.Enable(True)
            else:
                self.tooltip.Enable(False)

            
        if self.mouseDown and self.mouse_down_on_atom >=0 and event.button == self.DRAW_BUTTON:
            if not event.key == "control":
    
                self.lineXs = [self.mouseDownX,event.xdata]
                self.lineYs = [self.mouseDownY,event.ydata]
    
                if not self.line:
                    self.line = self.subplot.plot(self.lineXs,self.lineYs,'-r')[-1]
                else:
                    self.line.set_data(self.lineXs,self.lineYs)

            else:
                atom = self.molecule.atom_stack[self.mouse_down_on_atom]
                atom.setPos(event.xdata,event.ydata)
                for bond in self.molecule.bond_stack:
                    if atom in bond.atoms:
                        bond.refresh()
#                self.draw
            self.subplot.axis(self.AXIS_MIN_MAX)
            self.canvas.draw()
    #------------------------------------------------------------------------------------------------          
    def mouseWheel(self,event):

        if event.key == "control":
            if event.step < 0:
                self.rotate(-1,event.xdata,event.ydata)
            else:
                self.rotate(1,event.xdata,event.ydata)
            
        else:
            if event.step < 0:
                self.resize(-1,event.xdata,event.ydata)
            else:
                self.resize(1,event.xdata,event.ydata)
    #------------------------------------------------------------------------------------------------          
    def rotate(self,dir,xo=0,yo=0):
        axis = self.subplot.axis()
        xo =  (axis[1]+axis[0])/2.
        yo =  (axis[3]+axis[2])/2.
        x = numpy.arange(0.,2*numpy.pi,0.05)
        trans = matplotlib.transforms.blended_transform_factory(self.subplot.transData,self.subplot.transAxes)
        for atom in self.molecule.atom_stack:
            pos = numpy.array([atom.x-xo,atom.y-yo]).reshape((2,1))

            dt = 5*math.pi/180.*dir
            mat = numpy.mat([[math.cos(dt),-math.sin(dt)],[math.sin(dt),math.cos(dt)]],float)
            new_pos = mat*pos
            x_pos = new_pos[0,0]+xo
            y_pos = new_pos[1,0]+yo


            atom.setPos(x_pos,y_pos)

        for bond in self.molecule.bond_stack:
            bond.refresh()
            
        self.draw()        
        self.GetParent().results_display_2dmo._redrawFlag = True
        
    def _onIdle( self, evt ):

        PlotPanel._onIdle(self,evt)
        if self.solver_update_needed:
            
            
#            self.solver.setData(self.molecule.huckelMatrix())
            
            
            #self.solver.setNumElectrons(self.molecule.numAtoms)
            
            
            #self.GetParent().setLevelPointer(0)
            
            
            self.solver_update_needed = False
            
        if self._redrawFlag:
            
            self._redrawFlag = False
            self.draw()
            
        
    def _getRandPos(self,n):
        return [(math.cos(x*2*math.pi/n),math.sin(x*2*math.pi/n)) for x in range(0,n)]

    def createFromHuckel(self,data,minimize=True):

        if 0 not in data.shape and data.shape[1] == data.shape[0]:
            size = data.shape[0]

            init_pos = self._getRandPos(size)
                
            for ii in range(size):
                self.addNewAtom(init_pos[ii][0],init_pos[ii][1],data[ii,ii],sym = 'A')
                for jj in range(ii):
                    if abs(data[ii,jj])>settings.eps:
                        self.addNewBond((ii,jj),data[ii,jj])
            if minimize:
                self.molecule.minimizePositions()
                        
        
    #------------------------------------------------------------------------------------------------          
    def refreshFromHuckel(self):
#        if self._doRefreshFromHuckel:

        if self.solver.getSize()>0 and self.GetParent().visual_mode.IsChecked():
            huck = self.solver.data
            cur_huck = self.molecule.huckelMatrix()

            for ii in range(self.molecule.numAtoms):
                self.molecule.atom_stack[ii].hx = huck[ii,ii]
                for jj in range(ii):
                    if abs(cur_huck[ii,jj])<settings.eps and abs(huck[ii,jj])>settings.eps:
                        #new bond was created in the huckel matrix input
                        self.addNewBond((ii,jj))
                        
                    elif abs(cur_huck[ii,jj])>settings.eps and abs(huck[ii,jj])<settings.eps:
                        #bond was deleted from huckel matrix input
                        self.deleteAtomOrBond(self.molecule.getBond(ii,jj).artist)
    
                    if abs(huck[ii,jj])>settings.eps:
                        self.molecule.getBond(ii,jj).k_xy = huck[ii,jj]
            #self.molecule.updateFromHuckel(huck)
            self._redrawFlag = True
 #   self._doRefreshFromHuckel = True
    
    #------------------------------------------------------------------------------------------------              
    def initializePlot(self):

        self.subplot = self.figure.add_subplot( 111 )
        self.figure.suptitle('Molecule Sketch Pad')
        self.figure.subplots_adjust(left=0.05,right=0.95,top=0.925,bottom=0.05)
        self.canvas.SetToolTip(self.tooltip)
        self.canvas.mpl_connect('button_press_event', self.mouseDown)
        self.canvas.mpl_connect('button_release_event', self.mouseUp)
        self.canvas.mpl_connect('motion_notify_event', self.mouseMove)
        self.canvas.mpl_connect('scroll_event',self.mouseWheel)
        self.canvas.mpl_connect('pick_event',self.pickEvent)
        self.subplot.axis(self.AXIS_MIN_MAX)
        self.subplot.set_aspect('auto',adjustable="datalim")
        self.subplot.set_xticklabels([""])
        self.subplot.set_yticklabels([""])

        
    #------------------------------------------------------------------------------------------------              
    def draw( self ):
        if not hasattr( self, 'subplot' ):
            self.initializePlot()            

        if len(self.atom_bond_stack2)>=0:
            self.relable_atoms()
        self.subplot.axis('equal')            
        self.subplot.axis(self.AXIS_MIN_MAX)
        self.canvas.draw()
#------------------------------------------------------------------------------------------------                      
    def resize(self,zoom = 0,xo=0,yo=0):
        if len(self.molecule.atom_stack)>1:
            if zoom == 0 :
                # auto zoom out 
                scale = 1
                xs = [atom.x for atom in self.molecule.atom_stack]
                ys = [atom.y for atom in self.molecule.atom_stack]
                min_x,max_x = min(xs),max(xs)
                min_y,max_y = min(ys),max(ys)
                
                dx = abs(max_x-min_x)*scale
                dy = abs(max_y-min_y)*scale
                ax = [min_x-dx,max_x+dx,min_y-dy,max_y+dy]
    
            else:
                #zoom in or out by scale %
                scale = 0.1
                ax = self.subplot.axis()
                dx = abs(ax[0]-ax[1])*scale
                dy = abs(ax[2]-ax[3])*scale
                if zoom < 0:
                    ax = [ax[0]-dx,ax[1]+dx,ax[2]-dy,ax[3]+dy]
                else:
                    ax = [ax[0]+dx,ax[1]-dx,ax[2]+dy,ax[3]-dy]

            self.XTOL = abs(ax[1]-ax[0])*self.TOL_RATIO
            self.YTOL = self.XTOL
            
            self.subplot.axis('equal')
            self.subplot.axis(ax)
            self.AXIS_MIN_MAX = ax
            self._redrawFlag = True
        
#------------------------------------------------------------------------------------------------                      
    def relable_atoms(self):
        xs = []
        ys = []
        hxs = []
        syms = []
        
        for atom in self.molecule.atom_stack:
            xs.append(atom.x)
            ys.append(atom.y)
            hxs.append(atom.hx)
            syms.append(atom.sym)
        
        self.subplot.texts = []
        
        delta = 0.1
        for ii in range(len(xs)):
            x,y = xs[ii],ys[ii]
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
            self.subplot.text(x+dx,y+dy,"%s $_{i=%d}$" % (syms[ii],ii+1),ha=ha,va=va,size='large',zorder=20,color='black')              
