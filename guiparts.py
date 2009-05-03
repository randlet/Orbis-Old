import time
import wx
import settings
#import matplotlib
#matplotlib.use( 'WXAgg',warn=False )

#matplotlib.interactive( False )
#from simplehuckel import matplotlib,FigureCanvasWxAgg,Figure,pylab
#from matplotlib.figure import Figure
#import pylab
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
import matplotlib
#import pylab

from molecule import Atom
import numpy
import math
import molecule
from scipy.optimize import fmin_l_bfgs_b
from TimedToolTip import TimedToolTip
ID_SPIN_CTRL_NUME = 200

class PlotPanel (wx.Panel):
    """The PlotPanel has a Figure and a Canvas. OnSize events simply set a flag, and the actual resizing of the figure is triggered by an Idle event."""

    def __init__( self, parent, color=None, dpi=None, **kwargs ):

        # initialize Panel
        if 'id' not in kwargs.keys():
            kwargs['id'] = wx.ID_ANY
        if 'style' not in kwargs.keys():
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, parent, **kwargs )
        

        # initialize matplotlib stuff
        self.figure = Figure( None, dpi )

        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )

        self._SetSize()
        self.draw()
        self.SetColor()
        self._resizeflag = False

        self.Bind(wx.EVT_IDLE, self._onIdle)
        self.Bind(wx.EVT_SIZE, self._onSize)

    def SetColor( self, rgbtuple=None ):
        """Set figure and canvas colours to be the same."""
        if rgbtuple is None:
            rgbtuple = wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNFACE ).Get()
        clr = [c/255. for c in rgbtuple]
        self.figure.set_facecolor( clr )
        self.figure.set_edgecolor( clr )
        self.canvas.SetBackgroundColour( wx.Colour( *rgbtuple ) )

    def _onSize( self, event ):
        self._resizeflag = True
        event.Skip()
    def _onIdle( self, evt ):
        
        if self._resizeflag:
            self._resizeflag = False
            self._SetSize()
        evt.Skip()
    def _SetSize( self ):
        pixels = self.GetSize()
        self.SetSize( pixels )
        self.canvas.SetSize( pixels )
        self.figure.set_size_inches( float( pixels[0] )/self.figure.get_dpi(),
                                     float( pixels[1] )/self.figure.get_dpi() )

    def draw(self): pass # abstract, to be overridden by child classes


class ELDPlotPanel (PlotPanel):
    """Plots several lines in distinct colors."""
    PICK_TOLERANCE = 5

    def __init__( self, parent,**kwargs ):

        self.tooltip = TimedToolTip(settings.tool_tip_time,tip='tip with a long %s line and a newline\n' % (' '*100))
        self.tooltip.SetTip('Left-Click on energy level to see corresponding orbital diagram and eigenvector')
        self.tooltip.Enable(False)
        self._redrawFlag = True
        self.tooltip.SetDelay(500)
        self.show_tip = True
        self.solver = parent.huckel_solver
        self.levels = self.solver.populated_levels
        # initiate plotter
        PlotPanel.__init__( self, parent, **kwargs )
        self.Bind(wx.EVT_IDLE,self.tooltip._onIdle)
        
        
    def refreshFromHuckel(self):
        self.levels = self.solver.populated_levels
        self._redrawFlag = False
        self.draw()

    def pickEvent(self,event):

        #index of the artist that was picked
        try:
            idx = self.subplot.lines.index(event.artist)
            
            self.GetParent().setLevelPointer(idx)
        except:
            pass

    def mouseMove(self,event):
        
        if event.xdata != None and event.ydata != None and len(self.levels)>0 and self.show_tip == True:
            self.tooltip.SetTip('Left-Click on energy level to see corresponding orbital diagram and eigenvector')
            self.tooltip.Enable(True)
            self.show_tip = False
            
        #else:
        #    self.tooltip.Enable(False)
            
    def axesEnter(self,event):

        self.show_tip = True
        
    def _onIdle(self,evt):
        
        PlotPanel._onIdle(self,evt)
        if self._redrawFlag:
            
            self.draw()
            self._redrawFlag = False
        evt.Skip()
#            self.canvas.draw()
        
    def getDegenLevels(self):
        #if energies are within this % then they are considered splitting of a single level
        level_delta = 0.001
        
        last_energy = self.levels[0][0]

        degen_levels =[[self.levels[0]]]
        degen_level_idx = 0
        
        for ii,level in enumerate(self.levels[1:]):
            energy,vec,ne = level
            
            if abs(last_energy-energy)<settings.eps or abs(1.-energy/last_energy)<level_delta:
                degen_levels[-1].append(level)
            else:
                degen_levels.append([level])

            last_energy = energy
            
        return degen_levels
    
    def draw( self ):

        if not hasattr( self, 'subplot' ):
            self.subplot = self.figure.add_subplot( 111 )
            self.figure.suptitle('Energy Level Diagram')
            self.canvas.SetToolTip(self.tooltip)
            self.canvas.mpl_connect('pick_event',self.pickEvent)
            self.canvas.mpl_connect('motion_notify_event', self.mouseMove)            
            self.canvas.mpl_connect('axes_enter_event',self.axesEnter)
            self.subplot.set_ylabel("Energy")
            self.subplot.set_xticklabels([""])            
            self.subplot.set_xticks([])            
            self.drawLegend()

            self.width = 1.
            self.space = 0.5*self.width
            self.width_space = self.width+self.space
            self.steps = 4
            self.step = self.width/self.steps
            self.height_fac = 0.2
            
        else:
            for artist in self.subplot.get_children():
                if isinstance(artist,matplotlib.lines.Line2D):
                    artist.remove()


        if len(self.levels)>0:
            
            de = (self.levels[-1][0]-self.levels[0][0])*self.height_fac 
            min_energy = self.levels[0][0] - de
            max_energy = self.levels[-1][0] + de

            if abs(max_energy-min_energy)<settings.eps:
                min_energy -= 1.
                max_energy += 1.

            max_widths = [8]
            level_idx = 0
            level_pointer = self.GetParent().level_pointer
            degen_levels = self.getDegenLevels()           

            for levels in degen_levels:                    

                energy= levels[0][0]
                
                n = len(levels)
                total_width = n*(self.width_space)-self.space

                max_widths.append(total_width)

                for ii in range(n):
                    s = -total_width*0.5 + ii*(self.width_space)
                    f = s+self.width
                    ne = levels[ii][2]
                    
                    x = numpy.arange(s,f+self.step,self.step)
                    y = [energy]*len(x)

                    fmt = self.getFmt(ne)
                    
                    if level_idx == level_pointer:
                        fmt += 'o'
                        markersize = 8
                    else:
                        markersize = 1

                    self.subplot.plot( [x[0],x[-1]], [y[0],y[-1]],fmt,picker=self.PICK_TOLERANCE,linewidth = 1,markersize=markersize)

                    level_idx += 1
            
            max_width = max(max_widths)
            self.subplot.axis([-0.6*max_width,0.6*max_width,min_energy,max_energy])

        self.canvas.draw()

    def getFmt(self,ne=-1):
        fmts = (('.w','Filling:'),('k-','N=0'),('b--','0<N<1'),('b-','N=1'),('r--','1<N<2'),('r-','N=2'))                                        
        if ne <0:
            return fmts
        elif ne <settings.eps:
            return fmts[1][0]
        elif 0<ne < 1:
            return fmts[2][0]
        elif abs(ne-1)<settings.eps:
            return fmts[3][0]
        elif 1 < ne <2:
            return fmts[4][0]
        else :
            return fmts[5][0]

 

        
    
    def drawLegend(self):


        [self.subplot.plot([-999],[0],fmt[0],label=fmt[1],linewidth=3) for fmt in self.getFmt()]
            
        font = matplotlib.font_manager.FontProperties(size=10)
        legend = self.subplot.legend(loc=8,ncol=3,prop=font,columnspacing=1,markerscale=4)
        legend.draw_frame(False)
        


class HuckelMatrix(wx.grid.Grid):
    #TODO: consolidate HuckelMatrix and Results Matrix

    COL_WIDTH = 36
    INIT_BASIS_SIZE = 0
    COPY_DELIM = '\t'

    def __init__(self, parent, ID=-1, label="", pos=wx.DefaultPosition, size=(100, 25)): 

        wx.grid.Grid.__init__(self,parent,ID,pos,size,wx.RAISED_BORDER,label)

        self.Bind(wx.grid.EVT_GRID_EDITOR_CREATED, self.OnGridEditorCreated)
        self.Bind(wx.EVT_MENU_OPEN, self.OnMenuOpen)
        self.Bind(wx.EVT_KEY_UP,self.onKeyPress)
        self.solver = self.GetParent().huckel_solver
        

        self.tooltip = TimedToolTip(settings.tool_tip_time,tip='Click on a row to see corresponding orbital' )
        #self.tooltip.Enable(False)
        self.tooltip.SetDelay(500)
        #self.show_tip = True
        
        self.SetToolTip(self.tooltip)
        
        num_rows = self.solver.getSize()
        num_cols = num_rows
        self.CreateGrid(num_cols,num_rows)

        self.setLabels()

        self.SetMargins(5,5)
        self.SetDefaultColSize(self.COL_WIDTH)
        self.SetRowLabelSize(self.COL_WIDTH)

        init_data = self.solver.data
        
        self.setData(init_data)

    def onKeyPress(self,event):

        if event.ControlDown() and event.GetKeyCode() in [67,322]:
            self.copy()

        event.Skip()

    def copy(self):
        data = self.getData()
        paster = wx.TextDataObject()

        paste_data = ''
        num = self.GetNumberCols()

        for ii in range(num):
            paste_data += self.COPY_DELIM.join([str(data[ii,jj]) for jj in range(num)]+['\n'])
        paster.SetText(paste_data)

        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(paster)
            wx.TheClipboard.Close()

    def refreshFromHuckel(self):

        self.setSize(self.solver.getSize())
        
        self.setData(self.solver.data)

    def setData(self,data):
        assert type(data) == numpy.matrix
        for ii in range(data.shape[0]):
            for jj in range(ii+1):
                datum = data[ii,jj]
                if abs(datum) < settings.eps:
                    if ii == jj:
                        datum = "0.0"
                    else:
                        datum = ""
                self.SetCellValue(ii,jj,str(datum))
                self.SetCellValue(jj,ii,str(datum))

    def getData(self):
        num = self.GetNumberCols()
        data = numpy.mat(numpy.zeros((num,num),float))
        for ii in range(num):
            for jj in range(num):
                val = self.GetCellValue(ii,jj)
                if val =="":
                    val = '0'
                data[ii,jj] = float(val)
        return data

    def setSize(self,size):

        cur_size = self.GetNumberCols()

        diff = size - cur_size

        if diff > 0:
            self.AppendCols(diff,updateLabels=False)
            self.AppendRows(diff,updateLabels=False)
            for ii in range(cur_size,size):
                for jj in range(ii+1):
                    self.SetCellValue(ii,jj,"0.0")
                    self.SetCellValue(jj,ii,"0.0")

        if diff < 0 and size>-1:
            diff = abs(diff)
            self.DeleteCols(cur_size-diff,diff,updateLabels=False)
            self.DeleteRows(cur_size-diff,diff,updateLabels=False)

        self.setLabels()        

    def setLabels(self,custom=[]):
        size = self.GetNumberCols()
        for ii in range(size):
            if custom and len(custom)>=size:
                label = custom[ii]
            else:
                label = str(ii+1)

            self.SetColLabelValue(ii,label)


    def OnGridEditorCreated(self, event):
        """ Bind the kill focus event to the newly instantiated cell editor """
        editor = event.GetControl()
        editor.Bind(wx.EVT_KILL_FOCUS, self.OnKillFocus)


    def OnKillFocus(self, event):

        row,col = self.GetGridCursorRow(),self.GetGridCursorCol()
        val = self.GetCellValue(row,col)

        try :
            float(val)
        except:
            val = "0.0"

        self.SetCellValue(col,row,val) 


        self.SaveEditControlValue()
        self.HideCellEditControl()

        self.solver.setMatrixElement(row,col,float(val))
        self.solver.setMatrixElement(col,row,float(val))

    def OnMenuOpen(self, event):
        self.HideCellEditControl()    


    def OnCut(self, event): # Cut selection
        """ Cuts the selection """
        self.grid.Cut()
        self.PositionUpdate()

    def OnCopy(self, event): # Copy Selection
        """ Copies the selection """
        self.grid.Copy()
        self.PositionUpdate()
        event.skip()

    def OnPaste(self, event): # Paste Selection
        """ Paste the Cut or Copied elements from the clipboard """
        self.grid.Paste()
        self.PositionUpdate()

    def OnDelete(self, event): # Delete the selection
        """ Deletes the selected portion """
        fromHere, toHere = self.control.GetSelection()
        self.control.Remove(fromHere, toHere)	


class ResultsMatrix(wx.grid.Grid):

    COL_WIDTH = 40
    FMT = "%.4G"
    COPY_DELIM = '\t'

    def __init__(self, parent, ID=-1, label="", pos=wx.DefaultPosition, size=(-1, -1),row_labels=[],col_labels=[]): 

        wx.grid.Grid.__init__(self,parent,ID,pos,size,wx.RAISED_BORDER,label)
        self.solver = self.GetParent().GetGrandParent().huckel_solver

        num_rows = self.solver.getSize()

        num_cols = num_rows
        self.CreateGrid(num_cols,num_rows)

        self.setLabels(row_labels,col_labels)

        self.SetMargins(5,5)
        self.SetDefaultColSize(self.COL_WIDTH)
        #self.SetRowLabelSize(60)

        data = numpy.mat(numpy.zeros((num_rows,num_rows),float))

        self.setData(data)

        self.Bind(wx.EVT_KEY_UP,self.onKeyPress)

    def onKeyPress(self,event):

        if event.ControlDown() and event.GetKeyCode() in [67,322]:
            self.copy()

        event.Skip()

    def copy(self):
        data = self.getData()
        paster = wx.TextDataObject()

        paste_data = ''
        num_col = self.GetNumberCols()
        num_row = self.GetNumberRows()

        for ii in range(num_row):
            paste_data += self.COPY_DELIM.join([str(data[ii,jj]) for jj in range(num_col)]+['\n'])
        paster.SetText(paste_data)

        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(paster)
            wx.TheClipboard.Close()

    def getData(self):
        num_row = self.GetNumberRows()
        num_col = self.GetNumberCols()
        data = numpy.mat(numpy.zeros((num_row,num_col),float))
        for ii in range(num_row):
            for jj in range(num_col):
                val = self.GetCellValue(ii,jj)
                if val =="":
                    val = '0'
                data[ii,jj] = float(val)
        return data



    def setLabels(self,row_labels=[],col_labels=[],reverse=False):
        num_rows = self.GetNumberRows()
        num_cols = self.GetNumberCols()
        for ii in range(num_cols):
            if len(col_labels) >= num_cols:
                label = col_labels[ii]
            else:
                label = str(ii+1)
            self.SetColLabelValue(ii,label)

        for ii in range(num_rows):
            
            if len(row_labels) >= num_rows:
                label = row_labels[ii]
            else:
                label = str(ii+1)

            if reverse:
                self.SetRowLabelValue(num_rows-ii-1,label)
            else:
                self.SetRowLabelValue(ii,label)
            
            self.AutoSizeRowLabelSize(ii)
            
        self.ForceRefresh()


    def setData(self,data,reverse=False):
        
        size_ii,size_jj = data.shape
        
        
        for ii in range(size_ii):
            for jj in range(size_jj):
                val = data[ii,jj]
                if abs(val)<settings.eps:
                    val = 0.0

                if reverse:
                    self.SetCellValue(size_ii-ii-1,jj, self.FMT %(val))
                    self.SetReadOnly(size_ii-ii-1,jj,True)
                else:
                    self.SetCellValue(ii,jj, self.FMT %(val))
                    self.SetReadOnly(ii,jj,True)

        self.AutoSizeRows()
        self.AutoSizeColumns()
        
    def setSize(self,row_size,col_size=-1):
        
        if col_size<0:
            col_size = row_size
            
        cur_size = self.GetNumberCols()

        diff = col_size - cur_size
        if diff > 0:
            self.AppendCols(diff,updateLabels=False)

        if diff < 0 and col_size>-1:
            diff = abs(diff)
            self.DeleteCols(cur_size-diff,diff,updateLabels=False)

        cur_size = self.GetNumberRows()
        diff = row_size - cur_size

        if diff > 0:
            self.AppendRows(diff,updateLabels=False)
        if diff < 0 and row_size>-1:
            diff = abs(diff)
            self.DeleteRows(cur_size-diff,diff,updateLabels=False)
        
    def refreshFromHuckel(self):
        pass


class AtomAtomPolarizabilityMatrix(ResultsMatrix):

    def __init__(self, parent, ID=-1, label="", pos=wx.DefaultPosition, size=(100, 25),row_labels=[],col_labels=[]): 
        ResultsMatrix.__init__(self, parent, ID=ID, label=label, pos=pos, size=size,row_labels=row_labels,col_labels=col_labels)

    def refreshFromHuckel(self):
        self.setSize(self.solver.getSize())
        data = self.solver._calcAAPolarizability() #aa_polar

        if len(data)>0:
            self.setLabels()
            self.setData(data)


class NetChargeMatrix(ResultsMatrix):

    def __init__(self, parent, data_name,ID=-1, label="", pos=wx.DefaultPosition, size=(100, 25),row_labels=[],col_labels=[]): 
        self.data_name = data_name
        ResultsMatrix.__init__(self, parent, ID=ID, label=label, pos=pos, size=size,row_labels=row_labels,col_labels=col_labels)

        
    def refreshFromHuckel(self):

        data = self.solver._calcNetCharges()#net_charges
        na = self.solver.getSize()
        if na >= 1:
            self.setSize(1,na)
        else:
            self.setSize(0)

        if len(data)>0:
            disp_data = numpy.mat(numpy.zeros((1,na),float))
            
            col_labels = ["Atom\n%d" % (x+1)for x in range(len(data))]
            row_labels = ["Net Charge"]
            
            for ii in range(na):
                val = data[ii]
                disp_data[0,ii] = val 
                    
            self.setLabels(row_labels,col_labels)

            self.setData(disp_data)

            
class AtomBondPolarizabilityMatrix(ResultsMatrix):

    def __init__(self, parent, data_name,ID=-1, label="", pos=wx.DefaultPosition, size=(100, 25),row_labels=[],col_labels=[]): 
        self.data_name = data_name
        ResultsMatrix.__init__(self, parent, ID=ID, label=label, pos=pos, size=size,row_labels=row_labels,col_labels=col_labels)

    def createData(self):
    
        #calculate atom bond polarizability when required
        data = self.solver._calcABPolarizability()
#        data = self.solver.ab_polar
        na,nb = self.solver.getSize(),self.solver.getNumBonds()
        disp_data = numpy.mat(numpy.zeros((na,nb),float))
        if len(data)>0:
            
            col_labels = ["Bond\n(%d,%d)" % (x+1,y+1)for x,y,z in data[0]]
            row_labels = ["Atom %d" % (x+1)for x in range(len(data))]
            
            for ii in range(na):
                for jj in range(nb):
                    val = data[ii][jj][2]
                    disp_data[ii,jj] = val
                    
            self.setLabels(row_labels,col_labels)
        
            self.setData(disp_data)
        
    def refreshFromHuckel(self):
        self.setSize(self.solver.getSize(),self.solver.getNumBonds())
        self.createData()



class EigenMatrix(ResultsMatrix):

    def __init__(self, parent, ID=-1, label="", pos=wx.DefaultPosition, size=(-1, -1),row_labels=[],col_labels=[]): 

        ResultsMatrix.__init__(self, parent, ID=ID, label=label, pos=pos, size=size,row_labels=row_labels,col_labels=col_labels)

        self.Bind(wx.grid.EVT_GRID_SELECT_CELL, self.OnGridCellClicked)
        self.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK,self.OnGridCellClicked)

    def OnGridCellClicked(self,event):
        
        
        row = event.GetRow()

        if row >= 0:
            self.SelectRow(row)
            self.GetParent().GetGrandParent().setLevelPointer(self.solver.getSize()-row-1)

        #only want to propogate cursor events and not when col labels are clicked.
        if event.GetEventType() == wx.grid.EVT_GRID_SELECT_CELL.evtType[0]:
            event.Skip()

    def refreshFromHuckel(self):
        self.setSize(self.solver.getSize())
        vecs = self.solver.eigen_vecs
        count = len(vecs)
        if count>0:
            row_labels = ["E = " +x for x in map(lambda x : "%5.3f" % (x),self.solver.eigen_vals)]
            self.setLabels(row_labels=row_labels,reverse=True)
            
            self.setData(numpy.mat(vecs),reverse=True)
            lp = self.GetParent().GetGrandParent().level_pointer
            if len(vecs)>lp:
                self.SelectRow(count-lp-1)


class PiBondMatrix(ResultsMatrix):

    def __init__(self, parent, ID=-1, label="", pos=wx.DefaultPosition, size=(100, 25),row_labels=[],col_labels=[]): 

        ResultsMatrix.__init__(self, parent, ID=ID, label=label, pos=pos, size=size,row_labels=row_labels,col_labels=col_labels)

    def refreshFromHuckel(self):

        self.setSize(self.solver.getSize())
#        data = self.solver.bond_orders
        data = self.solver._calcBondOrders()
#        if data.any():
        self.setLabels()
        self.setData(data)


class ControlPanel(wx.Panel): 
    def __init__(self, parent, ID=-1, label="", pos=wx.DefaultPosition, size=(200,100)): 
        wx.Panel.__init__(self, parent, ID, pos, size, wx.RAISED_BORDER|wx.EXPAND, label) 

        self.label = label 
        sizer = wx.BoxSizer(wx.VERTICAL)

        self.solver = self.GetParent().huckel_solver
        self.sketch_pad = self.GetParent().sketch_pad
        self.eigen_plot = self.GetParent().results_display_2dmo
        self.huckel_matrix = self.GetParent().huckel_matrix
        size = self.solver.getSize()


        #set up number of electrons control
        self.num_e = wx.SpinCtrl(self,-1,"",min=settings.MIN_NUM_E,max=size*2,initial=size,name="num_e")
        num_e_box = wx.StaticBox(self,label="Number of electrons",style=wx.EXPAND)
        num_e_sizer = wx.StaticBoxSizer(num_e_box,wx.VERTICAL)
        num_e_sizer.Add(self.num_e,flag=wx.EXPAND)
        self.num_e.Bind(wx.EVT_KEY_DOWN,self.onKeyPressNumE)
        self.Bind(wx.EVT_SPINCTRL,self.onNumE,self.num_e)
        sizer.Add(num_e_sizer,flag=wx.EXPAND,border=20)

        #set up basis set size control
        self.basis_size = wx.SpinCtrl(self,-1,"",min=0,max=40,initial=size,name="basis_size")
        self.basis_size.Enable(False)
        basis_size_box = wx.StaticBox(self,label="Size of basis set")
        basis_size_sizer = wx.StaticBoxSizer(basis_size_box,wx.VERTICAL)
        basis_size_sizer.Add(self.basis_size,flag=wx.EXPAND)
        self.basis_size.Bind(wx.EVT_KEY_DOWN,self.onKeyPressBasis)
        self.Bind(wx.EVT_SPINCTRL,self.onBasisSize,self.basis_size)
        sizer.Add(basis_size_sizer,flag=wx.EXPAND)

        atom_types = Atom.ATOM_TYPES
        
        atom_list = ["%s - %s" % (x,atom_types[x]["description"]) for x in sorted(atom_types.keys())]
        if "C" in atom_types.keys():
            init = atom_list[sorted(atom_types.keys()).index("C")]
        else:
            init = atom_list[0]
        self.atom_type = wx.ComboBox(self,-1,init,wx.DefaultPosition,wx.DefaultSize,atom_list,wx.CB_DROPDOWN|wx.CB_READONLY)
        atom_box = wx.StaticBox(self,label="Atom Type",style=wx.EXPAND)
        atom_sizer = wx.StaticBoxSizer(atom_box,wx.VERTICAL)
        atom_sizer.Add(self.atom_type,wx.EXPAND)
        self.atom_type.Bind(wx.EVT_COMBOBOX,self.onAtomType)
        
        sizer.Add(atom_sizer)
        
        self.minimize = wx.Button(self,-1,"Redraw",size=(125,-1))
        self.minimize.SetToolTip(wx.ToolTip("Redraw the current molecule"))
        self.Bind(wx.EVT_BUTTON,self.onMinimize,self.minimize)
        sizer.Add(self.minimize)

        self.clear = wx.Button(self,-1,"Clear",size=(125,-1))
        self.clear.SetToolTip(wx.ToolTip("Clear the current session"))
        self.Bind(wx.EVT_BUTTON,self.onClear,self.clear)
        sizer.Add(self.clear)


        self.SetSizer(sizer)
        self.Layout()

    def onAtomType(self,event):
        atype = self.atom_type.GetValue().split('-')[0].strip()
        self.GetParent().sketch_pad.current_atom_type = atype
        
    def onClear(self,event):

        self.solver.reset()
        if self.GetParent().visual_mode.IsChecked():
            self.sketch_pad.reset()

    def onMinimize(self,event):
        if len(self.sketch_pad.molecule.atom_stack)>0:
            self.sketch_pad.molecule.minimizePositions()
            self.sketch_pad.resize()
            self.eigen_plot.draw()
            self.sketch_pad.SetFocus()

    def onKeyPressBasis(self,event):
        key = event.GetKeyCode()

        if key == wx.WXK_RETURN or key == wx.WXK_NUMPAD_ENTER:
            self.onBasisSize(event)
        event.Skip()

    def onKeyPressNumE(self,event):
        key = event.GetKeyCode()

        if key == wx.WXK_RETURN or key == wx.WXK_NUMPAD_ENTER:
            self.onNumE(event)
        event.Skip()

    def onBasisSize(self,event):
        size = self.basis_size.GetValue()
        self.huckel_matrix.setSize(size)
        data = self.huckel_matrix.getData()
        #data = numpy.matrix(numpy.zeros((size,size),float))
        
        self.solver.setData(data)
        #self.solver.setNumElectrons(size)
        #self.num_e.SetValue(size)
        
        event.Skip()

    def onNumE(self,event):
        self.solver.setNumElectrons(self.num_e.GetValue())
        event.Skip()
        
        
    def refreshFromHuckel(self):

        size = self.solver.getSize()
        self.basis_size.SetValue(size)
        num_e = self.solver.num_e
        self.num_e.SetRange(settings.MIN_NUM_E,size*2)
        self.num_e.SetValue(num_e)







