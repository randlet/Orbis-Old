#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
# generated by wxGlade 0.6.3 on Sun Mar 15 14:50:06 2009
import matplotlib
matplotlib.use( 'WXAgg',warn=False )
import traceback
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
import pylab
import copy
import csv
import wx
import wx.grid
import numpy
from editatombond import EditAtomTypes
from help import About
import os
import pickle
import sys
import settings


# begin wxGlade: extracode
# end wxGlade

ID_MENU_FILE_EXIT = 100
ID_MENU_FILE_NEW = 101
ID_MENU_FILE_OPEN = 102
ID_MENU_FILE_SAVE = 103

ID_CUT = 200
ID_COPY = 201
ID_PASTE = 202
ID_DEL = 203
ID_SELECTALL = 204


class MainFrame(wx.Frame):

    def __init__(self, fname,*args, **kwds):

        wx.Frame.__init__(self, *args, **kwds)

        self.level_pointer = 0
        #huckel solver and sketch_pad must be set up first since many of the controls and panels rely on it existing
        self.huckel_solver = huckelsolver.HuckelSolver()
        self.sketch_pad = MoleculePlotPanel(self, size=(500,300),style=wx.RAISED_BORDER)
        self.huckel_solver.addListener(self.sketch_pad.refreshFromHuckel)        


        self.results_display_2dmo = EigenPlotPanel(self,size=(400,300),style=wx.RAISED_BORDER)
        self.huckel_solver.addListener(self.results_display_2dmo.refreshFromHuckel)


        self.huckel_matrix = HuckelMatrix(self,size=(300, 300),label="Huckel Matrix")
        self.huckel_solver.addListener(self.huckel_matrix.refreshFromHuckel)
        
        self.controls = ControlPanel(self, label="Controls", size=(125,100))
        self.huckel_solver.addListener(self.controls.refreshFromHuckel)

        atype = self.controls.atom_type.GetValue().split('-')[0].strip()
        self.sketch_pad.current_atom_type = atype
        
        #create energy level diagram plot
        self.eld = ELDPlotPanel(self,size=(400,300),style=wx.RAISED_BORDER)
        self.huckel_solver.addListener(self.eld.refreshFromHuckel)

        #create results notebook
        self.results_display = wx.Notebook(self, -1, size=(400,300),style=wx.NO_BORDER)
        self.huckel_solver.addListener(self.handleResultsRefreshFromHuckel)


        self.panel_children = {}
        #eigenvalue/vectors tab 
        self.results_display_eig = wx.Panel(self.results_display,name="eigen panel") 
        self.eigen_matrix = EigenMatrix(self.results_display_eig,-1,size=(300,270),label="eigen matrix")
        self.panel_children["eigen panel"] = self.eigen_matrix
#        self.huckel_solver.addListener(self.eigen_matrix.refreshFromHuckel)

        #pibond orders
        self.results_display_pibond = wx.Panel(self.results_display,name="pibond panel")
        self.pibond_matrix = PiBondMatrix(self.results_display_pibond,-1,size=(300,270),label="pibond matrix")
        self.panel_children["pibond panel"] = self.pibond_matrix
#        self.huckel_solver.addListener(self.pibond_matrix.refreshFromHuckel)

        #net atomic charges
        self.results_display_charge = wx.Panel(self.results_display,name="charge panel")
        self.net_charge = NetChargeMatrix(self.results_display_charge,-1,size=(300,270),label="net charges")
        self.panel_children["charge panel"] = self.net_charge
#        self.huckel_solver.addListener(self.net_charge.refreshFromHuckel)
        
        #atom-atom and atom-bond polarizabilities tab
        self.results_display_pol = wx.Panel(self.results_display,name="atom atom panel")
        self.atom_atom_matrix = AtomAtomPolarizabilityMatrix(self.results_display_pol,-1,size=(300,270),label="atom atom")
        self.panel_children["atom atom panel"] = self.atom_atom_matrix
#        self.huckel_solver.addListener(self.atom_atom_matrix.refreshFromHuckel)

        #atom-bond polarizabilities tab
        self.results_display_atom_bond_pol = wx.Panel(self.results_display,name= "atom bond panel")
        self.atom_bond_matrix = AtomBondPolarizabilityMatrix(self.results_display_atom_bond_pol,-1,size=(300,270),label="atom bond")
        self.panel_children["atom bond panel"] = self.atom_bond_matrix
#        self.huckel_solver.addListener(self.atom_bond_matrix.refreshFromHuckel)

        self.huckel_solver.addListener(self.basisSizeChange)
        
        self.file_name = None
        self.save_required = False
        self.huckel_solver.addListener(self.setSaveRequired)
        
        self.SetTitle('Orbis - Simple Huckel Solver: Untitled')
        self.Bind(wx.EVT_KEY_DOWN,self.onKeyPressBasis)
#        self.Bind(wx.EVT_IDLE,self.huckel_solver.update)

        self.setupMenuBar()
        self.createToolBar()
        self.CreateStatusBar()
        
        self.doLayout()
        #self.checkTimeBomb()
        self.results_display.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED,self.OnNoteBookPage)        
        
        if fname != None:
            self.loadFile(fname)
        # end wxGlade
        
    def setSaveRequired(self):
        self.save_required = True
        
    def handleResultsRefreshFromHuckel(self):

        
        page_label = self.results_display.GetCurrentPage().GetName()

        self.panel_children[self.results_display.GetCurrentPage().GetName()].refreshFromHuckel()

    def OnNoteBookPage(self,event):

        page_label = self.results_display.GetPage(event.GetSelection()).GetName()

        page = self.panel_children[page_label]
        page.refreshFromHuckel()
        event.Skip()
#        self.eigen_matrix
#        if page_label == "eigen panel":
#            self.eigen_matrix.refreshFromHuckel()
            
#        refreshFromHuckel()
#        print page

#        self.results_display.GetPage(page).refreshFromHuckel()
        

    def checkTimeBomb(self):
        import datetime
        to = datetime.date.today()

        tf = datetime.date(*settings.timebomb)
        if int(str(tf - to).split(' ')[0])<0:
            wx.MessageBox("We're sorry but this Orbis beta version has expired. Please visit %s to purchase a full copy." % (settings.website),"Trial Expired",style=wx.OK)
            self.Close()

        
    def basisSizeChange(self):
        pass
        
    def setLevelPointer(self,level):
        if level != self.level_pointer:
            self.level_pointer = level
            self.results_display_2dmo.draw()
            
            self.eld.draw()
    
            self.eigen_matrix.SelectRow(self.huckel_solver.getSize()-level-1)

    def onKeyPressBasis(self,event):
        event.Skip()

    def createToolBar(self):   
        toolbar = self.CreateToolBar()                           
        for each in self.toolbarData():                          
            self.createSimpleTool(toolbar, *each)                
        toolbar.AddSeparator()                           
        toolbar.Realize()   

    def createSimpleTool(self, toolbar, label, filename, help, handler,id):                                  
        if not label:                                        
            toolbar.AddSeparator()                           
            return                  

        bmp = wx.Image(filename,wx.BITMAP_TYPE_PNG).Scale(20,20).ConvertToBitmap()

        tool = toolbar.AddSimpleTool(id, bmp, label, help)   
        self.Bind(wx.EVT_MENU, handler, tool,id=id)  


    def toolbarData(self): 
        return (("New", settings.image_dir+"/new.png", "Start a new session",self.OnNew,ID_MENU_FILE_NEW), 
                ("", "", "", "",-1), 
                ("Open", settings.image_dir+"/open.png", "Open previously saved session", self.OnLoad,ID_MENU_FILE_OPEN),
                ("Save", settings.image_dir+"/save.png", "Save current session", self.OnSave,ID_MENU_FILE_SAVE))

    
    def fileMenuData(self):
        return (('&New\tCtrl+N','Start a new session',self.OnNew),
                ('&Open\tCtrl+O','Open a previously saved session',self.OnLoad),
                ('&Save Session\tCtrl+S','Save the current Session',self.OnSave),
                ('Save Session &As...\tCtrl+Shift+S','Save the current session as...',self.OnSaveAs),
                (),
                ('Import &Geometry Data...\tCtrl+G','Import molecule geometry data from a file',self.OnImportGeom),
                ('Import &Huckel Determinant...\tCtrl+H','Import Huckel determinant data from a text file',self.OnImportHuck),
                ('&Export Results...\tCtrl+E','Export results of Huckel Calculation',self.OnExport),
                (),
                ('&Quit\tCtrl+Q' ,'Terminate Scimple Huckel Solver',self.OnQuit))

    def editMenuData(self):
        return (('&Copy\tCtrl+C','Copy selected data to the clipboard',self.OnCopy),)    

    def helpMenuData(self):
        return (('&About...','Information about Orbis',self.onAbout),
                ('&Documentation...','Information about the operation of Orbis',self.onDocs),
                ('&Visit the Orbis website...','Purchase the full version of Orbis or check the latest documentation',self.onWeb),                
                )    
    
    def viewMenuData(self):
        return (('&Visual Mode\tCtrl+I',"Switch between visual (molecule) mode and general determinant mode",self.onVisualMode,wx.ITEM_CHECK),
                ('&Redraw\tCtrl+R',"Redraw the current molecule",self.controls.onMinimize),
                ('&Zoom In\tCtrl+=',"Zoom in on molecule sketch pad",self.onZoomIn),
                ('&Zoom Out\tCtrl+-',"Zoom out of molecule sketch pad",self.onZoomOut),
                ('&Rotate Clockwise\tCtrl+Shift+=',"Rotate the molecule clockwiwe",self.onRotClock),
                ('&Rotate Counter Clockwise\tCtrl+Shift+-',"Rotate the molecule counter-clockwiwe",self.onRotCClock),
                )
    
    def onRotClock(self,event):
        self.sketch_pad.rotate(-1)
    def onRotCClock(self,event):
        self.sketch_pad.rotate(1)
    
    def onZoomIn(self,event):
        self.sketch_pad.resize(1)
    def onZoomOut(self,event):
        self.sketch_pad.resize(-1)

    def toolMenuData(self):
        return (('&Edit Atom Types...\tCtrl+A','Edit an existing atom type or add a new type of atom to the existing list',self.OnEditAtomType),
                ('&Edit Bond Types...\tCtrl+B','Edit an existing bond type or add a new type of bond to the existing list',self.OnEditBondType))
                
    def createMenu(self,data):
        menu = wx.Menu()
        for datum in data:
            if len(datum) == 0:
                menu.AppendSeparator()
            else:
                if len(datum) == 4:
                    item = wx.MenuItem(menu,-1,datum[0],datum[1],kind=datum[3])
                else:
                    item = wx.MenuItem(menu,-1,datum[0],datum[1])
                self.Bind(wx.EVT_MENU,datum[2],item)
                menu.AppendItem(item)
        return menu
    

    def setupMenuBar(self):
        menu_file = self.createMenu(self.fileMenuData())

        menu_edit = self.createMenu(self.editMenuData())
        
        menu_view = self.createMenu(self.viewMenuData())


        self.visual_mode = menu_view.GetMenuItems()[0]
        self.visual_mode.Check(True)
        self.redraw = menu_view.GetMenuItems()[1]
        self.zoomi = menu_view.GetMenuItems()[2]
        self.zoomo = menu_view.GetMenuItems()[3]
        self.rotc = menu_view.GetMenuItems()[4]
        self.rotcc = menu_view.GetMenuItems()[5]
        self.menu_view_items = [self.redraw,self.zoomi,self.zoomo,self.rotc,self.rotcc]

        menu_tools = self.createMenu(self.toolMenuData())
        menu_help = self.createMenu(self.helpMenuData())        
        
        menus = ((menu_file,"&File"),(menu_edit,"&Edit"),(menu_view,"&View"),(menu_tools,'&Tools'),(menu_help,'&Help'))

        menu_bar = wx.MenuBar()
        for menu,menu_str in menus:
            menu_bar.Append(menu,menu_str)

        self.SetMenuBar(menu_bar)
    
    def onAbout(self,event):
        About(self).ShowModal()

    def onDocs(self,event):
        try:
            import webbrowser
            webbrowser.open(settings.doc_root)
        except:
            wx.MessageBox('Unable to find the documentation files. Please visit www.simplehuckel.com/docs.html instead',"Can't find documentation",style=wx.ICON_ERROR)

    def onWeb(self,event):
        try:
            import webbrowser
            webbrowser.open(settings.website)
        except:
            wx.MessageBox('Unable to open %s' % (settings.website),"Can't find website",style=wx.ICON_ERROR)
            
    def OnEditAtomType(self,event):
        atom_dlg = EditAtomTypes()
        
        if atom_dlg.ShowModal() == wx.ID_OK:
            Atom.WriteAtomLib(atom_dlg.atom_types)
            prev = self.controls.atom_type.GetValue()
            new_names = Atom.AtomNames()
            self.controls.atom_type.SetItems(new_names)
            if prev in new_names:
            
                self.controls.atom_type.SetValue(prev)
            else:
                self.controls.atom_type.SetValue(new_names[0])
            
        atom_dlg.Destroy()
        
    def OnEditBondType(self,event):
        bond_dlg = EditBondTypes()
        
        if bond_dlg.ShowModal() == wx.ID_OK:
            Bond.WriteBondLib(bond_dlg.bond_types)
        bond_dlg.Destroy()

    def importHuck(self,fpath):
        try:

            csv_file = open(fpath,'r')
            row = csv_file.readline()
            dialect = csv.Sniffer().sniff(row)
            csv_file.seek(0)
            reader = csv.reader(csv_file,dialect)
            
        except:
            try:
                dialect.delimter = ' '
                csv_file.seek(0)
                reader = csv.reader(csv_file,dialect)
                
            except:
                return 0

        try:
            data = [map(float,row) for row in reader if len(row)>0]
            data = numpy.mat(data)
            
            if wx.MessageBox("Do you want to try to autodraw a molecule from this data?","Import as Molecule?",style=wx.ICON_QUESTION|wx.YES_NO) == wx.YES:
                self.setVisualMode()
    
                self.sketch_pad.reset()
                self.huckel_solver.reset()
                
                self.sketch_pad.createFromHuckel(data)

                self.huckel_solver.setData(data,data.shape[0])                                    
                
                self.sketch_pad.draw()
                self.eld.draw()
                self.results_display_2dmo.draw()
        
                self.sketch_pad.resize()                
            else:
                self.setVisualMode(False)
                self.huckel_solver.setData(data)


                
        except:
            return 0

        return 1            
    
    def importGeom(self,fname,parserType):
        try:
            self.visual_mode.Check(True)
            self.setVisualMode()                                        
    
            parser = parserType(fname)
            parser.parse()
            parser.close()
    
            data = numpy.mat(numpy.zeros((len(parser.atoms),len(parser.atoms)),float))
            
            for ii, atom in enumerate(parser.atoms):
                sym,x,y,z = atom
                self.sketch_pad.addNewAtom(x,y,sym=sym)
                atom = self.sketch_pad.molecule.atom_stack[-1]
                data[ii,ii] = atom.hx
    
            for a,b,t in parser.bonds:
                self.sketch_pad.addNewBond((a,b))
                k_xy = self.sketch_pad.molecule.bond_stack[-1].k_xy
                data[a,b] = k_xy
                data[b,a] = k_xy
    
            self.huckel_solver.setData(data,data.shape[0])                                    
            
            self.sketch_pad.draw()
            self.eld.draw()
            self.results_display_2dmo.draw()
    
            self.sketch_pad.resize()
            
            return 1
        except:
            return 0




    def OnImportGeom(self,event):
        import geomparsers
        file_type_handlers = {0:geomparsers.ParseCSV, 1:geomparsers.ParseMOL, 2:geomparsers.ParseXYZ, 3:geomparsers.ParseCSV}
        wildcard = "Comma Separated Value (*.csv)|*.csv|"\
                 "Molfile (*.mol)|*.mol|"\
                 "XYZ file (*.xyz)|*.xyz|"\
                 "All (*.*)|*.*"
        
        file_dlg = wx.FileDialog(self,"Import Molecule Geometry Data",wildcard=wildcard,style=wx.OPEN)

        if file_dlg.ShowModal() == wx.ID_OK:
            
            result = self.importGeom(file_dlg.GetPath(),file_type_handlers[file_dlg.GetFilterIndex()])
            if result == 0:
                err = "Unable to import data from %s: Please check the formatting of the file." %(file_dlg.GetFilename())                
                err += "\nIf you believe this is a bug in Orbis please send an email describing the problem to mail@simplehuckel.com"
                wx.MessageBox(err,"Import Error",style=wx.OK|wx.ICON_ERROR)                
                    

    def OnImportHuck(self,event):
        wildcard = "Comma Separated Value (*.csv)|*.csv|"\
                 "All (*.*)|*.*"
        
        file_dlg = wx.FileDialog(self,"Import Huckel Determinant Data",wildcard=wildcard,style=wx.OPEN)

        if file_dlg.ShowModal() == wx.ID_OK:
            
            result = self.importHuck(file_dlg.GetPath())
            if result == 0:
                err = "Unable to import data from %s: Please check the formatting of the file." %(file_dlg.GetFilename())
                err += "\nIf you believe this is a bug in Orbis please send an email describing the problem to mail@simplehuckel.com"
                wx.MessageBox(err,"Import Error",style=wx.OK|wx.ICON_ERROR)                
                    
    
                
    def exportData(self,fname):
        outfile = open(fname,'w')
        size = self.huckel_solver.getSize()
        fmt = "%.4G"
        outfile.write('Huckel Determinant\n')
        huckel_dat =self.huckel_solver.data.tolist()
        outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
        for ii,row in enumerate(huckel_dat):
            out = str(ii+1)+','+','.join([fmt % (x) for x in row])
            outfile.write(out+'\n')

        outfile.write('\n\nEigenvalues and vectors\n')
        outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
        for ii in range(size):
            
            out = 'E = %f' % (self.huckel_solver.eigen_vals[size-ii-1])
            out += ','+','.join([fmt % (x if abs(x)>settings.eps else 0) for x in self.huckel_solver.eigen_vecs[size-ii-1]])
            outfile.write(out+'\n')

        outfile.write('\n\nPi-Bond Orders\n')
        outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
        for ii,row in enumerate(self.huckel_solver.bond_orders.tolist()):
            out = str(ii+1)+','+','.join([fmt % (x if abs(x)>settings.eps else 0) for x in row])
            outfile.write(out+'\n')
            

        outfile.write('\n\nAtom-Atom Polarizabilities\n')
        outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
        for ii,row in enumerate(self.huckel_solver.aa_polar.tolist()):
            out = str(ii+1)+','+','.join([fmt % (x if abs(x)>settings.eps else 0) for x in row])
            outfile.write(out+'\n')

        
        outfile.write('\n\nAtom-Bond Polarizabilities for\n' )
        data = self.huckel_solver.ab_polar

        outfile.write(','+','.join(["Bond (%d-%d)" % (x+1,y+1)for x,y,z in data[0]])+'\n')

        for ii,row in enumerate(data):
            out = 'Atom '+str(ii+1)+','+','.join([fmt % (x[2] if abs(x[2])>settings.eps else 0 ) for x in row])
            outfile.write(out+'\n')
            
        outfile.close()
        
        
        
    def OnExport(self,event):

        file_dlg = wx.FileDialog(self,"Export to file",wildcard = "Comma Separated Values (*.csv)|*.csv|All (*.*)|*.*",style=wx.SAVE)
        
        fname = None
        result = None
        while fname == None and result != wx.ID_CANCEL:
            result = file_dlg.ShowModal()
            if result == wx.ID_OK:
                if os.path.exists(file_dlg.GetPath()):
                    overwrite_dlg = wx.MessageBox("This filename alread exists...Overwrite?","Warning",style=wx.ICON_WARNING|wx.YES_NO|wx.NO_DEFAULT)
                    if overwrite_dlg == wx.YES:
                        fname = file_dlg.GetPath()
                else:
                    fname = file_dlg.GetPath()
                 
                if fname != None:
                    self.exportData(fname)
                    
        file_dlg.Destroy()

    def loadFile(self,fpath):
#        try:
            f= open(fpath,'rb')
            session_dict = pickle.load(f)
            f.close()
            self.file_name = fpath

            if session_dict["mode"] == "visual":

                self.visual_mode.Check(True)
                self.setVisualMode()                                        

                n_a = len(session_dict["atoms"])
                n_b = len(session_dict["bonds"])
                data = numpy.mat(numpy.zeros((n_a,n_a)),float)
                for ii, atom in enumerate(session_dict["atoms"]):
                    x,y,hx,sym = atom
                    data[ii,ii] = hx
                        #make sure all atom types are present
#                        if sym not in Atom.ATOM_TYPES.keys():
#                            Atom.ATOM_TYPES[sym] = {settings.h_delim:hx,"description":sym}
                            
#                        self.sketch_pad.addNewAtom(x,y,hx,sym)

                for bond in session_dict["bonds"]:
                    con,k_xy = bond
                    data[con[0],con[1]] = k_xy
                    data[con[1],con[0]] = k_xy


                self.sketch_pad.createFromHuckel(data,False)

                self.huckel_solver.setData(data,session_dict["num_e"])                        
                
                for ii,props in enumerate(session_dict["atoms"]):
                    atom = self.sketch_pad.molecule.atom_stack[ii]
                    atom.setData(*props)

                    
                for ii,props in enumerate(session_dict["bonds"]):
                    con,k_xy = props
                    bond = self.sketch_pad.molecule.getBond(con[0],con[1])
                    bond.k_xy = k_xy 
                    bond.refresh()
                    
##                        make sure this bond exists int he bond types
                    #Bond.AddBond(a,b,k_xy)

                self.sketch_pad.draw()
                self.eld.draw()
                self.results_display_2dmo.draw()

                self.sketch_pad.resize()
            else:
                
                self.huckel_solver.setData(session_dict["data"],session_dict["num_e"])
                #self.visual_mode.Check(False)
                self.setVisualMode(False)

            self.SetTitle('Orbis - Simple Huckel Solver: %s' % (self.file_name))
#                self.huckel_solver.setNumElectrons(session_dict["num_e"])
                
#        except:
#            wx.MessageBox("Error while loading %s" % fpath,style = wx.ICON_ERROR)
                
        
    def OnLoad(self,event):
        
        file_dlg = wx.FileDialog(self,"Choose a saved session file",wildcard="*.huc")
        result = file_dlg.ShowModal()
        self.level_pointer = 0
#        self.setLevelPointer(0)

        if result == wx.ID_OK:
            self.loadFile(file_dlg.GetPath())

        file_dlg.Destroy()

        
    def setVisualMode(self,vis_mode=True):
        
        if vis_mode:
            self.visual_mode.Check(True)            
            self.sketch_pad.reset()


            self.huckel_solver.reset()


#            if self.sketch_pad.refreshFromHuckel not in self.huckel_solver.listeners:
#                self.huckel_solver.addListener(self.sketch_pad.refreshFromHuckel)

#            if self.results_display_2dmo.refreshFromHuckel not in self.huckel_solver.listeners:
#                self.huckel_solver.addListener(self.results_display_2dmo.refreshFromHuckel)
                
            for x in self.menu_view_items:
                x.Enable(True)
#            self.redraw.Enable(True)
            self.controls.minimize.Enable(True)
            self.main_sizer.Show(self.sketch_pad,recursive=True)
            self.results_sizer.Show(self.results_display_2dmo,recursive=True)
            self.controls.basis_size.Enable(False)
            self.controls.atom_type.Enable(True)
            if self.results_display.GetPageCount()<4:
                self.results_display.AddPage(self.results_display_pibond,u"\u03A0 - Bond Orders")
                self.results_display.AddPage(self.results_display_charge,"Charge Densities")
                self.results_display.AddPage(self.results_display_pol,"A-A Polarizabilities")
                self.results_display.AddPage(self.results_display_atom_bond_pol,"A-B Polarizabilities")
        else:
            self.visual_mode.Check(False)
 #           if self.sketch_pad.refreshFromHuckel in self.huckel_solver.listeners:
 #               self.huckel_solver.removeListener(self.sketch_pad.refreshFromHuckel)

  #          if self.results_display_2dmo.refreshFromHuckel in self.huckel_solver.listeners:
  #              self.huckel_solver.removeListener(self.results_display_2dmo.refreshFromHuckel)
                
#            self.redraw.Enable(False)
            for x in self.menu_view_items:
                x.Enable(False)

            self.controls.minimize.Enable(False)
            self.main_sizer.Hide(self.sketch_pad,recursive=True)
            self.results_sizer.Hide(self.results_display_2dmo,recursive=True)
            
            self.controls.basis_size.Enable(True)
            self.controls.atom_type.Enable(False)
            
            self.results_display.RemovePage(4)            
            self.results_display.RemovePage(3)
            self.results_display.RemovePage(2)
            self.results_display.RemovePage(1)
#            self.results_display_atom_bond_pol.Hide()
            self.results_display_pol.Hide()

            
        self.Layout()            

    def checkSaveRequired(self):
        if self.save_required:
            result = wx.MessageBox('Any unsaved changes you have made will be lost. Save now?',style=wx.YES_NO|wx.CANCEL)
            if  result == wx.YES:
                r2 = self.OnSaveAs(None)
                if  r2 == wx.ID_CANCEL:
                    result = wx.CANCEL

            if result in (wx.YES, wx.NO):
                self.save_required = False
        else:
            result = wx.NO

        return result
    
    def onVisualMode(self,event):
        result = self.checkSaveRequired()

        self.OnNew(None)        
        if result != wx.CANCEL:
            #print self.visual_mode.IsChecked()
            
            if self.visual_mode.IsChecked():
                self.setVisualMode()
                
            else:
                self.setVisualMode(False)

    
        
    def OnGridEditorCreated(self, event):
        """ Bind the kill focus event to the newly instantiated cell editor """
        editor = event.GetControl()
        editor.Bind(wx.EVT_KILL_FOCUS, self.OnKillFocus)
        event.Skip()


    def OnCopy(self,event):
        event.Skip()

       
    def saveFile(self,fname):
        try:
            import pickle
            session_dict ={"version":settings.version}
            session_dict["num_e"] = self.controls.num_e.GetValue()
    
            if self.visual_mode.IsChecked():
                session_dict["mode"] = "visual"
                atoms = []
                bonds = []
                molecule = self.sketch_pad.molecule            
                for atom in molecule.atom_stack:
                    atoms.append([atom.x,atom.y,atom.hx,atom.sym])
                for bond in molecule.bond_stack:
                    ai,bi = molecule.atom_stack.index(bond.a),molecule.atom_stack.index(bond.b)
                    bonds.append([(ai,bi),bond.k_xy])
            
                session_dict["atoms"] = atoms
                session_dict["bonds"] = bonds
                
            else:
                session_dict["mode"] = "general"
                session_dict["data"] = self.huckel_solver.data
                
            output = open(fname,'wb')
            pickle.dump(session_dict,output)
            output.close()
            
            self.file_name = fname
            self.save_required=False
        except:
            wx.MessageBox('There was a problem writing to %s. Please try again.' %(fname),style=wx.ICON_ERROR)
        

    def OnSaveAs(self,event):
        save_dlg = wx.FileDialog(self,"Save Session As",wildcard="*.huc",style=wx.FD_SAVE)
        result = None
        fname = None
        while fname == None and result != wx.ID_CANCEL:
            result = save_dlg.ShowModal()
            if result == wx.ID_OK:
                if os.path.exists(save_dlg.GetPath()):
                    overwrite_dlg = wx.MessageBox("This filename alread exists...Overwrite?","Warning",style=wx.ICON_WARNING|wx.YES_NO|wx.NO_DEFAULT)
                    if overwrite_dlg == wx.YES:
                        fname = save_dlg.GetPath()
                else:
                    fname = save_dlg.GetPath()
                 
                if fname != None:
                    self.saveFile(fname)
                    
        save_dlg.Destroy()
        if result != wx.ID_CANCEL:
            self.SetTitle('Orbis - Simple Huckel Solver: %s' % (self.file_name))
        return result
        
    
    def OnSave(self,event):
        if self.file_name == None:
            result = self.OnSaveAs(event)
        else:
            result = self.saveFile(self.file_name)
            if result != wx.CANCEL:
                self.SetTitle('Orbis - Simple Huckel Solver: %s' % (self.file_name))
        return result
    
    def OnNew(self,event):
        result = self.checkSaveRequired()
        if result in (wx.YES,wx.NO):
            self.file_name = None
            self.controls.onClear(None)
            self.SetTitle('Orbis - Simple Huckel Solver: Untitled')
        
    def OnQuit(self,event):
        
        self.Close()

    def doLayout(self):
        # begin wxGlade: MainFrame.__do_layout

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        input_sizer = wx.BoxSizer(wx.HORIZONTAL)
        input_sizer.Add(self.controls, 0.,flag=wx.EXPAND)
        input_sizer.Add(self.sketch_pad, 1,flag=wx.EXPAND)
        input_sizer.Add(self.huckel_matrix,1, flag=wx.EXPAND)
        main_sizer.Add(input_sizer, 1,flag=wx.EXPAND)

        results_sizer = wx.BoxSizer(wx.HORIZONTAL)
        results_sizer.Add(self.eld,1,wx.EXPAND,5)
        results_sizer.Add(self.results_display_2dmo,1,wx.EXPAND,5)

        self.results_display.AddPage(self.results_display_eig, "Eigenvectors")
        eig_sizer = wx.BoxSizer(wx.HORIZONTAL)
        eig_sizer.Add(self.eigen_matrix,1,wx.EXPAND)
        self.results_display_eig.SetSizer(eig_sizer)
        self.eigen_matrix.SetSize(self.results_display.GetPage(0).GetSize())

        self.results_display.AddPage(self.results_display_pibond,"Pi Bond Orders")
        pib_sizer = wx.BoxSizer(wx.HORIZONTAL)
        pib_sizer.Add(self.pibond_matrix,1,wx.EXPAND)
        self.results_display_pibond.SetSizer(pib_sizer)
        self.pibond_matrix.SetSize(self.results_display.GetPage(0).GetSize())

        self.results_display.AddPage(self.results_display_charge,"Charge Densities")
        nc_sizer = wx.BoxSizer(wx.HORIZONTAL)
        nc_sizer.Add(self.net_charge,1,wx.EXPAND)
        self.results_display_charge.SetSizer(nc_sizer)
        self.net_charge.SetSize(self.results_display.GetPage(0).GetSize())

        
        self.results_display.AddPage(self.results_display_pol, "A-A Polarizabilities")
        pol_sizer = wx.BoxSizer(wx.HORIZONTAL)
        pol_sizer.Add(self.atom_atom_matrix,1,wx.EXPAND)
        self.results_display_pol.SetSizer(pol_sizer)
        self.atom_atom_matrix.SetSize(self.results_display.GetPage(0).GetSize())        
        
        pol_sizer2 = wx.BoxSizer(wx.VERTICAL)

        pol_sizer2.Add(self.atom_bond_matrix,1,wx.EXPAND)

        self.results_display_atom_bond_pol.SetSizer(pol_sizer2)
        self.atom_bond_matrix.SetSize(self.results_display.GetPage(0).GetSize())        

        self.results_display.AddPage(self.results_display_atom_bond_pol, "A-B Polarizabilities")
        
        results_sizer.Add(self.results_display, 1, wx.EXPAND)


        main_sizer.Add(results_sizer, 1, wx.EXPAND, 0)

        self.SetSizer(main_sizer)
        main_sizer.Fit(self)
        self.main_sizer = main_sizer
        self.results_sizer = results_sizer
        self.setVisualMode(True) 
        self.Layout()
        self.Maximize()
        icon1 = wx.Icon(settings.icon_file, wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon1)        

def check_num_atoms_exceeded(num):
    if num  >= settings.max_num_atoms: 
        
        wx.MessageBox("We're sorry but this Orbis trial version only allows %d atoms. Please visit %s to purchase a full copy." % (settings.max_num_atoms,settings.website),"Max Atoms Exceeded",style=wx.OK)
        return True
    else:
        return False
    
if __name__ == "__main__":

    try:

        import settings
        import huckelsolver

        from molecule import Atom, Bond
        from sketchpad import MoleculePlotPanel
        from eigenplot import EigenPlotPanel
        from guiparts import *
        from editatombond import *

        import psyco
        psyco.full()
        
        a = Atom(0,0) #preload atom types
        b = Atom(1,1)
        
        bond = Bond(a,b)

        app = wx.PySimpleApp(0)
        wx.InitAllImageHandlers()
        if len(sys.argv)>1:
            fpath = sys.argv[1].strip()
        else:
            fpath = None

        main_frame = MainFrame(fpath,None, -1, "Orbis - Simple Huckel Solver")

        app.SetTopWindow(main_frame)
        
            
        main_frame.Show()
        app.MainLoop()
    except:
        print traceback.print_exc()        
        print sys.exc_info()

    finally:
        settings.logfile.close()