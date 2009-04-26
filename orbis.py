#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
# generated by wxGlade 0.6.3 on Sun Mar 15 14:50:06 2009
import matplotlib
matplotlib.use( 'WXAgg',warn=False )

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure
import pylab
import copy
import csv
import wx
import wx.grid
import numpy


import os
import sys


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

    def __init__(self, *args, **kwds):

        wx.Frame.__init__(self, *args, **kwds)

        self.level_pointer = 0
        #huckel solver and sketch_pad must be set up first since many of the controls and panels rely on it existing
        self.huckel_solver = huckelsolver.HuckelSolver()
        self.sketch_pad = MoleculePlotPanel(self, size=(400,300))
        self.huckel_solver.addListener(self.sketch_pad.refreshFromHuckel)        


        self.results_display_2dmo = EigenPlotPanel(self,size=(400,300))
        self.huckel_solver.addListener(self.results_display_2dmo.refreshFromHuckel)


        self.huckel_matrix = HuckelMatrix(self,size=(400, 300),label="Huckel Matrix")
        self.huckel_solver.addListener(self.huckel_matrix.refreshFromHuckel)
        
        self.controls = ControlPanel(self, label="Controls", size=(125,100))
        self.huckel_solver.addListener(self.controls.refreshFromHuckel)

        atype = self.controls.atom_type.GetValue().split('-')[0].strip()
        self.sketch_pad.current_atom_type = atype
        
        #create energy level diagram plot
        self.eld = ELDPlotPanel(self,size=(300,300))
        self.huckel_solver.addListener(self.eld.refreshFromHuckel)

        #create results notebook
        self.results_display = wx.Notebook(self, -1, size=(400,300),style=wx.NB_TOP)

        #eigenvalue/vectors tab 
        self.results_display_eig = wx.Panel(self.results_display) 
        self.eigen_matrix = EigenMatrix(self.results_display_eig,-1,size=(300,270),label="eigen matrix")


        self.huckel_solver.addListener(self.eigen_matrix.refreshFromHuckel)

        self.results_display_pibond = wx.Panel(self.results_display)
        self.pibond_matrix = PiBondMatrix(self.results_display_pibond,-1,size=(300,270),label="pibond matrix")
        self.huckel_solver.addListener(self.pibond_matrix.refreshFromHuckel)

        #atom-atom and atom-bond polarizabilities tab
        self.results_display_pol = wx.Panel(self.results_display)
        self.atom_atom_matrix = AtomAtomPolarizabilityMatrix(self.results_display_pol,-1,size=(300,270),label="atom atom")
        self.huckel_solver.addListener(self.atom_atom_matrix.refreshFromHuckel)

        #atom-bond polarizabilities tab
        self.results_display_atom_bond_pol = wx.Panel(self.results_display)

        self.atom_bond_matrix = AtomBondPolarizabilityMatrix(self.results_display_atom_bond_pol,-1,size=(300,270),label="atom bond")
        self.huckel_solver.addListener(self.atom_bond_matrix.refreshFromHuckel)


        
        self.atom_num = wx.SpinCtrl(self.results_display_atom_bond_pol,-1,"",min=1,max=1,initial=1,name="atom_num")
        
        self.atom_num.Bind(wx.EVT_KEY_UP,self.onKeyPressAtomNum)
        self.Bind(wx.EVT_SPINCTRL,self.onAtomNumChange,self.atom_num)
        self.huckel_solver.addListener(self.basisSizeChange)
        
        self.file_name = None
        
        self.Bind(wx.EVT_KEY_DOWN,self.onKeyPressBasis)
        self.Bind(wx.EVT_IDLE,self.huckel_solver.update)

        self.setupMenuBar()
        self.createToolBar()
        self.CreateStatusBar()
        
        self.doLayout()
        # end wxGlade

    def onKeyPressAtomNum(self,event):
        key = event.GetKeyCode()
        
        if key == wx.WXK_RETURN or key == wx.WXK_NUMPAD_ENTER:
            self.atom_bond_matrix.setAtomNum(self.atom_num.GetValue()-1)
        event.Skip()

    def basisSizeChange(self,event):
        
        self.atom_num.SetRange(1,self.huckel_solver.getSize())
        
        if self.atom_bond_matrix.atom_num>self.huckel_solver.getSize()-1:
            
            self.atom_bond_matrix.setAtomNum(0)
        
        
    def onAtomNumChange(self,event):
        self.atom_bond_matrix.setAtomNum(self.atom_num.GetValue()-1)
        event.Skip()
        
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
                ('&Import Data...','Import Huckel matrix data from a text file',self.OnImport),
                ('&Export Results...','Export results of Huckel Calculation',self.OnExport),
                (),
                ('&Quit\tCtrl+Q' ,'Terminate Scimple Huckel Solver',self.OnQuit))

    def editMenuData(self):
        return (('&Copy\tCtrl+C','Copy selected data to the clipboard',self.OnCopy),)    
    
    def viewMenuData(self):
        return (('&Visual Mode\tCtrl+I',"Switch between visual (molecule) mode and general matrix mode",self.onVisualMode,wx.ITEM_CHECK),
                ('&Redraw\tCtrl+R',"Redraw the current molecule",self.controls.onMinimize),)
                
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

        menu_tools = self.createMenu(self.toolMenuData())
        
        menus = ((menu_file,"&File"),(menu_edit,"&Edit"),(menu_view,"&View"),(menu_tools,'&Tools'))

        menu_bar = wx.MenuBar()
        for menu,menu_str in menus:
            menu_bar.Append(menu,menu_str)

        self.SetMenuBar(menu_bar)
    
        
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

        
    def OnImport(self,event):
        wildcard = "Comma Separated Value (*.csv)|*.csv|"\
                 "All (*.*)|*.*"
        
        file_dlg = wx.FileDialog(self,"Import Data",wildcard=wildcard,style=wx.OPEN)
        if file_dlg.ShowModal() == wx.ID_OK:
            
            try:
                
                csv_file = open(file_dlg.GetPath(),'r')
                row = csv_file.readline()
                dialect = csv.Sniffer().sniff(row)
                csv_file.seek(0)
                reader = csv.reader(csv_file,dialect)
                
                
                data = [map(float,row) for row in reader if len(row)>0]
                data = numpy.mat(data)
                
                if wx.MessageBox("Do you want to try to autodraw a molecule from this data?","Import as Molecule?",style=wx.ICON_QUESTION|wx.YES_NO) == wx.YES:
                    self.setVisualMode()

                    self.sketch_pad.reset()
                    self.huckel_solver.reset()
                    
                    self.sketch_pad.createFromHuckel(data)
                    
                else:
                    self.setVisualMode(False)
                    self.huckel_solver.setData(data)
                
            except:
                wx.MessageBox("Unable to import data from %s" %(file_dlg.GetFilename()),"Import Error",style=wx.OK|wx.ICON_ERROR)

                
    def exportData(self,fname):
        outfile = open(fname,'w')
        size = self.huckel_solver.getSize()
        
        outfile.write('Huckel Matrix\n')
        huckel_dat =self.huckel_solver.data.tolist()
        outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
        for ii,row in enumerate(huckel_dat):
            out = str(ii+1)+','+','.join(["%e" % (x) for x in row])
            outfile.write(out+'\n')

        outfile.write('\n\nEigenvalues and vectors\n')
        outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
        for ii in range(size):
            
            out = 'E = %f' % (self.huckel_solver.eigen_vals[size-ii-1])
            out += ','+','.join(["%e" % (x) for x in self.huckel_solver.eigen_vecs[size-ii-1]])
            outfile.write(out+'\n')

        outfile.write('\n\nPi-Bond Orders\n')
        outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
        for ii,row in enumerate(self.huckel_solver.bond_orders.tolist()):
            out = str(ii+1)+','+','.join(["%f" % (x) for x in row])
            outfile.write(out+'\n')
            

        outfile.write('\n\nAtom-Atom Polarizabilities\n')
        outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
        for ii,row in enumerate(self.huckel_solver.aa_polar.tolist()):
            out = str(ii+1)+','+','.join(["%e" % (x) for x in row])
            outfile.write(out+'\n')

        for jj in range(size):
            outfile.write('\n\nAtom-Bond Polarizabilities for Atom %d\n' %(jj+1))
            outfile.write(','+','.join(["%d" % (x+1) for x in range(size)])+'\n')
            for ii,row in enumerate(self.huckel_solver._calcSingleABPolarizability(jj).tolist()):
                out = str(ii+1)+','+','.join(["%e" % (x) for x in row])
                outfile.write(out+'\n')
            
            
        outfile.close()
        #for ii in range(self.huckel_solver.getSize()):
         #   csv_writer.
        
        
        
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

    def OnLoad(self,event):
        import pickle
        file_dlg = wx.FileDialog(self,"Choose a saved session file",wildcard="*.huc")
        result = file_dlg.ShowModal()

        if result == wx.ID_OK:

            #try:
                f= open(file_dlg.GetPath(),'rb')
                session_dict = pickle.load(f)
                f.close()
                self.file_name = file_dlg.GetPath()
                self.controls.onClear(None)
                if session_dict["mode"] == "visual":
                    
                    self.visual_mode.Check(True)
                    for atom in session_dict["atoms"]:
                        x,y,hx,sym = atom
                        #make sure all atom types are present
                        if sym not in Atom.ATOM_TYPES.keys():
                            Atom.ATOM_TYPES[sym] = {settings.h_delim:hx,"description":sym}
                            
                        self.sketch_pad.addNewAtom(x,y,hx,sym)
                    for bond in session_dict["bonds"]:
                        con,k_xy = bond
                        a = self.sketch_pad.molecule.atom_stack[con[0]].sym
                        b = self.sketch_pad.molecule.atom_stack[con[1]].sym
                        
                        Bond.AddBond(a,b,k_xy)
                        
                        self.sketch_pad.addNewBond(con,k_xy)

                    self.setVisualMode()
                    self.sketch_pad.resize()
                else:
                    
                    self.huckel_solver.setData(session_dict["data"])
                    #self.visual_mode.Check(False)
                    self.setVisualMode(False)

                self.huckel_solver.setNumElectrons(session_dict["num_e"])
                
            #except:
             #   wx.MessageBox("Error while loading %s" % file_dlg.GetFilename(),style = wx.ICON_ERROR)
                
        file_dlg.Destroy()
        
    def setVisualMode(self,vis_mode=True):
        if vis_mode:
        
            self.visual_mode.Check(True)
            self.huckel_solver.addListener(self.sketch_pad.refreshFromHuckel)
            self.huckel_solver.addListener(self.results_display_2dmo.refreshFromHuckel)
            self.redraw.Enable(True)
            self.controls.minimize.Enable(True)
            self.main_sizer.Show(self.sketch_pad,recursive=True)
            self.results_sizer.Show(self.results_display_2dmo,recursive=True)
            self.controls.basis_size.Enable(False)
            self.controls.atom_type.Enable(True)

        else:
            self.visual_mode.Check(False)
            self.huckel_solver.removeListener(self.sketch_pad.refreshFromHuckel)
            self.huckel_solver.removeListener(self.results_display_2dmo.refreshFromHuckel)
            self.redraw.Enable(False)
            self.controls.minimize.Enable(False)
            self.main_sizer.Hide(self.sketch_pad,recursive=True)
            self.results_sizer.Hide(self.results_display_2dmo,recursive=True)
            self.controls.basis_size.Enable(True)
            self.controls.atom_type.Enable(False)
        self.Layout()            
        
    def onVisualMode(self,event):

        if self.visual_mode.IsChecked():
            if wx.MessageBox("Switching to visual mode will clear the current session. Are you sure you want to continue?","Visual mode",style=wx.YES_NO|wx.ICON_QUESTION) == wx.YES:

                self.sketch_pad.reset()
                self.huckel_solver.reset()
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
        
    def OnSave(self,event):
        if self.file_name == None:
            self.OnSaveAs(event)
        else:
            self.saveFile(self.file_name)
            
    def OnNew(self,event):
        do_new = wx.MessageBox("Any unsaved work will be lost.  Do you want to start a new session?","New Session?",style=wx.YES_NO)
        if do_new == wx.YES:
            self.file_name = None
            self.controls.onClear(None)
        
        
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

        self.results_display.AddPage(self.results_display_pibond,u"\u03A0 - Bond Orders")
        pib_sizer = wx.BoxSizer(wx.HORIZONTAL)
        pib_sizer.Add(self.pibond_matrix,1,wx.EXPAND)
        self.results_display_pibond.SetSizer(pib_sizer)
        self.pibond_matrix.SetSize(self.results_display.GetPage(0).GetSize())

        self.results_display.AddPage(self.results_display_pol, "Atom-Atom Polarizabilities")
        pol_sizer = wx.BoxSizer(wx.HORIZONTAL)
        pol_sizer.Add(self.atom_atom_matrix,1,wx.EXPAND)
        self.results_display_pol.SetSizer(pol_sizer)
        self.atom_atom_matrix.SetSize(self.results_display.GetPage(0).GetSize())        
        
        pol_sizer2 = wx.BoxSizer(wx.VERTICAL)

        atom_num_box = wx.StaticBox(self.results_display_atom_bond_pol,label="Atom Number")
        atom_num_sizer = wx.StaticBoxSizer(atom_num_box,wx.VERTICAL)
        atom_num_sizer.Add(self.atom_num,flag=wx.EXPAND)
        
        pol_sizer2.Add(atom_num_sizer,0,wx.EXPAND)        
        pol_sizer2.Add(self.atom_bond_matrix,1,wx.EXPAND)

        self.results_display_atom_bond_pol.SetSizer(pol_sizer2)
        self.atom_bond_matrix.SetSize(self.results_display.GetPage(0).GetSize())        

        self.results_display.AddPage(self.results_display_atom_bond_pol, "Atom-Bond Polarizabilities")
        
        results_sizer.Add(self.results_display, 1, wx.EXPAND, 5)


        main_sizer.Add(results_sizer, 1, wx.EXPAND, 0)

        self.SetSizer(main_sizer)
        main_sizer.Fit(self)
        self.main_sizer = main_sizer
        self.results_sizer = results_sizer
        self.setVisualMode(True)
        self.Layout()
        #self.Maximize()
        
        # end wxGlade

# end of class MainFrame


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
        app = wx.PySimpleApp(0)
        wx.InitAllImageHandlers()
        main_frame = MainFrame(None, -1, "Simple Huckel Solver")
        app.SetTopWindow(main_frame)
        main_frame.Show()
        app.MainLoop()
        
    except:
        print sys.exc_info()
    finally:
        settings.logfile.close()