import wx
import copy
from molecule import Atom,Bond
import settings

from validators import NumberValidator

    
class EditAtomTypes(wx.Dialog):
    
    def __init__(self):
        width = 300
        height = 300
        vpitch = 50
        size = (width,height)
        wx.Dialog.__init__(self,None,-1, 'Edit Atom Types',size =size)
        
        cancel = wx.Button(self,wx.ID_CANCEL, "Cancel", pos = (width-85,height-60))
        ok = wx.Button(self,wx.ID_OK,"OK",pos=(width - 165,height-60))
        ok.SetDefault()
        
        self.atom_types = copy.deepcopy(Atom.ATOM_TYPES)
        
        self.atom_list = Atom.AtomNames(self.atom_types)
        
        self.atom_type_ctrl = wx.ListBox(self,-1,(10,10),choices = self.atom_list,style=wx.LB_SINGLE)
        self.atom_type_ctrl.Select(0)
        
        edit = wx.Button(self,-1,"Edit...",pos=(width-85,10))
        edit.Bind(wx.EVT_BUTTON,self.OnEdit)
        new = wx.Button(self,-1,"New...",pos=(width-85,40))
        new.Bind(wx.EVT_BUTTON,self.OnNew)

        delete = wx.Button(self,-1,"Delete...",pos=(width-85,70))
        delete.Bind(wx.EVT_BUTTON,self.OnDelete)

        restore = wx.Button(self,-1,"Restore...",pos=(width-85,100))
        restore.Bind(wx.EVT_BUTTON,self.OnRestore)

    
    def OnRestore(self,event):
        atom_type = self.atomTypeFromSelection()
        result = wx.MessageBox("Restoring will erase all custom atoms defined. Do you want to restore atom types to the defaults?",style=wx.ICON_WARNING|wx.YES_NO)
        if result == wx.YES:
            Atom._readAtomTypes(settings.default_atomic_data_file)
            self.atom_types = copy.deepcopy(Atom.ATOM_TYPES)
            self.refreshChoices()

            
    def showAtomType(self,name,atom_type,hx,des):
        edit = wx.Dialog(self,-1,'Select Properties for %s' % (name) ,size =(335,160))

        txt = wx.StaticText(edit, -1,u"Atom Type:",pos = (12,10))
        self.atom_type = wx.TextCtrl(edit,-1,atom_type,pos=(70,10),size=(25,-1))
        self.atom_type.SetMaxLength(2)
        
        
        txt = wx.StaticText(edit, -1,u"Set hx for atom type %s:" % (atom_type), style=wx.ALIGN_CENTRE,pos = (10,40))
        txt_size = txt.GetSize()
        pos = (txt_size[0]+15,txt.GetPosition()[1]-3)
        self.hx = wx.TextCtrl(edit, -1,pos=pos,validator = NumberValidator())
        
        self.hx.SetValue(str(hx))
        
        txt = wx.StaticText(edit, -1,u"Description:",pos = (12,70))
        self.des = wx.TextCtrl(edit,-1,des,pos=(70,70),size=(250,-1))
        self.des.SetMaxLength(30)

        ok_button = wx.Button(edit,wx.ID_OK,"OK",pos=(170,100))
        ok_button.SetDefault()
        cancel_button = wx.Button(edit,wx.ID_CANCEL, "Cancel", pos = (250,100))

        return edit

    def atomTypeFromSelection(self):
        return self.atom_list[self.atom_type_ctrl.GetSelection()].split('-')[0].strip()
    
    def OnEdit(self,event):
        name = self.atom_list[self.atom_type_ctrl.GetSelection()]
        atom_type = self.atomTypeFromSelection()
        hx = self.atom_types[atom_type][settings.h_delim]
        des = self.atom_types[atom_type]["description"]
        edit = self.showAtomType(name,atom_type,hx,des)
        self.atom_type.SetEditable(False)
        if edit.ShowModal() ==wx.ID_OK:
            self.atom_types[atom_type][settings.h_delim] = self.hx.GetValue()
            self.atom_types[atom_type]["description"] = self.des.GetValue().encode('ascii')
            
        edit.Destroy()

    def refreshChoices(self):
        self.atom_list = Atom.AtomNames(self.atom_types)
        self.atom_type_ctrl.Set(self.atom_list)
        self.atom_type_ctrl.SetSelection(0)

    def OnNew(self,event):
        new = self.showAtomType("New","",0.0,"New Atom Type")
        self.atom_type.SetEditable(True)
        result = None
        
        while result == None :
            result = new.ShowModal()
            if result == wx.ID_OK:
                val = self.atom_type.GetValue().strip()
                if val not in self.atom_types.keys() and len(val)>0:
                    
                    self.atom_types[self.atom_type.GetValue().encode('ascii')] = {settings.h_delim:float(self.hx.GetValue()),"description":self.des.GetValue().encode('ascii')}
                    self.refreshChoices()
                else:
                    if len(val) == 0:
                        msg = "Please specify an atom type!"
                    else:
                        msg = "An Atom of that type already exists! Please choose another type name"
                    result = None
                    wx.MessageBox(msg,'Create Atom Problem',style=wx.ICON_WARNING)

        return result
    def OnDelete(self,event):
        atom_type = self.atomTypeFromSelection()
        result = wx.MessageBox("Are you sure you want to delete this atom type?",'Delete Atom Type',style=wx.ICON_WARNING|wx.YES_NO)
        if result == wx.YES:
            self.atom_types.pop(atom_type)
            self.refreshChoices()


            
class EditBondTypes(wx.Dialog):
    
    def __init__(self):
        width = 300
        height = 300
        vpitch = 50
        size = (width,height)
        wx.Dialog.__init__(self,None,-1, 'Edit Bond Types',size =size)
        
        cancel = wx.Button(self,wx.ID_CANCEL, "Cancel", pos = (width-85,height-60))
        ok = wx.Button(self,wx.ID_OK,"OK",pos=(width - 165,height-60))
        ok.SetDefault()
        
        self.bond_types = copy.deepcopy(Bond.BOND_TYPES)
        
        
        self.bond_list = Bond.BondNames(self.bond_types)
        
        self.bond_type_ctrl = wx.ListBox(self,-1,(10,10),choices = self.bond_list,style=wx.LB_SINGLE)
        self.bond_type_ctrl.Select(0)
        
        edit = wx.Button(self,-1,"Edit...",pos=(width-85,10))
        edit.Bind(wx.EVT_BUTTON,self.OnEdit)
        new = wx.Button(self,-1,"New...",pos=(width-85,40))
        new.Bind(wx.EVT_BUTTON,self.OnNew)

        delete = wx.Button(self,-1,"Delete...",pos=(width-85,70))
        delete.Bind(wx.EVT_BUTTON,self.OnDelete)

        restore = wx.Button(self,-1,"Restore...",pos=(width-85,100))
        restore.Bind(wx.EVT_BUTTON,self.OnRestore)

#    def bondTypeChoices(self):
#        return ["%s - %s" % (x,self.bond_types[x]["description"]) for x in sorted(self.bond_types.keys())]
    
    def OnRestore(self,event):
        bond_type = self.bondTypeFromSelection()
        result = wx.MessageBox("Restoring will erase all custom bonds defined. Do you want to restore bond types to the defaults?",style=wx.ICON_WARNING|wx.YES_NO)
        if result == wx.YES:
            Bond._readBondTypes(settings.default_bond_data_file)
            self.bond_types = copy.deepcopy(Bond.BOND_TYPES)
            self.refreshChoices()

        
    def showBondType(self,name,bond_type,kxy):
        edit = wx.Dialog(self,-1,'Select Properties for %s' % (name) ,size =(335,160))
        atom_list = sorted(Atom.ATOM_TYPES.keys())

        txt = wx.StaticText(edit, -1,u"Atom Type 1:",pos = (12,10))
        self.atom_type1 = wx.ComboBox(edit,-1,atom_list[0],pos=(80,10),choices=atom_list,style=wx.CB_DROPDOWN|wx.CB_READONLY)

        txt = wx.StaticText(edit, -1,u"Atom Type 2:",pos = (130,10))
        self.atom_type2 = wx.ComboBox(edit,-1,atom_list[0],pos=(210,10),choices=atom_list,style=wx.CB_DROPDOWN|wx.CB_READONLY)
        
        txt = wx.StaticText(edit, -1,u"Set k_xy for bond type", style=wx.ALIGN_CENTRE,pos = (10,40))
        txt_size = txt.GetSize()
        pos = (txt_size[0]+15,txt.GetPosition()[1]-3)
        self.kxy = wx.TextCtrl(edit, -1,pos=pos,validator=NumberValidator())
        
        self.kxy.SetValue(str(kxy))
        
        ok_button = wx.Button(edit,wx.ID_OK,"OK",pos=(170,100))
        ok_button.SetDefault()
        cancel_button = wx.Button(edit,wx.ID_CANCEL, "Cancel", pos = (250,100))

        return edit

    def bondTypeFromSelection(self):
        return self.bond_list[self.bond_type_ctrl.GetSelection()]
    
    def OnEdit(self,event):
        name = self.bond_list[self.bond_type_ctrl.GetSelection()]
        bond_type = name
        kxy = self.bond_types[bond_type][settings.k_delim]
        
        edit = self.showBondType(name,bond_type,kxy)
        a,b = bond_type.split('-')
        
        self.atom_type1.SetItems([a])
        self.atom_type2.SetItems([b])
        self.atom_type1.SetValue(a)
        self.atom_type2.SetValue(b)
        
        if edit.ShowModal() ==wx.ID_OK:
            self.bond_types[bond_type][settings.k_delim] = self.kxy.GetValue()
            
        edit.Destroy()

    def refreshChoices(self):
        self.bond_list = Bond.BondNames(self.bond_types)
        self.bond_type_ctrl.Set(self.bond_list)
        self.bond_type_ctrl.SetSelection(0)

    def sketchNew(self,a,b):
        new = self.showBondType("New","",0.0)
        self.atom_type1.SetItems([a])
        self.atom_type2.SetItems([b])
        self.atom_type1.SetValue(a)
        self.atom_type2.SetValue(b)
        return self.createNew(new)
    
    def createNew(self,new):
        
        result = None
        
        while result == None :
            result = new.ShowModal()
            if result == wx.ID_OK:
                b1 = '%s-%s' % (self.atom_type1.GetValue(),self.atom_type2.GetValue())
                b2 = '%s-%s' % (self.atom_type2.GetValue(),self.atom_type1.GetValue())
                b1 = b1.encode('ascii')
                b2 = b2.encode('ascii')
                
                if b1 not in self.bond_types.keys() and b2 not in self.bond_types.keys():
                    self.bond_type = b1
                    
                    self.bond_types[self.bond_type] = {settings.k_delim:float(self.kxy.GetValue())}
                    self.refreshChoices()
                else:
                    result = None
                    wx.MessageBox("A bond of that type already exists! Please choose another type name",'Bond Already Exists',style=wx.ICON_WARNING)
        new.Destroy()
        return result

    def OnNew(self,event):
        new = self.showBondType("New","",0.0)
        self.createNew(new)
        
    def OnDelete(self,event):
        bond_type = self.bondTypeFromSelection()
        
        result = wx.MessageBox("Are you sure you want to delete this bond type?",'Delete Bond Type',style=wx.ICON_WARNING|wx.YES_NO)
        if result == wx.YES:
            self.bond_types.pop(bond_type)
            self.refreshChoices()

            