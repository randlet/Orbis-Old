import wx

class NumberValidator(wx.PyValidator):
    def __init__(self):
        wx.PyValidator.__init__(self)
        
    def Clone(self):
        return NumberValidator()
    
    def Validate(self,win):
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()
        
        try:
            val = float(text)
            textCtrl.Refresh() 
            return True
        except:
            
            wx.MessageBox("This field must be a number!", "Error") 
            textCtrl.SetBackgroundColour("pink") 
            textCtrl.SetFocus() 
            textCtrl.SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW)) 
            textCtrl.Refresh() 
            return False
    def TransferToWindow(self): 
        return True 
    def TransferFromWindow(self): 
        return True