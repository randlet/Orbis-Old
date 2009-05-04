import settings
import wx
import wx.html

class About(wx.Dialog): 
    text = ''' 
<html> 
<body bgcolor="#FFFFFF"> 
<center><table bgcolor="#DCE4F6" width="100%" cellspacing="0" 
cellpadding="0" border="1"> 
<tr> 
    <td align="center"><h1>Orbis - Simple Huckel Solver</h1></td> 
</tr> 
</table> 
</center> 
<p><b>Orbis</b> is a program for solving Simple Huckel Molecular Orbital (SHMO) problems.
For more information see www.simplehuckel.com.
</p> 
<p><b>Orbis</b> is brought to you by 
Heaviside Software, Copyright 
&copy; 2009.</p> 
</body> 
</html> 
''' 
    def __init__(self, parent): 
        wx.Dialog.__init__(self, parent, -1, 'About Orbis v%s' %(settings.version), 
                          size=(440, 400) ) 
        html = wx.html.HtmlWindow(self) 
        html.SetPage(self.text) 
        button = wx.Button(self, wx.ID_OK, "Okay") 
        sizer = wx.BoxSizer(wx.VERTICAL) 
        sizer.Add(html, 1, wx.EXPAND|wx.ALL, 5) 
        sizer.Add(button, 0, wx.ALIGN_CENTER|wx.ALL, 5) 
        self.SetSizer(sizer) 
        self.Layout()