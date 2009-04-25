
from wx import ToolTip,EVT_IDLE
import time

class TimedToolTip(ToolTip):
    
    def __init__(self,duration,*args,**kwargs):
        
        ToolTip.__init__(self,*args,**kwargs)
        self.duration = duration
        self.start_time = None
        
        
    def Enable(self,*args,**kwargs):
        
        ToolTip.Enable(*args,**kwargs)
        if args[0] == True:
            self.start_time = time.time()
        else:
            self.start_time = None
            
        
    def _onIdle(self,event):

        if self.start_time != None:
            if time.time() - self.start_time > self.duration:
                self.Enable(False)
        event.Skip()