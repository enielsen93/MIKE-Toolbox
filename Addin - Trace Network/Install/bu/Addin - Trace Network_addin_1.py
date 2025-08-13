import arcpy
import pythonaddins
import mikegraph

class GroupComboBoxClass1(object):
    """Implementation for Addin - Trace Network_addin.combobox (ComboBox)"""
    def __init__(self):
        self.items = ["item1", "item2"]
        self.editable = True
        self.enabled = True
        self.dropdownWidth = 'WWWWWWWWWWWWWWWWWW'
        self.width = 'WWWWWWWWWWWWWWWWWW'
    def onSelChange(self, selection):
        self.selected_group = selection
        layers_in_group = [layer for layer in arcpy.mapping.ListLayers(mxd) if GroupComboBoxClass1.selected_group + "\\" in layer.longName]
        
        
        
        self.msm_Node = [layer for layer in layers_in_group if u"Brønd" in layer]
        self.msm_Link = [layer for layer in layers_in_group if "Ledning" in layer]
        self.ms_Catchment = [layer for layer in layers_in_group if "Delopland" in layer]
        
        self.graph = mikegraph.MikeNetwork(self.msm_Node.workspacePath)
        pass
    def onEditChange(self, text):
        pass
    def onFocus(self, focused):
        if focused:
            self.mxd = arcpy.mapping.MapDocument("Current")
            self.df = arcpy.mapping.ListDataFrames(mxd)[0]
            
            group_layers = [layer for layer in arcpy.mapping.ListLayers(mxd) if layer.isGroupLayer]
            self.items = []
            if len(layers) != 0:
                for layer in group_layers:
                    self.items.append(layer.longName)
        pass
    def onEnter(self):
        pass
    def refresh(self):
        pass

class TraceDownstreamClass3(object):
    """Implementation for Addin - Trace Network_addin.button_1 (Button)"""
    def __init__(self):
        self.enabled = True
        self.checked = False
    def onClick(self):
        
        
        
        pass

class TraceUpstreamClass2(object):
    """Implementation for Addin - Trace Network_addin.button (Button)"""
    def __init__(self):
        self.enabled = True
        self.checked = False
    def onClick(self):
        pass