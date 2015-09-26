#!/usr/bin/python
# -*- coding: utf-8

"""GUI application for exact eigen-energy of the Kitaev model."""

__author__ = "Satoshi Morita"
__version__ = "2.0.0"

import wx
import kitaev

class Kitaev(wx.Frame):
    def __init__(self, *arg, **kwargs):
        super(Kitaev,self).__init__(*arg, **kwargs)
        self.InitUI()

    def InitUI(self):
        menu_bar = wx.MenuBar()
        sub_menu = wx.Menu()
        item = sub_menu.Append(wx.ID_EXIT, '&Quit', 'Quit application')
        self.Bind(wx.EVT_MENU, self.OnQuit, item)
        menu_bar.Append(sub_menu, '&File')
        self.SetMenuBar(menu_bar)

        self.status_bar = self.CreateStatusBar()

        panel = wx.Panel(self)

        panel1 = wx.Panel(panel)
        text0 = wx.StaticText(panel1, label="L1 =")
        text1 = wx.StaticText(panel1, label="L2 =")
        text2 = wx.StaticText(panel1, label="M  =")
        self.label0 = wx.StaticText(panel1, label=str(param[0]), name="L1")
        self.label1 = wx.StaticText(panel1, label=str(param[1]))
        self.label2 = wx.StaticText(panel1, label=str(param[2]))
        slider0 = wx.Slider(panel1, value=param[0], minValue=0, maxValue=24, name="L1")
        slider1 = wx.Slider(panel1, value=param[1], minValue=0, maxValue=24, name="L2")
        slider2 = wx.Slider(panel1, value=param[2], minValue=0, maxValue=24, name="M")
        slider0.Bind(wx.EVT_SCROLL, self.OnSliderScroll)
        slider1.Bind(wx.EVT_SCROLL, self.OnSliderScroll)
        slider2.Bind(wx.EVT_SCROLL, self.OnSliderScroll)
        grid1 = wx.FlexGridSizer(3,3,10,20)
        grid1.AddMany([(text0),(self.label0),(slider0,1,wx.EXPAND),
                       (text1),(self.label1),(slider1,1,wx.EXPAND),
                       (text2),(self.label2),(slider2,1,wx.EXPAND)])
        grid1.AddGrowableCol(2, 1)
        panel1.SetSizer(grid1)

        panel2 = wx.Panel(panel)
        text4 = wx.StaticText(panel2, label="Jz/Jx =")
        tctrl0 = wx.TextCtrl(panel2, value=str(param[3]))
        tctrl0.Bind(wx.EVT_TEXT, self.OnText)
        grid2 = wx.BoxSizer(wx.HORIZONTAL)
        grid2.Add(text4, 0, wx.RIGHT | wx.ALIGN_CENTER_HORIZONTAL, 10)
        grid2.Add(tctrl0, 1, wx.EXPAND)
        panel2.SetSizer(grid2)

        panel3 = wx.Panel(panel)
        text5 = wx.StaticText(panel3, label="Result")
        self.tc_result = wx.TextCtrl(panel3, value="", style=wx.TE_MULTILINE | wx.TE_READONLY)
        grid3 = wx.BoxSizer(wx.VERTICAL)
        grid3.Add(text5, 0, wx.BOTTOM, 5)
        grid3.Add(self.tc_result, 1, wx.EXPAND)
        panel3.SetSizer(grid3)

        panel4 = wx.Panel(panel)
        button0 = wx.Button(panel4, label='Calculate')
        button0.Bind(wx.EVT_BUTTON, self.OnButton)
        grid4 = wx.BoxSizer(wx.VERTICAL)
        grid4.Add(button0, 0, wx.ALIGN_RIGHT, 0)
        panel4.SetSizer(grid4)

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(panel1, 0, wx.BOTTOM | wx.EXPAND, 15)
        vbox.Add(panel2, 0, wx.BOTTOM | wx.EXPAND, 15)
        vbox.Add(panel3, 1, wx.BOTTOM | wx.EXPAND, 15)
        vbox.Add(panel4, 0, wx.BOTTOM | wx.EXPAND, 0)
        vbox_top = wx.BoxSizer(wx.VERTICAL)
        vbox_top.Add(vbox, 1, wx.ALL | wx.EXPAND, 20)
        panel.SetSizer(vbox_top)

        self.SetSize((440, 560))
        self.SetTitle('Kitaev model')
        self.Show(True)

    def OnQuit(self, e):
        self.Close()

    def OnSliderScroll(self, e):
        obj = e.GetEventObject()
        val = obj.GetValue()
        name = obj.GetName()
        if name=="L1": param[0] = int(val)
        elif name=="L2": param[1] = int(val)
        elif name=="M": param[2] = int(val)
        self.UpdateLabel()

    def OnText(self, e):
        obj = e.GetEventObject()
        try:
            val = float(obj.GetValue())
            param[3] = val
            self.status_bar.SetStatusText("Jz/Jx= {0}".format(val))
        except:
            self.status_bar.SetStatusText("could not convert string to float")

    def OnButton(self, e):
        obj = e.GetEventObject()
        label = obj.GetLabel()
        if label=="Calculate": self.PrintResult(param)

    def UpdateLabel(self):
        self.label0.SetLabel(str(param[0]))
        self.label1.SetLabel(str(param[1]))
        self.label2.SetLabel(str(param[2]))

    def PrintResult(self,param):
        if param[0]<1: param[0]=1
        if param[1]<1: param[1]=1
        param[2] = param[2]%param[0]
        self.status_bar.SetStatusText("Calculating... [l1,L2,M,Jz]={0}".format(param))
        try:
            ene = GetEnergyList(param)
        except:
            self.status_bar.SetStatusText("Unexpected error")
            return
        output = []
        output.append("# L1= {0}".format(param[0]))
        output.append("# L2= {0}".format(param[1]))
        output.append("# M= {0}".format(param[2]))
        output.append("# Jz/Jx= {0:.6e}".format(param[3]))
        output.append("# loop\tenergy")
        for l,val in enumerate(ene):
            output.append("{0}\t{1:+.12e}".format(l,val))
        self.tc_result.SetValue("\n".join(output))
        self.status_bar.SetStatusText("Done!")


def GetEnergyList(param):
    from kitaev import min_energy, Bond
    kitaev.L1 = param[0]
    kitaev.L2 = param[1]
    kitaev.M = param[2]
    kitaev.Jz = param[3]
    return [ min_energy(Bond(loop))[0] for loop in range(4) ]

def main():
    app = wx.App()
    Kitaev(None)
    app.MainLoop()

param = [6, 2, 2, 1.0]

if __name__ == '__main__':
    main()
