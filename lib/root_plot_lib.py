"""@package docstring
Python package to make nice-looking ROOT plots
"""

from array import array
#from os import path
import ROOT
import ctypes

def get_palette_hig21001():
  '''Returns list of colors used in HIG-21-001
  '''
  return [ROOT.TColor.GetColor('#de5a6a'), ROOT.TColor.GetColor('#ffcc66'), 
          ROOT.TColor.GetColor('#c0e684'), ROOT.TColor.GetColor('#64c0e8'), 
          ROOT.TColor.GetColor('#9999cc'), ROOT.TColor.GetColor('#ffccff')]

def get_palette_official(nplots):
  '''Returns list of colors to be used for plots

  @params
  nplots - number of different overlayed plots
  '''
  if (nplots <= 6):
    return [ROOT.TColor.GetColor('#5790fc'), ROOT.TColor.GetColor('#f89c20'), 
            ROOT.TColor.GetColor('#e42536'), ROOT.TColor.GetColor('#964a8b'), 
            ROOT.TColor.GetColor('#9c9ca1'), ROOT.TColor.GetColor('#7a21dd')]
  elif (nplots <= 10):
    return [ROOT.TColor.GetColor('#3f90da'), ROOT.TColor.GetColor('#ffa90e'), 
            ROOT.TColor.GetColor('#bd1f01'), ROOT.TColor.GetColor('#832db6'), 
            ROOT.TColor.GetColor('#94a4a2'), ROOT.TColor.GetColor('#a96b59'),
            ROOT.TColor.GetColor('#e76300'), ROOT.TColor.GetColor('#b9ac70'),
            ROOT.TColor.GetColor('#717581'), ROOT.TColor.GetColor('#92dadd')]
  raise ValueError('More colors than allowed in official colors.')

def get_palette_lines(nplots):
  '''Returns list of colors to be used for plots

  @params
  nplots - number of different overlayed plots
  '''
  if (nplots <= 6):
    return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ffbd38'), 
            ROOT.TColor.GetColor('#aee82a'), ROOT.TColor.GetColor('#28aee8'), 
            ROOT.TColor.GetColor('#446bcc'), ROOT.TColor.GetColor('#7346cc')]
  if (nplots <= 8):
    return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ff9c38'), 
            ROOT.TColor.GetColor('#f6c24b'), ROOT.TColor.GetColor('#aee82a'), 
            ROOT.TColor.GetColor('#28aee8'), ROOT.TColor.GetColor('#446bcc'), 
            ROOT.TColor.GetColor('#7346cc'), ROOT.TColor.GetColor('#ee66ac')]
  return [ROOT.TColor.GetColor('#de394d'), ROOT.TColor.GetColor('#ff9c38'), 
          ROOT.TColor.GetColor('#f6c24b'), ROOT.TColor.GetColor('#aee82a'), 
          ROOT.TColor.GetColor('#1b9e48'), ROOT.TColor.GetColor('#28aee8'), 
          ROOT.TColor.GetColor('#446bcc'), ROOT.TColor.GetColor('#50378f'),
          ROOT.TColor.GetColor('#7346cc'), ROOT.TColor.GetColor('#ee66ac')]

def get_graph_edges(graph):
  '''Returns tuple (low_edge, high_edge) of graph

  @params
  graph - TGraph to check
  '''
  lo_edge = 999999.9
  hi_edge = -999999.9
  for point in range(graph.GetN()):
    x = graph.GetPointX(point)
    exlo = graph.GetErrorXlow(point)
    exhi = graph.GetErrorXhigh(point)
    if (x-exlo)<lo_edge:
      lo_edge = x-exlo
    if (x+exhi)>hi_edge:
      hi_edge = x+exhi
  return (lo_edge, hi_edge)

class RplPlot:
  '''Simple class to hold plot options and apply to ROOT plots
  '''
  def __init__(self):
    '''Initializer to set defaults
    '''
    self.plot_bottom = False
    self.y_title = 'Events/bin'
    self.y_title_lower = 'Data/MC'
    self.x_title = 'x variable'
    self.z_title = 'Events/bin'
    self.x_min = -999.0
    self.x_max = -999.0
    self.y_min = -999.0
    self.y_max = -999.0
    self.y_min_lower = 0.75
    self.y_max_lower = 1.25
    self.log_x = False
    self.log_y = False
    self.log_y_bottom = False
    self.lumi_data = [(138,13)] #list of tuples in format (lumi [fbinv], energy [TeV])
    self.title_type = 'cms work in progress'
    self.legend_xlo = 0.19
    self.legend_xhi = 0.91
    self.legend_ylo = 0.78
    self.legend_yhi = 0.90
    self.legend_customsize = False
    self.legend_ncolumns = -1
    self.hists = []
    self.hist_color = []
    self.hist_style = [] #point, outline, filled
    self.color_index = 0
    self.palette = get_palette_official(6)
    self.graphs = []
    self.graph_color = []
    self.bottom_plots = []
    self.bottom_is_ratio = True
    self.bottom_plot_color = []
    self.n_plots = 0
    self.is_2d = False

  def plot_hist(self, hist, color=None, style='point'):
    '''plots a TH1

    @params
    hist - root TH1 to plot, title is used as legend entry
    color - custom color to use for plot, None uses a default
    style - point, outline, filled. Sets style of plot
    '''
    if self.is_2d:
      raise RuntimeError('Cannot add both 2D and 1D plots')
    if (color==None):
      if (self.color_index >= 6):
        raise RuntimeError('Error: too many plots for selected palette')
      self.hist_color.append(self.palette[self.color_index])
      self.color_index += 1
    else:
      self.hist_color.append(color)
    if (style not in ['point','outline','filled']):
      raise ValueError('Unsupported plot style')
    self.hist_style.append(style)
    self.hists.append(hist)
    if self.n_plots == 0:
      xaxis_title = hist.GetXaxis().GetTitle()
      if self.x_title == 'x variable' and xaxis_title != '':
        self.x_title = xaxis_title
    self.n_plots += 1
    return self

  def plot_points(self, hist, color=None):
    '''plots a TH1 as data points

    @params
    hist - root TH1 to plot, title is used as legend entry
    color - custom color to use for plot, None uses a default
    '''
    return self.plot_hist(hist, color, 'point')

  def plot_outline(self, hist, color=None):
    '''plots a TH1 as an outline

    @params
    hist - root TH1 to plot, title is used as legend entry
    color - custom color to use for plot, None uses a default
    '''
    return self.plot_hist(hist, color, 'outline')

  def plot_filled(self, hist, color=None):
    return self.plot_hist(hist, color, 'filled')

  def plot_graph(self, graph, color=None):
    '''plots a TGraph

    @params
    graph - root TGraph to plot, title is used as legend entry
    color - custom color to use for plot, None uses a default
    '''
    if self.is_2d:
      raise RuntimeError('Cannot add both 2D and 1D plots')
    if (color==None):
      if (self.color_index >= 6):
        raise RuntimeError('Error: too many plots for selected palette')
      self.graph_color.append(self.palette[self.color_index])
      self.color_index += 1
    else:
      self.graph_color.append(color)
    self.graphs.append(graph)
    if self.n_plots == 0:
      xaxis_title = graph.GetXaxis().GetTitle()
      if self.x_title == 'x variable' and xaxis_title != '':
        self.x_title = xaxis_title
    self.n_plots += 1
    return self

  def plot_colormap(self,hist):
    '''plots a TH2 as a heatmap
    hist  TH2 to plot
    '''
    if len(self.hists) != 0 or len(self.graphs) != 0:
      raise RuntimeError('Cannot add both 2D and 1D plots')
    self.hists.append(hist)
    self.is_2d = True
    self.hist_color.append(0)
    self.hist_style.append('colormap')
    if self.n_plots == 0:
      xaxis_title = hist.GetXaxis().GetTitle()
      if self.x_title == 'x variable' and xaxis_title != '':
        self.x_title = xaxis_title
      yaxis_title = hist.GetYaxis().GetTitle()
      if self.y_title == 'Events/bin' and yaxis_title != '':
        self.y_title = yaxis_title
      zaxis_title = hist.GetZaxis().GetTitle()
      if self.z_title == 'Events/bin' and zaxis_title != '':
        self.z_title = zaxis_title
    self.n_plots += 1

  def add_ratio(self, numerator_name, denominator_name, switch_colors=False):
    '''adds a ratio of two plot elements to the plot

    @params
    numerator_name    name of plot element to use as numerator
    denominator_name  name of plot element to use as denominator
    switch_colors     switch colors of numerator and denominator
    '''
    numerator_hist = None
    denominator_hist = None
    numerator_color = 1
    denominator_color = 1
    for hist, color in zip(self.hists, self.hist_color):
      if (hist.GetName()==numerator_name):
        numerator_hist = hist
        numerator_color = color
      if (hist.GetName()==denominator_name):
        denominator_hist = hist
        denominator_color = color
    if (numerator_hist != None and denominator_hist != None):
      #add uncertainties for reference on first plot
      if len(self.bottom_plots)==0:
        self.bottom_plots.append(denominator_hist.Clone())
        self.bottom_plots[-1].Divide(denominator_hist)
        if not switch_colors:
          self.bottom_plot_color.append(numerator_color)
        else:
          self.bottom_plot_color.append(denominator_color)
      self.bottom_plots.append(numerator_hist.Clone())
      self.bottom_plots[-1].Divide(denominator_hist)
      if not switch_colors:
        self.bottom_plot_color.append(denominator_color)
      else:
        self.bottom_plot_color.append(numerator_color)
      self.plot_bottom = True
    else:
      raise ValueError('Could not find numerator and denominator plots requested for ratio.')
    #TODO implement ratios of graphs
    return self

  def add_difference(self, minuend_name, subtrahend_name, switch_colors=False):
    '''adds a ratio of two plot elements to the plot

    @params
    minuend_name      name of plot element to use as minuend
    subtrahend_name   name of plot element to use as subtrahend
    switch_colors     switch colors of numerator and denominator
    '''
    self.bottom_is_ratio = False
    self.y_title_lower = 'Data-MC'
    minuend_hist = None
    subtrahend_hist = None
    minuend_color = 1
    subtrahend_color = 1
    for hist, color in zip(self.hists, self.hist_color):
      if (hist.GetName()==minuend_name):
        minuend_hist = hist
        minuend_color = color
      if (hist.GetName()==subrahend_name):
        subrahend_hist = hist
        subrahend_color = color
    if (minuend_hist != None and subtrahend_hist != None):
      #add uncertainties for reference on first plot
      if len(self.bottom_plots)==0:
        self.bottom_plots.append(strahend_hist.Clone())
        self.bottom_plots[-1].Add(subtrahend_hist,-1.0)
        if not switch_colors:
          self.bottom_plot_color.append(minuend_color)
        else:
          self.bottom_plot_color.append(subtrahend_color)
      self.bottom_plots.append(minuend_hist.Clone())
      self.bottom_plots[-1].Add(subtrahend_hist, -1.0)
      if not switch_colors:
        self.bottom_plot_color.append(subtrahend_color)
      else:
        self.bottom_plot_color.append(minuend_color)
      self.plot_bottom = True
    else:
      raise ValueError('Could not find minuend and subtrahend plots requested for difference.')
    #TODO implement differences of graphs
    return self

  def draw(self, filename='my_plot.pdf'):
    '''draws plot and saves to output file

    @params
    filename - name of plot to save file
    '''
    ROOT.gStyle.SetOptStat(0)

    if not self.is_2d:
      #find maxima and minima
      if (self.y_max == -999.0 or self.y_min == -999.0):
        self.y_min = 0
        if self.log_y:
          self.y_min = 0.01
        self.y_max = 0 
        for hist in self.hists:
          if hist.GetMaximum()>self.y_max:
            self.y_max = hist.GetMaximum()
        for graph in self.graphs:
          graph_max = ROOT.TMath.MaxElement(graph.GetN(), graph.GetY())
          if graph_max>self.y_max:
            self.y_max = graph_max
        if self.log_y:
          self.y_max = ((self.y_max/self.y_min)**1.45)*self.y_min
        else:
          self.y_max = self.y_max*1.45
      if (self.x_max == -999.0 or self.x_min == -999.0):
        self.x_min = 999999.0
        self.x_max = -999999.0
        #TODO implement something more robust, for now, just use first hist
        for hist in self.hists:
          lo_edge = hist.GetXaxis().GetBinLowEdge(1)
          hi_edge = hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins())
          if (lo_edge < self.x_min):
            self.x_min = lo_edge
          if (hi_edge > self.x_max):
            self.x_max = hi_edge
        for graph in self.graphs:
          lo_edge, hi_edge = get_graph_edges(graph)
          #lo_edge = ROOT.TMath.MinElement(graph.GetN(), graph.GetX())
          #hi_edge = ROOT.TMath.MaxElement(graph.GetN(), graph.GetX())
          if (lo_edge < self.x_min):
            self.x_min = lo_edge
          if (hi_edge > self.x_max):
            self.x_max = hi_edge

      #set up dummy histograms
      dummy_hist_upper = ROOT.TH1D('','',100,self.x_min,self.x_max)
      dummy_hist_upper.SetMinimum(self.y_min)
      dummy_hist_upper.SetMaximum(self.y_max)
      dummy_hist_upper.SetLabelSize(0.028,'y')
      dummy_hist_upper.SetTitleSize(0.032,'y')
      dummy_hist_upper.GetYaxis().SetTitle(self.y_title)
      ROOT.TGaxis.SetExponentOffset(-0.05,0.0,'y')
      dummy_hist_upper.GetYaxis().SetNdivisions(606)
      dummy_hist_upper.GetXaxis().SetNdivisions(606)
      if (not self.plot_bottom):
        dummy_hist_upper.GetXaxis().SetTitle(self.x_title)
        dummy_hist_upper.SetTitleSize(0.032,'x')
      else:
        dummy_hist_upper.SetLabelSize(0,'x')

      dummy_hist_lower = ROOT.TH1D('','',100,self.x_min,self.x_max)
      dummy_hist_lower.SetMinimum(self.y_min_lower)
      dummy_hist_lower.SetMaximum(self.y_max_lower)
      dummy_hist_lower.SetLabelSize(0.028,'x')
      dummy_hist_lower.SetLabelSize(0.028,'y')
      dummy_hist_lower.SetTitleOffset(1.5,'y')
      #dummy_hist_lower.SetTitleOffset(2.0,'y')
      #dummy_hist_lower.SetTitleSize(0.45/float(len(self.y_title_lower)),'y')
      dummy_hist_lower.SetTitleSize(0.032,'x')
      dummy_hist_lower.GetYaxis().SetTitle(self.y_title_lower)
      dummy_hist_lower.GetYaxis().SetNdivisions(606)
      dummy_hist_lower.GetYaxis().SetLimits(self.y_min_lower,self.y_max_lower)
      dummy_hist_lower.GetXaxis().SetTitle(self.x_title)
      dummy_hist_lower.GetXaxis().SetNdivisions(606)

    #set up pads
    canvas_name = filename[:filename.rfind('.')]
    if canvas_name.rfind('/') != -1:
      canvas_name = canvas_name[canvas_name.rfind('/')+1:]
    can = ROOT.TCanvas('c_'+canvas_name,'c',600,600)
    top_pad = ROOT.TPad('top_pad','',0.0,0.0,1.0,1.0)
    top_pad.SetTicks(1,1)
    right_margin = 0.06
    left_margin = 0.14
    bottom_margin = 0.1
    top_margin = 0.05
    if (self.is_2d):
      bottom_margin = 0.2
      right_margin = 0.18
    if (self.plot_bottom):
      bottom_margin = 0.31
    top_pad.SetMargin(left_margin,right_margin,bottom_margin,top_margin)
    top_pad.SetFillStyle(4000)
    top_pad.SetLogy(self.log_y)
    top_pad.SetLogx(self.log_x)

    bot_pad = ROOT.TPad('bot_pad','',0.0,0.0,1.0,1.0)
    bot_pad.SetTicks(1,1)
    bottom_margin = 0.0
    top_margin = 1.0
    if (self.plot_bottom):
      bottom_margin = 0.1
      top_margin = 0.71
    bot_pad.SetMargin(left_margin,right_margin,bottom_margin,top_margin)
    bot_pad.SetFillStyle(4000)
    bot_pad.SetLogy(self.log_y_bottom)
    bot_pad.SetLogx(self.log_x)

    #draw upper plot
    top_pad.Draw()
    top_pad.cd()

    #draw plots and legend
    if not self.is_2d:
      dummy_hist_upper.Draw()
      if (self.legend_ncolumns == -1):
        self.legend_ncolumns = self.n_plots//4+1
      if (not self.legend_customsize):
        if (self.n_plots < self.legend_ncolumns*4):
          self.legend_ylo = 0.9-0.03*(self.n_plots//self.legend_ncolumns+1)
      leg = ROOT.TLegend(self.legend_xlo,self.legend_ylo,self.legend_xhi,self.legend_yhi)
      leg.SetEntrySeparation(0)
      leg.SetTextSize(0.03)
      if (self.legend_ncolumns > -1):
        leg.SetNColumns(self.legend_ncolumns)
      else:
        n_columns = self.n_plots//4+1
        leg.SetNColumns(n_columns)
      color_index = 0
      for hist, color, style in zip(self.hists, self.hist_color, self.hist_style):
        hist.SetLineWidth(3)
        hist.SetLineColor(color)
        if style=='point':
          hist.SetMarkerColor(color)
          hist.SetMarkerStyle(ROOT.kFullCircle)
          hist.Draw('same P')
          leg.AddEntry(hist, hist.GetTitle(), 'LP')
        elif style=='outline':
          hist.Draw('same hist')
          leg.AddEntry(hist, hist.GetTitle(), 'F')
        elif style=='filled':
          hist.SetFillColor(color)
          hist.Draw('same hist')
          leg.AddEntry(hist, hist.GetTitle(), 'F')
      for graph, color in zip(self.graphs, self.graph_color):
        graph.SetLineWidth(3)
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(ROOT.kFullCircle)
        graph.Draw('same P')
        leg.AddEntry(graph, graph.GetTitle(), 'LP')
      leg.SetBorderSize(0)
      leg.Draw('same')
    else:
      ROOT.TGaxis.SetExponentOffset(-0.05,0.0,'y')
      self.hists[0].SetTitleSize(0.032,'y')
      self.hists[0].GetYaxis().SetTitle(self.y_title)
      self.hists[0].GetYaxis().SetNdivisions(606)
      self.hists[0].SetTitleSize(0.032,'x')
      self.hists[0].GetXaxis().SetTitle(self.x_title)
      self.hists[0].GetXaxis().SetNdivisions(606)
      self.hists[0].SetTitleSize(0.032,'z')
      self.hists[0].SetTitleOffset(2.0,'z')
      self.hists[0].GetZaxis().SetTitle(self.z_title)
      self.hists[0].Draw('colz')

    #draw CMS and lumi labels
    label = ROOT.TLatex()
    label.SetTextSize(0.032)
    label.SetNDC(ROOT.kTRUE)
    label.SetTextAlign(11)
    if (self.title_type == 'cms preliminary'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}')
    elif (self.title_type == 'cms work in progress'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Work in Progress}}')
    elif (self.title_type == 'cms supplementary'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Supplementary}}')
    elif (self.title_type == 'cms simulation'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}')
    elif (self.title_type == 'cms simulation supplementary'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Supplementary}}')
    elif (self.title_type == 'cms simulation preliminary'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Preliminary}}')
    elif (self.title_type == 'cms private work'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS} #scale[0.8]{#font[52]{Private Work}}')
    elif (self.title_type == 'cms'):
      label.DrawLatex(0.15,0.96,'#font[62]{CMS}')
    label.SetTextAlign(31)
    label.SetTextSize(0.03)
    lumi_energy_string = '#font[42]{'
    first = True
    for lumi_datum in self.lumi_data:
      if first:
        first = False
      else:
        lumi_energy_string += ' + '
      lumi_energy_string += str(lumi_datum[0])+' fb^{-1} ('+str(lumi_datum[1])+' TeV)'
    lumi_energy_string += '}'
    if not (self.title_type == 'cms simulation'):
      label.DrawLatex(1.0-right_margin-0.01,0.96,lumi_energy_string)
    top_pad.Modified()

    #draw lower plot
    if (self.plot_bottom):
      can.cd()
      bot_pad.Draw('same')
      bot_pad.cd()
      dummy_hist_lower.Draw()
      bottom_ref_value = 1.0
      if not self.bottom_is_ratio:
        bottom_ref_value = 0.0
      line = ROOT.TLine(self.x_min,bottom_ref_value,self.x_max,bottom_ref_value)
      line.SetNDC(ROOT.kFALSE)
      line.SetLineStyle(2)
      line.SetLineColor(ROOT.kBlack)
      line.SetLineWidth(2)
      line.Draw('SAME')
      first_bottom = True
      for bottom_plot, color in zip(self.bottom_plots, self.bottom_plot_color):
        if first_bottom:
          first_bottom = False
          bottom_plot.SetMarkerSize(0)
          bottom_plot.SetLineWidth(0)
          bottom_plot.SetFillColorAlpha(color,0.33)
          bottom_plot.Draw('same E2')
        else:
          bottom_plot.SetLineWidth(3)
          bottom_plot.SetLineColor(color)
          bottom_plot.SetMarkerColor(color)
          bottom_plot.SetMarkerStyle(ROOT.kFullCircle)
          bottom_plot.Draw('same P')
      bot_pad.Modified()

    #draw everything and save
    can.Draw()
    can.SaveAs(filename)
    return self

