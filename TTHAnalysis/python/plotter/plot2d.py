import CMGTools.TTHAnalysis.plotter.CMS_lumi as CMS_lumi
from optparse import OptionParser
import ROOT, sys, array
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)    

parser = OptionParser(usage="%prog rootfile histogram")
parser.add_option("-e"    , "--energy", dest="energy"  , type="float", default=13, help="Energy to plot");
parser.add_option("-l"    , "--lumi", dest="lumi"  , type="float", default=35.9, help="Luminosity to plot");
parser.add_option("-x"    , dest="xlabel", type="string", default=None        , help="Label of the x axis");
parser.add_option("-y"    , dest="ylabel", type="string", default=None        , help="Label of the y axis");
parser.add_option("-z"    , dest="zlabel", type="string", default="fake ratio", help="Label of the z axis");
parser.add_option("--exts", dest="exts"  , type="string", action="append", default=["png","pdf","C","root"], help="Output types");
parser.add_option("--lspam", dest="lspam"  , type="string", default="Preliminary", help="Stuff to write onto the plot (Preliminary, Simulation, Internal)");
parser.add_option("--plotmode", dest="plotmode", type="string", default="colz text e", help="Option to use when drawing the plot");
parser.add_option("--zrange"  , dest="zrange"  , type="float" , nargs=2, default=[0,1], help="Range of the z axis");
parser.add_option("--digits"  , dest="digits"  , type="float" , default=1.3, help="Format of digits for text on plots");
parser.add_option("--logx"    , dest="logx"    , action="store_true", default=False, help="Logx");
parser.add_option("--logy"    , dest="logy"    , action="store_true", default=False, help="Logy");
parser.add_option("--logz"    , dest="logz"    , action="store_true", default=False, help="Logz");
(options, args) = parser.parse_args()


ROOT.gStyle.SetPaintTextFormat("%sf"%options.digits)    
canvas = ROOT.TCanvas("c", "c", 600, 600)

f = ROOT.TFile.Open(args[0],"read")
if not f: sys.exit()

histo = f.Get(args[1])
if not histo: sys.exit()

canvas.cd()

histo.Draw()
ROOT.gPad.Update()

canvas.SetTickx(1)
canvas.SetTicky(1)
canvas.SetRightMargin (0.19)
canvas.SetTopMargin   (0.06)
canvas.SetLeftMargin  (0.14)
canvas.SetBottomMargin(0.14)


if options.xlabel: histo.GetXaxis().SetTitle(options.xlabel)
if options.ylabel: histo.GetYaxis().SetTitle(options.ylabel)
if options.zlabel: histo.GetZaxis().SetTitle(options.zlabel)

histo.GetXaxis().SetNdivisions(505)
histo.GetXaxis().SetLabelFont(42)
histo.GetXaxis().SetLabelSize(0.035)
histo.GetXaxis().SetTitleFont(42)
histo.GetXaxis().SetTitleSize(0.05)
histo.GetXaxis().SetLabelOffset(0.007)
histo.GetXaxis().SetTitleOffset(1.2)

histo.GetYaxis().SetNdivisions(505)
histo.GetYaxis().SetLabelFont(42)
histo.GetYaxis().SetLabelSize(0.035)
histo.GetYaxis().SetTitleFont(42)
histo.GetYaxis().SetTitleSize(0.05)
histo.GetYaxis().SetLabelOffset(0.007)
histo.GetYaxis().SetTitleOffset(1.30)

histo.GetZaxis().SetLabelFont(42)
histo.GetZaxis().SetTitleFont(42)
histo.GetZaxis().SetLabelSize(0.035)
histo.GetZaxis().SetTitleSize(0.035)
histo.GetZaxis().SetTitleOffset(1.7)
histo.SetMinimum(options.zrange[0])
histo.SetMaximum(options.zrange[1])
NRGBs = 5
NCont = 255
stops = array.array("d",[0.00, 0.34, 0.61, 0.84, 1.00])
red   = array.array("d",[0.50, 0.50, 1.00, 1.00, 1.00])
green = array.array("d",[0.50, 1.00, 1.00, 0.60, 0.50])
blue  = array.array("d",[1.00, 1.00, 0.50, 0.40, 0.50])
ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
ROOT.gStyle.SetNumberContours(NCont)
canvas.cd()
histo.Draw("colz")
ROOT.gPad.Update()
palette = histo.GetListOfFunctions().FindObject("palette")
palette.SetX1NDC(1.-0.18)
palette.SetY1NDC(0.14)
palette.SetX2NDC(1.-0.13)
palette.SetY2NDC(1.-0.06)
palette.SetLabelFont(42)
palette.SetLabelSize(0.035)

histo.Draw(options.plotmode)

CMS_lumi.writeExtraText = 1
CMS_lumi.lumi_13TeV = "%.1f fb^{-1}" % options.lumi
CMS_lumi.extraText  = options.lspam
CMS_lumi.lumi_sqrtS = options.energy
CMS_lumi.CMS_lumi(canvas, 4, 0, 0.05)

if options.logx: canvas.SetLogx()
if options.logy: canvas.SetLogy()
if options.logz: canvas.SetLogz()

for ext in options.exts:
	canvas.SaveAs(args[1]+"."+ext)


