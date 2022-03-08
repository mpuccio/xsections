import json
from array import array
import numpy as np
from mendeleev import element
from ROOT import TGraphErrors, TFile
import re

regex = r"[0-9]*"
inFile = open("inel.json")
data = json.load(inFile)

targets = [ds['reaction']['Target'] for ds in data['datasets']]
xsections = {}
labels = ['EN', 'Data', 'dData']
for item in targets:
  xsections[item] = {}
  for lab in labels:
    xsections[item][lab] = []

for ds in data['datasets']:
  target = ds['reaction']['Target']
  for point in ds['dataPoints']:
    skip = False
    for lab in labels:
      if lab not in point:
        skip = True
    if skip:
      break
    for lab in labels:
      if lab in point:
        xsections[target][lab].append(point[lab])
print(xsections)

outFile = TFile("xsections.root","recreate")
for lab, xs in xsections.items():
  if len(xs['Data'])==0:
    continue
  x = array('f',[ x / 1000 for x in xs['EN']])
  ex = array('f', [0 for _ in xs['EN']])
  y = array('f', xs['Data'])
  ey = array('f', xs['dData'])
  graph = TGraphErrors(len(xs['Data']), x, y, ex, ey)
  graph.SetMarkerStyle(20)
  matches = list(filter(None, re.findall(regex, lab)))
  print(matches)
  if int(matches[1]) == 0:
    abund = 0
    el = element(int(matches[0]))
    for iso in el.isotopes:
      if iso.abundance is not None and iso.abundance > abund:
        abund = iso.abundance
        massN = iso.mass_number
    print(f"missing mass number, using the one from mendelev for element {el.name}: {massN}")
    graph.Write(f'{matches[0]}-{massN}')
  else:
    graph.Write(f'{matches[0]}-{matches[1]}')