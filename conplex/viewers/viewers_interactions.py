# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol  import params

import pwem.viewers.views as views
import pwem.viewers.showj as showj

from ..protocols import ProtConPLexPrediction


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
  """
  Create a heatmap from a numpy array and two lists of labels.

  Parameters
  ----------
  data
      A 2D numpy array of shape (M, N).
  row_labels
      A list or array of length M with the labels for the rows.
  col_labels
      A list or array of length N with the labels for the columns.
  ax
      A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
      not provided, use current axes or create a new one.  Optional.
  cbar_kw
      A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
  cbarlabel
      The label for the colorbar.  Optional.
  **kwargs
      All other arguments are forwarded to `imshow`.
  """

  if ax is None:
    ax = plt.gca()

  if cbar_kw is None:
    cbar_kw = {}

  # Plot the heatmap
  im = ax.imshow(data, **kwargs)

  # Create colorbar
  cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
  cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

  # Show all ticks and label them with the respective list entries.
  ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
  ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

  # Let the horizontal axes labeling appear on top.
  ax.tick_params(top=True, bottom=False,
                 labeltop=True, labelbottom=False)

  # Rotate the tick labels and set their alignment.
  plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
           rotation_mode="anchor")

  # Turn spines off and create white grid.
  ax.spines[:].set_visible(False)

  ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
  ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
  ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
  ax.tick_params(which="minor", bottom=False, left=False)

  return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
  """
  A function to annotate a heatmap.

  Parameters
  ----------
  im
      The AxesImage to be labeled.
  data
      Data used to annotate.  If None, the image's data is used.  Optional.
  valfmt
      The format of the annotations inside the heatmap.  This should either
      use the string format method, e.g. "$ {x:.2f}", or be a
      `matplotlib.ticker.Formatter`.  Optional.
  textcolors
      A pair of colors.  The first is used for values below a threshold,
      the second for those above.  Optional.
  threshold
      Value in data units according to which the colors from textcolors are
      applied.  If None (the default) uses the middle of the colormap as
      separation.  Optional.
  **kwargs
      All other arguments are forwarded to each call to `text` used to create
      the text labels.
  """

  if not isinstance(data, (list, np.ndarray)):
    data = im.get_array()

  # Normalize the threshold to the images color range.
  if threshold is not None:
    threshold = im.norm(threshold)
  else:
    threshold = im.norm(data.max()) / 2.

  # Set default alignment to center, but allow it to be
  # overwritten by textkw.
  kw = dict(horizontalalignment="center",
            verticalalignment="center")
  kw.update(textkw)

  # Get the formatter in case a string is supplied
  if isinstance(valfmt, str):
    valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

  # Loop over the data and create a `Text` for each "pixel".
  # Change the text's color depending on the data.
  texts = []
  for i in range(data.shape[0]):
    for j in range(data.shape[1]):
      kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
      text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
      texts.append(text)

  return texts

class ConPLexViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer ConPLex prediction'
  _targets = [ProtConPLexPrediction]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    intDic, protIds, molIds = self.parseInteractionsFile(self.getInteractionsFile())

    form.addSection(label='ConPLex viewer')
    sGroup = form.addGroup('Score filter')
    sGroup.addParam('chooseProt', params.EnumParam, label='Display results for protein: ',
                    choices=['All'] + protIds, default=0,
                    help='Display the selected results only for the specific protein')
    sGroup.addParam('chooseMol', params.EnumParam, label='Display results for molecules: ',
                    choices=['All'] + molIds, default=0,
                    help='Display the selected results only for the specific molecule')
    sGroup.addParam('scoreThres', params.FloatParam, label='Score threshold: ', default=0,
                    help='Display the interaction results over the selected score threshold')

    hGroup = form.addGroup('HeatMap View')
    hGroup.addParam('displayHeatMap', params.LabelParam, label='Display DTI heatmap: ',
                    help='Display a heatmap showing the scores for each protein-molecule pair')

  def _getVisualizeDict(self):
    return {
      'displayHeatMap': self._viewHeatMap,
    }


################# DOCKING VIEWS ###################
  def _viewHeatMap(self, paramName=None):
      intDic, protIds, molIds = self.parseInteractionsFile(self.getInteractionsFile())
      protIds, molIds = self.filterIds(protIds, molIds)
      intAr = self.formatInteractionsArray(intDic, protIds, molIds)
      intAr, protIds, molIds = self.filterScores(intAr, protIds, molIds)

      fig, ax = plt.subplots()
      im, cbar = heatmap(intAr, protIds, molIds, ax=ax,
                         cmap="YlGn", cbarlabel="ConPLex interaction score")
      texts = annotate_heatmap(im, valfmt="{x:.2f}")
      fig.tight_layout()
      plt.show()

#####################################################
  def getProtocol(self):
    return self.protocol

  def filterIds(self, protIds, molIds):
    if self.getEnumText('chooseProt') != 'All':
      protIds = [protId for protId in protIds if protId == self.getEnumText('chooseProt')]

    if self.getEnumText('chooseMol') != 'All':
      molIds = [molId for molId in molIds if molId == self.getEnumText('chooseMol')]

    return protIds, molIds

  def filterScores(self, intAr, protIds, molIds):
    ips, ims = [], []
    scThres = self.scoreThres.get()

    for ip, protId in enumerate(protIds):
      if any(intAr[ip,:] > scThres):
        ips.append(ip)

    for im, molId in enumerate(molIds):
      if any(intAr[:, im] > scThres):
        ims.append(im)

    if not len(protIds) == len(ips):
      protIds = list(np.array(protIds)[ips])
      intAr = intAr[ips, :]

    if not len(molIds) == len(ims):
      molIds = list(np.array(molIds)[ims])
      intAr = intAr[:, ims]

    return intAr, protIds, molIds

  def getInteractionsFile(self):
    return self.getProtocol().getPath('results.tsv')

  def formatInteractionsArray(self, intDic, protIds, molIds):
    intAr = np.zeros((len(protIds), len(molIds)))
    for i, protId in enumerate(protIds):
      for j, molId in enumerate(molIds):
        intAr[i, j] = intDic[protId][molId]
    return intAr


  def parseInteractionsFile(self, iFile):
    '''Return a dictionary of the form {protId: {molId: score}}'''
    intDic, molIds = {}, set([])
    with open(iFile) as f:
      for line in f:
        molId, protId, score = line.split('\t')
        molIds.add(molId)
        if protId in intDic:
          intDic[protId][molId] = score
        else:
          intDic[protId] = {molId: score}

    protIds = list(intDic.keys())
    molIds = list(molIds)
    protIds.sort(), molIds.sort()

    return intDic, protIds, molIds


