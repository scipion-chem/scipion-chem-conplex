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
from tkinter.messagebox import askokcancel

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol  import params

from pwem.protocols import ProtSubSet

from ..protocols import ProtConPLexPrediction, ProtExtractInteractingMols


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
    intDic, seqNames, molNames = self.getProtocol().parseInteractionsFile(self.getProtocol().getInteractionsFile())

    form.addSection(label='ConPLex viewer')
    sGroup = form.addGroup('Score filter')
    sGroup.addParam('chooseSeq', params.EnumParam, label='Display results for protein: ',
                    choices=['All'] + seqNames, default=0,
                    help='Display the selected results only for the specific protein')
    sGroup.addParam('chooseMol', params.EnumParam, label='Display results for molecules: ',
                    choices=['All'] + molNames, default=0,
                    help='Display the selected results only for the specific molecule')
    sGroup.addParam('scoreThres', params.FloatParam, label='Score threshold: ', default=0,
                    help='Display the interaction results over the selected score threshold')

    hGroup = form.addGroup('HeatMap View')
    hGroup.addParam('displayHeatMap', params.LabelParam, label='Display DTI heatmap: ',
                    help='Display a heatmap showing the scores for each protein-molecule pair')

    oGroup = form.addGroup('Generate output')
    oGroup.addParam('genProts', params.LabelParam, label='Generate protein sequences: ',
                    help='Generate an output of the filtered sequences of proteins'
                         '\ne.g: proteins passing the score filter specified above')
    oGroup.addParam('genMols', params.LabelParam, label='Generate small molecules: ',
                    help='Generate an output of the filtered molecules. '
                         '\ne.g: molecules passing the score filter specified above')

  def _getVisualizeDict(self):
    return {
      'displayHeatMap': self._viewHeatMap,

      'genProts': self._generateProts,
      'genMols': self._generateMols,
    }


################# DOCKING VIEWS ###################
  def _viewHeatMap(self, paramName=None):
    outSeqs = self.getOutSeqs()
    seqNames, molNames = outSeqs.getSequenceNames(), outSeqs.getInteractMolNames()
    if len(seqNames) * len(molNames) > 1000:
      try:
        msg = f"A heatmap with the scores of the selected sequence-molecule pairs will be displayed.\n" \
              f"Big sets of sequences/molecules might take a while to process and display.\n" \
              f"Do you want to continue?"
        answer = askokcancel("Display heatmap", msg, parent=None)
      except:
        answer = True
    else:
      answer = True

    if answer:
      intAr, seqNames, molNames = self.getFilteredOutput()

      fig, ax = plt.subplots()
      im, cbar = heatmap(intAr, seqNames, molNames, ax=ax,
                         cmap="YlGn", cbarlabel="ConPLex interaction score")
      texts = annotate_heatmap(im, valfmt="{x:.2f}")
      fig.tight_layout()
      plt.show()

  def _generateProts(self, paramName=None):
    prot = self.getProtocol()
    project = prot.getProject()
    _, seqNames, _ = self.getFilteredOutput()

    objIds = []
    for seq in prot.outputSequences:
      if seq.getSeqName() in seqNames:
        objIds.append(str(seq.getObjId()))

    try:
      msg = f"An output of the following {len(seqNames)} sequences will be generated:\n{','.join(seqNames)}"
      answer = askokcancel("Generate sequences output",  msg,  parent=None)
    except:
      answer = True

    if answer:
      protFilter = project.newProtocol(
        ProtSubSet,
        selectIds=True, range=','.join(objIds)
      )

      protFilter.setObjLabel('Filtered sequences')
      protFilter.inputFullSet.set(prot)
      protFilter.inputFullSet.setExtended('outputSequences')

      project.launchProtocol(protFilter, wait=True)

  def _generateMols(self, paramName=None):
    prot = self.getProtocol()
    project = prot.getProject()

    _, _, filtMolNames = self.getFilteredOutput()

    try:
      msg = f"An output of the following {len(filtMolNames)} molecules will be generated:\n{','.join(filtMolNames)}"
      answer = askokcancel("Generate sequences output",  msg,  parent=None)
    except:
      answer = True

    if answer:
      seqNames, molNames = self.getEnumText('chooseSeq'), self.getEnumText('chooseMol')
      scThres = self.scoreThres.get()
      protFilter = project.newProtocol(
        ProtExtractInteractingMols,
        chooseSeq=seqNames, chooseMol=molNames,
        scThres=scThres
      )
      protFilter.inputSequences.set(prot)
      protFilter.inputSequences.setExtended('outputSequences')

      project.launchProtocol(protFilter, wait=True)

#####################################################
  def getProtocol(self):
    return self.protocol

  def getOutSeqs(self):
    return self.protocol.outputSequences

  def getFilteredOutput(self):
    outSeqs = self.getOutSeqs()
    intDic = outSeqs.getInteractScoresDic()

    seqNames, molNames = outSeqs.getSequenceNames(), outSeqs.getInteractMolNames()
    seqNames, molNames = self.filterIds(seqNames, molNames)

    intAr = self.formatInteractionsArray(intDic, seqNames, molNames)
    intAr, seqNames, molNames = self.filterScores(intAr, seqNames, molNames)
    return intAr, seqNames, molNames

  def filterIds(self, seqNames, molNames):
    if self.getEnumText('chooseSeq') != 'All':
      seqNames = [seqName for seqName in seqNames if seqName == self.getEnumText('chooseSeq')]

    if self.getEnumText('chooseMol') != 'All':
      molNames = [molName for molName in molNames if molName == self.getEnumText('chooseMol')]

    return seqNames, molNames

  def filterScores(self, intAr, seqNames, molNames):
    ips, ims = [], []
    scThres = self.scoreThres.get()

    for ip, seqName in enumerate(seqNames):
      if any(intAr[ip,:] > scThres):
        ips.append(ip)

    for im, molName in enumerate(molNames):
      if any(intAr[:, im] > scThres):
        ims.append(im)

    if not len(seqNames) == len(ips):
      seqNames = list(np.array(seqNames)[ips])
      intAr = intAr[ips, :]

    if not len(molNames) == len(ims):
      molNames = list(np.array(molNames)[ims])
      intAr = intAr[:, ims]

    return intAr, seqNames, molNames


  def formatInteractionsArray(self, intDic, seqNames, molNames):
    intAr = np.zeros((len(seqNames), len(molNames)))
    for i, seqName in enumerate(seqNames):
      for j, molName in enumerate(molNames):
        intAr[i, j] = intDic[seqName][molName]
    return intAr





