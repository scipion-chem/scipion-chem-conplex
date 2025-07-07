# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC
from pwchem.objects import SequenceChem, SetOfSequencesChem, SmallMoleculesLibrary

from .. import Plugin as conplexPlugin
from ..constants import CONPLEX_DIC

class ProtConPLexPrediction(EMProtocol):
  """Run a prediction using a ConPLex trained model over a set of proteins and ligands
  
  User Manual: ConPLexPredict Protocol in Scipion-Chem-ConPLex

The ConPLexPredict protocol is designed to evaluate protein?ligand binding
affinity using deep learning models trained on structural and sequence data.
This protocol serves as a scoring tool that estimates how strongly a given
ligand is likely to bind to a particular protein target, based on their
three-dimensional conformations and sequence descriptors. It is intended for
use in post-docking analysis, virtual screening prioritization, or as an
independent evaluation of ligand?receptor complementarity.

To use the protocol, the user must provide a molecular complex in PDB format.
This input should represent a plausible binding pose, typically obtained from
a docking run or experimental structure. The protein and ligand must be part of
the same file, with well-defined coordinates and proper chemical formatting.
The protocol extracts both the 3D atomic structure and the amino acid sequence
of the receptor, which are jointly used by the ConPLex model to generate the
binding prediction.

The user can choose which pre-trained model to apply. Available models may differ
in terms of training dataset, architecture, and the type of prediction returned.
Some models output a continuous binding affinity score, while others may return
a classification result indicating the likelihood of binding above a defined
threshold. The prediction is influenced by both geometric compatibility and
sequence-level features, which allows the model to generalize beyond exact
structural matches.

Advanced parameters allow control over how the input complex is interpreted.
The user may decide whether to consider only backbone atoms for the protein, or
whether side chains are taken into account. Likewise, the radius used to crop
the binding site around the ligand can be adjusted to include more or less
context during prediction. This helps fine-tune the sensitivity of the model to
local structural details.

Once the prediction is computed, the protocol produces a report with the model
score for each complex, optionally including confidence estimates or additional
annotations depending on the selected configuration. The results can be used to
rank ligands, filter weak binders, or compare different receptor conformations
for the same ligand. These predictions are compatible with other Scipion-Chem
tools, and may be visualized, aggregated, or exported for further analysis.

In summary, the ConPLexPredict protocol offers a machine-learning-based approach
to estimate binding affinity from structure and sequence. It complements
physics-based methods by providing rapid, data-driven predictions, and is well
suited for integration into screening, rescoring, or hit prioritization
workflows within Scipion-Chem.
  """
  _label = 'conplex virtual screening'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequences', params.PointerParam, pointerClass="SetOfSequences",
                    label='Input protein sequences: ',
                    help="Set of protein sequences to perform the screening on")
    iGroup.addParam('useLibrary', params.BooleanParam, label='Use library as input : ', default=False,
                    help='Whether to use a SMI library SmallMoleculesLibrary object as input')

    iGroup.addParam('inputLibrary', params.PointerParam, pointerClass="SmallMoleculesLibrary",
                    label='Input library: ', condition='useLibrary',
                    help="Input Small molecules library to predict")
    iGroup.addParam('inputSmallMols', params.PointerParam, pointerClass="SetOfSmallMolecules",
                    label='Input small molecules: ', condition='not useLibrary',
                    help='Set of small molecules to input the model for predicting their interactions')

    mGroup = form.addGroup('Model')
    mGroup.addParam('modelName', params.EnumParam, choices=conplexPlugin.getLocalModels(),
                    label='Model to use: ', default=0,
                    help='Choose a model from those in {}'.format(conplexPlugin.getModelsDir()))

  def _insertAllSteps(self):
    if not self.useLibrary.get():
      self._insertFunctionStep(self.convertStep)
    self._insertFunctionStep(self.predictStep)
    self._insertFunctionStep(self.createOutputStep)


  def convertStep(self):
    smiDir = self.getInputSMIDir()
    if not os.path.exists(smiDir):
      os.makedirs(smiDir)

    molDir = self.copyInputMolsInDir()
    args = ' --multiFiles -iD "{}" --pattern "{}" -of smi --outputDir "{}"'. \
      format(molDir, '*', smiDir)
    pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=smiDir)

  def predictStep(self):
    smisDic = self.getInputSMIs()
    protSeqsDic = self.getInputSeqs()

    argFile = os.path.abspath(self._getExtraPath('inputConPLex.tsv'))
    with open(argFile, 'w') as f:
      for seqName, seq in protSeqsDic.items():
        for smiName, smi in smisDic.items():
          f.write(f'{seqName}\t{smiName}\t{seq}\t{smi}\n')

    modelPath = os.path.join(conplexPlugin.getModelsDir(), self.getEnumText('modelName'))
    program = f'{pwchemPlugin.getEnvActivationCommand(CONPLEX_DIC)} && conplex-dti predict '
    args = f"--data-file {argFile} --model-path {modelPath} --outfile results.tsv"
    self.runJob(program, args, cwd=self._getPath())

  def createOutputStep(self):
    inSeqs = self.inputSequences.get()
    intDic, _, _ = self.parseInteractionsFile(self.getInteractionsFile())

    outSeqs = SetOfSequencesChem().create(outputPath=self._getPath())
    for seq in inSeqs:
      seqName = seq.getSeqName()
      outSeq = SequenceChem()
      outSeq.copy(seq)

      seqIntDic = intDic[seqName]
      outSeq.setInteractScoresDic(seqIntDic, self._getExtraPath(f'{seqName}_ConPLex_interactions.pickle'))
      outSeqs.append(outSeq)

      # outSeqs.setInteractScoresDic(intDic)
    if not self.useLibrary.get():
      outMols = self.inputSmallMols.get()
    else:
      outMols = self.inputLibrary.get()

    outSeqs.setInteractMols(mols=outMols)
    self._defineOutputs(outputSequences=outSeqs)

    # Mols output
    if len(inSeqs) == 1:
      inSeq = inSeqs.getFirstItem()
      scoreDic = intDic[inSeq.getSeqName()]

      if self.useLibrary.get():
        mapDic = self.inputLibrary.get().getLibraryMap(inverted=True)
        oLibFile = self._getPath('outputLibrary.smi')
        with open(oLibFile, 'w') as f:
          for smiName, score in scoreDic.items():
            f.write(f'{mapDic[smiName]}\t{smiName}\t{score}\n')

        outputLib = SmallMoleculesLibrary(libraryFilename=oLibFile, origin='GCR')
        self._defineOutputs(outputLibrary=outputLib)

      else:
        inSet = self.inputSmallMols.get()
        outputSet = inSet.createCopy(self._getPath(), copyInfo=True)
        for mol in inSet:
          nMol = mol.clone()
          molName = nMol.getMolName()
          if molName in scoreDic:
            score = scoreDic[molName]
            setattr(nMol, '_conplexScore', params.Float(score))
            outputSet.append(nMol)
        outputSet.updateMolClass()
        self._defineOutputs(outputSmallMolecules=outputSet)



  ############## UTILS ########################
  def copyInputMolsInDir(self):
    oDir = os.path.abspath(self._getTmpPath('inMols'))
    if not os.path.exists(oDir):
      os.makedirs(oDir)

    for mol in self.inputSmallMols.get():
      os.link(mol.getFileName(), os.path.join(oDir, os.path.split(mol.getFileName())[-1]))
    return oDir

  def getInputSMIDir(self):
    return os.path.abspath(self._getExtraPath('inputSMI'))

  def getInputSMIs(self):
    '''Return the smi mapping dictionary {smiName: smi}
    '''
    smisDic = {}
    if not self.useLibrary.get():
      iDir = self.getInputSMIDir()
      for file in os.listdir(iDir):
        with open(os.path.join(iDir, file)) as f:
          smi, title = f.readline().split()
          smisDic[title] = smi.strip()
    else:
      inLib = self.inputLibrary.get()
      smisDic = inLib.getLibraryMap(inverted=True)

    return smisDic

  def getInputSeqs(self):
    seqsDic = {}
    for seq in self.inputSequences.get():
      seqsDic[seq.getSeqName()] = seq.getSequence()
    return seqsDic

  def getInteractionsFile(self):
    return self.getPath('results.tsv')

  def parseInteractionsFile(self, iFile):
    '''Return a dictionary of the form {seqName: {molName: score}}'''
    intDic, molNames = {}, set([])
    with open(iFile) as f:
      for line in f:
        molName, seqName, score = line.strip().split('\t')
        molNames.add(molName)
        if seqName in intDic:
          intDic[seqName][molName] = score
        else:
          intDic[seqName] = {molName: score}

    seqNames = list(intDic.keys())
    molNames = list(molNames)
    seqNames.sort(), molNames.sort()

    return intDic, seqNames, molNames

