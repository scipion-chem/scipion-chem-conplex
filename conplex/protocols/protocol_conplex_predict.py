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

   AI Generated:

      ProtConPLexPrediction - User Manual

      Overview
      --------
      The ProtConPLexPrediction protocol integrates the ConPLex deep-learning
      framework into Scipion-Chem for predicting protein–ligand interactions.
      It evaluates binding affinity by combining protein sequence information
      with ligand chemical representations, enabling fast, data-driven
      virtual screening.

      This protocol is designed to score large sets of ligands against one or
      multiple protein targets, providing an efficient alternative or complement
      to docking-based approaches.

      Input Requirements
      ------------------
      1. **Protein Sequences**:
         - Input must be provided as a SetOfSequences.
         - Each sequence represents a protein target for screening.

      2. **Ligands**:
         - Either:
           - A SetOfSmallMolecules, or
           - A SmallMoleculesLibrary (SMILES-based library).
         - Molecules must have valid chemical representations.

      Workflow
      --------
      1. **Ligand Preparation (if needed)**:
         - Conversion of molecular structures to SMILES format using OpenBabel.
         - Organization into a unified input directory.

      2. **Input Pairing**:
         - All protein–ligand combinations are generated.
         - Each pair is encoded into a tabular format for model input.

      3. **Model Selection**:
         - User selects a pre-trained ConPLex model.
         - Models may differ in architecture, training data, or prediction type.

      4. **Prediction Execution**:
         - ConPLex processes protein sequences and ligand SMILES.
         - Outputs interaction scores representing predicted binding affinity.

      Outputs
      -------
      - **Output Sequences (SetOfSequencesChem)**:
        - Each protein annotated with interaction scores for all ligands.
        - Scores stored as dictionaries and auxiliary files.

      - **Optional Molecule Outputs**:
        - If a single protein is provided:
          - Small molecules annotated with predicted scores, or
          - A ranked SmallMoleculesLibrary with scores.

      - **Interaction File**:
        - TSV file containing all protein–ligand predictions.

      Advanced Options
      ----------------
      - Selection of different trained ConPLex models.
      - Support for both explicit molecule sets and SMILES libraries.
      - Automatic handling of large combinatorial protein–ligand spaces.

      Validation & Warnings
      ---------------------
      - Predictions are data-driven and depend on model training scope.
      - Results may not generalize to proteins or chemotypes far from training data.
      - Large input sets may increase runtime due to combinatorial expansion.

      Practical Recommendations
      -------------------------
      - Use for rapid pre-screening before docking.
      - Combine with docking or rescoring protocols for improved reliability.
      - Validate model choice based on target domain when possible.
      - Limit input size for initial testing before scaling.

      Final Perspective
      -----------------
      ProtConPLexPrediction enables fast and scalable protein–ligand interaction
      prediction using deep learning. By leveraging sequence and chemical
      representations, it provides an efficient tool for virtual screening,
      rescoring, and prioritization workflows within Scipion-Chem.
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

