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
import json, os

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pwem.convert.atom_struct import AtomicStructHandler

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC
from pwchem.objects import SetOfSequencesChem, SmallMoleculesLibrary, SequenceChem

from .. import Plugin as conplexPlugin
from ..constants import CONPLEX_DIC

SEQ, AS, SEQS = 0, 1, 2

class ProtConPLexPrediction(EMProtocol):
  """Run a prediction using a ConPLex trained model over a set of proteins and ligands"""
  _label = 'conplex virtual screening'
  scoreName = 'ConPlex_score'
  stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                   label="Use GPU for execution: ",
                   help="This protocol has both CPU and GPU implementation.\
                                                 Select the one you want to use.")
    form.addHidden(params.GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                   help="Add a list of GPU devices that can be used")

    form.addSection(label='Input')
    iGroup = form.addGroup('Input Sequence')
    iGroup.addParam('inSeqForm', params.EnumParam, label='Input sequence(s) as: ', default=SEQS,
                    choices=['Sequence', 'AtomStruct', 'SetOfSequences'],
                    help='How to input the input sequence(s)')
    iGroup.addParam('inputSequence', params.PointerParam, pointerClass="Sequence", allowsNull=True,
                    label='Input protein sequence: ', condition=f'inSeqForm=={SEQ}',
                    help="Protein sequence to perform the screening on")
    iGroup.addParam('inputAS', params.PointerParam, pointerClass="AtomStruct", allowsNull=True,
                    label='Input protein structure: ', condition=f'inSeqForm=={AS}',
                    help="Protein structure to perform the screening on")
    iGroup.addParam('inChain', params.StringParam, label='Chain for input: ', condition=f'inSeqForm=={AS}',
                    help='Specify the protein chain to use as input')
    iGroup.addParam('inputSequences', params.PointerParam, pointerClass="SetOfSequences", allowsNull=True,
                    label='Input protein sequences: ', condition=f'inSeqForm=={SEQS}',
                    help="Set of protein sequences to perform the screening on")

    iGroup = form.addGroup('Input Ligands')
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
    args = f"--data-file {argFile} --model-path {modelPath} --outfile results.tsv "
    args += ' --device %(GPU)s'
    self.runJob(program, args, cwd=self._getPath())

  def createOutputStep(self):
      protSeqsDic = self.getInputSeqs()
      intDic, _, _ = self.parseInteractionsFile(self.getInteractionsFile())

      outSeqs = SetOfSequencesChem().create(outputPath=self._getPath())
      if self.inSeqForm.get() == SEQS:
          outSeqs.copyInfo(self.inputSequences.get())

      outputFile = self._getExtraPath("scoresFile.json")
      prevFile = outSeqs.getInteractScoresFile()
      if prevFile and os.path.exists(prevFile):
        shutil.copy(prevFile, outputFile)
      outSeqs.setInteractScoresFile(outputFile)

      data = {}
      for seqName, seq in protSeqsDic.items():
          outSeq = SequenceChem(name=seqName, sequence=seq)
          outSeqs.append(outSeq)

          if seqName not in data:
              data[seqName] = {}

          for mol, score in intDic[seqName].items():
              data[seqName][mol] = {self.scoreName: score}

      outSeqs.setInteractScoresDic(data)
      outSeqs.updateScoreTypes()

      outMols = self.inputLibrary.get() if self.useLibrary.get() else self.inputSmallMols.get()
      outSeqs.setInteractMols(mols=outMols)

      self._defineOutputs(outputSequences=outSeqs)

      if len(protSeqsDic) == 1:
          seqName = list(protSeqsDic.keys())[0]
          scoreDic = intDic[seqName]

          if self.useLibrary.get():
            inLib = self.inputLibrary.get()
            mapDic = inLib.getLibraryMap(inverted=True, fullLine=True)

            oLibFile = self._getPath('outputLibrary.smi')
            with open(oLibFile, 'w') as f:
              for smiName, score in scoreDic.items():
                f.write(f'{mapDic[smiName]}\t{score}\n')

            prevHeaders = inLib.getHeaders()
            outputLib = inLib.clone()
            outputLib.setFileName(oLibFile)
            outputLib.setHeaders(prevHeaders + ['Conplex_score'])
            self._defineOutputs(outputLibrary=outputLib)

          else:
              inSet = self.inputSmallMols.get()
              outputSet = inSet.createCopy(self._getPath(), copyInfo=True)

              for mol in inSet:
                  nMol = mol.clone()
                  molName = nMol.getMolName()

                  if molName in scoreDic:
                      score = scoreDic[molName]
                      setattr(nMol, self.scoreName, params.Float(score))
                      outputSet.append(nMol)

              outputSet.updateMolClass()
              self._defineOutputs(outputSmallMolecules=outputSet)


  ############## UTILS ########################
  def getDevices(self):
    if getattr(self, params.USE_GPU).get():
      gpuIdxs = getattr(self, params.GPU_LIST).get()
      if not gpuIdxs.strip():
        gpuIdxs = [0]
      else:
        gpuIdxs = [idx.strip() for idx in gpuIdxs.split(',')]
    else:
      gpuIdxs = ['cpu']
    return gpuIdxs

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
    if self.inSeqForm.get() == SEQ:
      seq = self.inputSequence.get()
      seqsDic[seq.getSeqName()] = seq.getSequence()
    elif self.inSeqForm.get() == AS:
      inAS = self.inputAS.get()
      seqName = os.path.basename(inAS.getFileName())
      handler = AtomicStructHandler(inAS.getFileName())
      struct = json.loads(getattr(self, 'inChain').get())  # From wizard dictionary
      chain_id, modelId = struct["chain"].upper().strip(), int(struct["model"])
      seq = str(handler.getSequenceFromChain(modelID=modelId, chainID=chain_id))

      seqsDic[seqName] = seq
    elif self.inSeqForm.get() == SEQS:
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

