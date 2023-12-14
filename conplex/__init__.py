# **************************************************************************
# *
# * Authors:	Carlos Oscar Sorzano (coss@cnb.csic.es)
# *			 	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *			 	Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains protocols for creating and using ConPLex models for virtual screening
"""

# General imports
import os, subprocess, json

# Scipion em imports
import pwem
from scipion.install.funcs import InstallHelper

# Plugin imports
from pwchem import Plugin as pwchemPlugin
from .bibtex import _bibtexStr
from .constants import *

# Pluging variables
_logo = 'mit_logo.png'

class Plugin(pwchemPlugin):
	"""
	"""
	_cplHome = os.path.join(pwem.Config.EM_ROOT, CONPLEX_DIC['name'] + '-' + CONPLEX_DIC['version'])

	@classmethod
	def _defineVariables(cls):
		cls._defineEmVar(CONPLEX_DIC['home'], cls._cplHome)

	@classmethod
	def defineBinaries(cls, env):
		"""
        This function defines the binaries for each package.
        """
		cls.addConPLexPackage(env)

	@classmethod
	def addConPLexPackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing AutoDock. """
		# Instantiating the install helper
		installer = InstallHelper(CONPLEX_DIC['name'], packageHome=cls.getVar(CONPLEX_DIC['home']),
															packageVersion=CONPLEX_DIC['version'])

		# Installing package
		installer.getCloneCommand('https://github.com/samsledje/ConPLex.git', targeName='CONPLEX_CLONED') \
			.getCondaEnvCommand(pythonVersion='3.9', requirementsFile=False) \
			.addCommand(f'{cls.getEnvActivationCommand(CONPLEX_DIC)} && cd ConPLex && '
									f'pip install conplex-dti && '
									# f'conplex-dti download --to datasets --benchmarks {cls.getConPLexBenchmarks()} && '
									f'wget {cls.getConPLexModelUrl()} --no-check-certificate -P models', 'CONPLEX_INSTALLED') \
			.addPackage(env, ['git', 'conda', 'make', 'wget'], default=default)


	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def runScript(cls, protocol, scriptName, args, envDict, cwd=None, popen=False):
		""" Run rdkit command from a given protocol. """
		scriptName = cls.getScriptsDir(scriptName)
		fullProgram = '%s && %s %s' % (cls.getEnvActivationCommand(envDict), 'python', scriptName)
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

	@classmethod
	def getModelsDir(cls):
		return os.path.abspath(os.path.join(cls.getVar(CONPLEX_DIC['home']), 'ConPLex/models'))

	@classmethod
	def getLocalModels(cls):
		models = set([])
		for file in os.listdir(cls.getModelsDir()):
			models.add(file.strip())
		models = list(models)
		if models == []:
			models = ['None found']
		models.sort()
		return models
		
	# ---------------------------------- Utils functions-----------------------
	@classmethod
	def getConPLexModelUrl(cls):
		return 'https://cb.csail.mit.edu/cb/conplex/data/models/BindingDB_ExperimentalValidModel.pt'

	@classmethod
	def getConPLexBenchmarks(cls):
		return 'davis bindingdb biosnap biosnap_prot biosnap_mol dude'
