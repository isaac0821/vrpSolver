  
import setuptools
from vrpSolver import __version__

long_description = """
	Working project: Learn VRP by coding
"""

setuptools.setup(
	name="vrpSolver",
	version=__version__,
	author="Lan Peng",
	author_email="lanpeng@buffalo.edu",
	description="VRP Solver",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/isaac0821/vrpSolver",
	license='MIT', 
	packages=['vrpSolver'], 
	download_url = "https://github.com/isaac0821/vrpSolver",
	project_urls={
		"Bug Tracker": "https://github.com/isaac0821/vrpSolver/issues",
		"Source Code": "https://github.com/isaac0821/vrpSolver",
	},
	python_requires='>=3',
	install_requires=[
		'numpy',  
		'scipy',
		'geopy',
		'matplotlib'
	],
	classifiers=[
		"Development Status :: 2 - Pre-Alpha",
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	]
)