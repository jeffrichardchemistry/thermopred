from setuptools import setup, find_packages

with open("README.md", 'r') as fr:
	description = fr.read()

setup(
    name='Thermopred',
    version='1.0.0',
    url='https://github.com/jeffrichardchemistry/thermopred',
    license='GNU GPL',
    author='Jefferson Richard; Diullio P',
    author_email='jrichardquimica@gmail.com',
    keywords='Cheminformatics, Chemistry, QSAR, QSPR, Fingerprint, Spectroscopy',
    description='A python package to predict some thermochemical properties.',
    long_description = description,
    long_description_content_type = "text/markdown",
    packages=['Thermopred'],
    include_package_data=True,
    install_requires=['scikit-learn==1.3.2', 'pandas<=2.1.4', 'numpy<=1.26.4', 'rdkit==2023.9.6','xgboost==2.0.3'],#use this line in case of error in instalation
	classifiers = [
		'Intended Audience :: Developers',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Chemistry',
		'Topic :: Scientific/Engineering :: Physics',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Natural Language :: English',
		'Operating System :: Microsoft :: Windows',
		'Operating System :: POSIX :: Linux',
		'Environment :: MacOS X',
		'Programming Language :: Python :: 3.8',
		'Programming Language :: Python :: 3.9',
		'Programming Language :: Python :: 3.10']
)
