from setuptools import setup

with open('README.md','r') as fh:
    long_description = fh.read()

setup(
    name='fitseq',
    version='1.2.0-rc1',
    description=(
        'A utility for fitting lineage fitnesses within a pooled competition '
        'experiment, achieved by iteratively optimizing models of '
        'individual and population-average lineage fitness.'
        ),
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=['fitseq'],
    scripts=['fitseq/evo_simulator.py','fitseq/fitseq.py'],
    install_requires=['numpy>=1.17.3','pandas>=0.25.3','scipy>=1.3.1','tqdm'],
    license='MIT',
    long_description=long_description,
    long_description_content_type='text/markdown'
    )
