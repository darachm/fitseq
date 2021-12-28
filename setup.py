from setuptools import setup

with open('README.md','r') as fh:
    long_description = fh.read()

setup(
    name='fitseq',
    version='v1.5.0',
    author='Fangfei Li, Darach Miller',
    url='https://github.com/darachm/fitseq',
    description=(
        'A utility for fitting lineage fitnesses within a pooled competition '
        'experiment, achieved by iteratively optimizing models of '
        'individual and population-average lineage fitness.'
        ),
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
        ],
    python_requires=">=3.6",
    packages=['fitseq'],
    scripts=['fitseq/evo_simulator.py','fitseq/fitseq.py'],
    install_requires=['numpy>=1.17.3','pandas>=0.25.3','scipy>=1.3.1','tqdm'],
    license='MIT'
    )
