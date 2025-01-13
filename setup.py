from setuptools import setup, find_packages

setup(
    name='RedeTIL',
    version='0.1.0',
    author='Shengqiao Gao',
    author_email='292966713@qq.com',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'seaborn',
        'scanpy',
        'matplotlib',
        'plotly',
        'statsmodels',
    ],
    package_data={
        'RedeTIL': ['data/*.txt'],
    },
    include_package_data=True,
    description='A package for analyzing single-cell features that are associated with immune-checkpoint blockade therapy outcomes.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/gaoshengqiao/RedeTIL',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)
