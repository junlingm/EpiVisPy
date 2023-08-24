from setuptools import setup, find_packages

setup(
    name='Kerwin_RA',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        
        
        'pandas',
        'geopandas',
        'pyreadr',
        'matplotlib',
        'shapely.geometry',
        'mpl_toolkits.mplot3d',
        'mpl_toolkits.mplot3d.art3d',
        'numpy',
        'matplotlib.animation',
        'IPython.display'
        
        
        
        ],
    author='Jiangquan Dai',
    author_email='m15914863239@163.com',
    url='https://github.com/junlingm/EpiVisPy',
    description='A package to generate a 3D animation over a 2D map.',
    license='MIT',
    keywords="animation geospatial",
)

