import wheel
import setuptools
 
with open("README.rst", "r") as fh:
    long_description = fh.read()
  
setuptools.setup(
    name = 'bacteria',         
    packages = ['bacteria'],   
    version = '0.1.9',      
    license='MIT',       
    description = 'Super Segger Analysis in Python.',
    long_description_content_type="text/x-rst",
    url = 'https://bacteria.readthedocs.io',  
    download_url = 'https://github.com/tuliofalmeida/bacteria',    
    keywords = ['Data analysis', 'Cell analysis', 'Bacteria', 'SuperSegger'],   
    install_requires=[           
            'numpy',
            'pandas',
            'graphviz',
            'scikit-learn',
            'seaborn',
            'natsort',
            'anytree',
            'tqdm',
            # 'matlabengine',
            'scipy'
        ],
    classifiers=[
      'Development Status :: 4 - Beta',      
      'Intended Audience :: Developers',      
      'Topic :: Software Development :: Build Tools',
      'License :: OSI Approved :: MIT License',   
      'Programming Language :: Python :: 3.9',
    ],   
)