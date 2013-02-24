"""Setup file to generate a distribution of njord

usage:    python setup.py sdist
          python setup.py install
"""


from distutils.core import setup

setup(name = 'satnrt',
      version = '0.5',
      description = 'Download and merge NRT satellite fields from NASA',
      long_description = "README.md",
      #long_description=open('docs/README.rst', 'rt').read()

      author = 'Bror Jonsson',
      author_email = 'brorfred@gmail.com',
      url = 'http://github.com/brorfred/neartime',
      requires = ["numpy(>=1.5)",
                  "matplotlib(>=1.1.0)",
                  "mpl_toolkits(>=1.0)",
                  "requests(>=1.1.0)",
                  "projmap(>=0.5)"],
      packages = ['satnrt'],
      scripts = ["slippy.py",],
      #package_data = {'njord': ['njord.cfg']},
     )
