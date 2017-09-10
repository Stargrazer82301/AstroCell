# Import smorgasbord
import setuptools

# Configure setup
setuptools.setup(name = 'astrocell',
                 version = '0.1',
                 description = 'Package to locate, count, and classify cells in microscopy; based on astronomical source-extaction techniques.',
                 url = 'https://github.com/Stargrazer82301/AstroCell',
                 author = 'Chris Clark (github.com/Stargrazer82301)',
                 author_email = 'cjrc88@gmail.com',
                 license = 'MIT',
                 classifiers = ['Programming Language :: Python :: 3.4',
                                'License :: OSI Approved :: MIT License',
                                'Intended Audience :: Science/Research',
                                'Development Status :: 3 - Alpha'],
                 packages = setuptools.find_packages(),
			  #install_requires = ['chrisfuncs==0.1','PIL>=4.2.1','scikit-image>=0.13.0','scikit-learn>=0.17.1','photutils>=0.3.2','joblib>=0.11','webcolors','pyamg>=3.3.0'],
			  #dependency_links=['https://github.com/Stargrazer82301/ChrisFuncs.git#egg=chrisfuncs-1.0'],
                 zip_safe = False)