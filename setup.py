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
                 zip_safe = False)