# -*- mode: python -*-

# To stop recursion error
import sys
sys.setrecursionlimit(5000)

# Yet more astropy-wrangling
import astropy
astropy_base_dir = os.path.dirname(astropy.__file__)
astropy_tree = Tree(astropy_base_dir, prefix='astropy')

block_cipher = None

a = Analysis(['GUI.py'],
             pathex=['/home/chris/Dropbox/Work/Scripts/AstroCell/AstroCell'],
             binaries=[],
             datas=[],
             hiddenimports=['six','packaging','packaging.version','packaging.specifiers'], # To make sure that the six module is included (needed by astropy)
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          astropy_tree,
          exclude_binaries=True,
          name='GUI',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               astropy_tree,
               strip=False,
               upx=True,
               name='GUI')
