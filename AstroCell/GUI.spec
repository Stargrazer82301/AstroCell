# -*- mode: python -*-

# To stop recursion error
import sys
sys.setrecursionlimit(5000)

# To make sure that the six module is included (needed by astropy)
hiddenimports=['six','packaging','packaging.version','packaging.specifiers']

block_cipher = None

a = Analysis(['GUI.py'],
             pathex=['/home/chris/Dropbox/Work/Scripts/AstroCell/AstroCell'],
             binaries=[],
             datas=[],
             hiddenimports=[],
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
               strip=False,
               upx=True,
               name='GUI')
