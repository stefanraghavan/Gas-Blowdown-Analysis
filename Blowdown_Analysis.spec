# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['Blowdown_Analysis.py'],
    pathex=[],
    binaries=[],
    datas=[('diagram.jpg', '.'), ('logo2.jpg', '.'), ('Theory & Assumptions Real Gas.pdf', '.'), ('gasProperties.xlsx', '.')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='Blowdown_Analysis',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['grid_icon.ico'],
)
