; Beavyr installer for Windows.
;
; Built by .github/workflows/release.yml on a Windows machine, since NSIS and
; the compiler both live there. Produces a setup.exe that installs into
; Program Files, puts Beavyr in the Start menu, and registers an uninstaller so
; it appears in "Apps & features" like any other program.
;
; VERSION is passed in with /DVERSION=x.y.z

!include "MUI2.nsh"
!include "FileFunc.nsh"

!ifndef VERSION
  !define VERSION "0.0.0"
!endif

Name "Beavyr"
OutFile "beavyr-${VERSION}-setup.exe"
Unicode true
InstallDir "$PROGRAMFILES64\Beavyr"
InstallDirRegKey HKLM "Software\Beavyr" "InstallDir"
RequestExecutionLevel admin

VIProductVersion "${VERSION}.0"
VIAddVersionKey "ProductName"     "Beavyr"
VIAddVersionKey "FileDescription" "Molecular viewer and quantum-chemistry front-end"
VIAddVersionKey "FileVersion"     "${VERSION}"
VIAddVersionKey "ProductVersion"  "${VERSION}"
VIAddVersionKey "LegalCopyright"  "GPL-3.0-or-later"

!define MUI_ABORTWARNING
!define MUI_ICON "beavyr.ico"
!define MUI_UNICON "beavyr.ico"

!insertmacro MUI_PAGE_LICENSE "LICENSE.txt"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!define MUI_FINISHPAGE_RUN "$INSTDIR\beavyr.exe"
!define MUI_FINISHPAGE_RUN_TEXT "Start Beavyr"
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

!insertmacro MUI_LANGUAGE "English"

Section "Beavyr" SecMain
  SetOutPath "$INSTDIR"
  File "beavyr.exe"
  File "beavyr.ico"
  File "LICENSE.txt"
  File "README.md"

  WriteRegStr HKLM "Software\Beavyr" "InstallDir" "$INSTDIR"

  CreateDirectory "$SMPROGRAMS\Beavyr"
  CreateShortcut "$SMPROGRAMS\Beavyr\Beavyr.lnk" "$INSTDIR\beavyr.exe" "" "$INSTDIR\beavyr.ico"
  CreateShortcut "$SMPROGRAMS\Beavyr\Uninstall Beavyr.lnk" "$INSTDIR\uninstall.exe"

  WriteUninstaller "$INSTDIR\uninstall.exe"

  ; This block is what puts Beavyr in "Apps & features", so it can be removed
  ; the ordinary way rather than by deleting a folder.
  !define UNINST_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\Beavyr"
  WriteRegStr   HKLM "${UNINST_KEY}" "DisplayName"     "Beavyr"
  WriteRegStr   HKLM "${UNINST_KEY}" "DisplayVersion"  "${VERSION}"
  WriteRegStr   HKLM "${UNINST_KEY}" "DisplayIcon"     "$INSTDIR\beavyr.ico"
  WriteRegStr   HKLM "${UNINST_KEY}" "Publisher"       "Peter Szabo, KU Leuven"
  WriteRegStr   HKLM "${UNINST_KEY}" "UninstallString" "$\"$INSTDIR\uninstall.exe$\""
  WriteRegStr   HKLM "${UNINST_KEY}" "InstallLocation" "$INSTDIR"
  WriteRegDWORD HKLM "${UNINST_KEY}" "NoModify" 1
  WriteRegDWORD HKLM "${UNINST_KEY}" "NoRepair" 1

  ${GetSize} "$INSTDIR" "/S=0K" $0 $1 $2
  WriteRegDWORD HKLM "${UNINST_KEY}" "EstimatedSize" "$0"
SectionEnd

Section "Uninstall"
  Delete "$INSTDIR\beavyr.exe"
  Delete "$INSTDIR\beavyr.ico"
  Delete "$INSTDIR\LICENSE.txt"
  Delete "$INSTDIR\README.md"
  Delete "$INSTDIR\uninstall.exe"
  RMDir  "$INSTDIR"

  Delete "$SMPROGRAMS\Beavyr\Beavyr.lnk"
  Delete "$SMPROGRAMS\Beavyr\Uninstall Beavyr.lnk"
  RMDir  "$SMPROGRAMS\Beavyr"

  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Beavyr"
  DeleteRegKey HKLM "Software\Beavyr"
SectionEnd
