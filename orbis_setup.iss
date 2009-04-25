; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{70C1BF62-DE5E-4402-A317-C8043AB2FDE4}
AppName=Orbis
AppVerName=Orbis 0.1.1
AppPublisher=Randle Taylor
AppPublisherURL=http://www.simplehuckel.com
AppSupportURL=http://www.simplehuckel.com
AppUpdatesURL=http://www.simplehuckel.com
DefaultDirName={pf}\Orbis
DefaultGroupName=Orbis
AllowNoIcons=yes
LicenseFile=C:\Python25\code\huckel\dist\license.txt
OutputBaseFilename=install_orbis_0.1.1
Compression=lzma
SolidCompression=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked
Name: "quicklaunchicon"; Description: "{cm:CreateQuickLaunchIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "C:\Python25\code\huckel\dist\orbis.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "C:\Python25\code\huckel\dist\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Icons]
Name: "{group}\Orbis"; Filename: "{app}\orbis.exe"
Name: "{group}\{cm:UninstallProgram,Orbis}"; Filename: "{uninstallexe}"
Name: "{commondesktop}\Orbis"; Filename: "{app}\orbis.exe"; Tasks: desktopicon
Name: "{userappdata}\Microsoft\Internet Explorer\Quick Launch\Orbis"; Filename: "{app}\orbis.exe"; Tasks: quicklaunchicon

[Run]
Filename: "{app}\orbis.exe"; Description: "{cm:LaunchProgram,Orbis}"; Flags: nowait postinstall skipifsilent
