# リリース手順

- Visual Studio 2019のソリューションエクスプローラーで「SELDLA」を右クリックし、「発行」からWin, Mac, Linux用をそれぞれ発行する。
- publish先のSELDLA/bin/Release/netcoreapp3.1/publishフォルダを開き、WSLを起動し、後半のコマンドを打つが、簡略化して次で良い。
  cd C:\Users\c2997\source\repos\SELDLA\bin\Release\netcoreapp3.1\publish
  wsl
  bash ../../../../publish.sh
- WSLでzipコマンドを使うので、インストールしてなければ「sudo apt install zip」でインストールしておく。

```
#cd /mnt/c/Users/yoshi_000/Documents/Visual\ Studio\ Code/github/SELDLA/
#dotnet.exe publish -c Release -f netcoreapp2.0 -r linux-x64 -o SELDLA/linux-x64;
#dotnet.exe publish -c Release -f netcoreapp2.0 -r win-x64 -o SELDLA/win-x64;
#dotnet.exe publish -c Release -f netcoreapp2.0 -r osx-x64 -o SELDLA/osx-x64

ver=`cat ../../../../Program.cs|grep "string version ="|sed 's/.*= "//; s/".*//'`;
rm -rf ~/temp-SELDLA
mkdir -p ~/temp-SELDLA #NTFS上だとchmodが効かない
mv win mac linux ~/temp-SELDLA
cd ~/temp-SELDLA
chmod 644 linux/* mac/*
chmod 755 linux/SELDLA mac/SELDLA
mv win SELDLA-${ver}-win-x64
mv mac SELDLA-${ver}-mac-x64
mv linux SELDLA-${ver}-linux-x64
zip -r SELDLA-${ver}-win-x64.zip SELDLA-${ver}-win-x64 &
tar zvcf SELDLA-${ver}-linux-x64.tar.gz SELDLA-${ver}-linux-x64 &
tar zvcf SELDLA-${ver}-mac-x64.tar.gz SELDLA-${ver}-mac-x64 &
wait

scp SELDLA-${ver}-win-x64.zip SELDLA-${ver}-linux-x64.tar.gz SELDLA-${ver}-mac-x64.tar.gz sakura:www/software/SELDLA
```
- https://www.wix.com/account/sites?referralAdditionalInfo=Dashboard に行って編集。