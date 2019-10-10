# リリース手順
# WSLでbash note.sh 2.0.8などと使う

ver=$1 #2.0.8

cd /mnt/c/Users/yoshi_000/Documents/Visual\ Studio\ Code/github/SELDLA/
dotnet.exe publish -c Release -f netcoreapp2.0 -r linux-x64 -o SELDLA/linux-x64;
dotnet.exe publish -c Release -f netcoreapp2.0 -r win-x64 -o SELDLA/win-x64;
dotnet.exe publish -c Release -f netcoreapp2.0 -r osx-x64 -o SELDLA/osx-x64

mkdir -p ~/test/SELDLA
mv SELDLA/win-x64 ~/test/SELDLA/SELDLA-${ver}-win-x64
mv SELDLA/osx-x64 ~/test/SELDLA/SELDLA-${ver}-mac-x64
mv SELDLA/linux-x64 ~/test/SELDLA/SELDLA-${ver}-linux-x64
cd ~/test/SELDLA
zip -r SELDLA-${ver}-win-x64.zip SELDLA-${ver}-win-x64 &
cd ~/test/SELDLA/SELDLA-${ver}-linux-x64
chmod 644 *
cd ~/test/SELDLA/SELDLA-${ver}-mac-x64
chmod 644 *
cd ~/test/SELDLA
chmod 755 SELDLA-${ver}-linux-x64/SELDLA SELDLA-${ver}-mac-x64/SELDLA
tar zvcf SELDLA-${ver}-linux-x64.tar.gz SELDLA-${ver}-linux-x64 &
tar zvcf SELDLA-${ver}-mac-x64.tar.gz SELDLA-${ver}-mac-x64 &
wait

scp ~/test/SELDLA/*.zip ~/test/SELDLA/*.tar.gz sakura:www/software/SELDLA

# https://www.wix.com/account/sites?referralAdditionalInfo=Dashboard に行って編集。