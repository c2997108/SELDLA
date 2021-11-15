# About

SELDLA (Scaffold Extender with Low Depth Linkage Analysis) is the tool for generating linkage maps and draft genomes with low depth (<1x) sequencing data. SELDLA is written by C# with .NET Core, and it is executable on Windows (7-), Mac (10.12-) and Linux (Ubuntu 14.04-, CentOS 6-).

![image](https://user-images.githubusercontent.com/5350508/141719091-951be7eb-d027-433e-b75b-7e6a204d89e2.png)

# How to install

## Windows

Download the zip file from [here](https://c29979108.wixsite.com/seldla) and unzip.
Then, open PowerShell or Command prompt, and run SELDLA.exe.

## Mac or Linux

Install .NET Core Runtime.

https://www.microsoft.com/net/download

In summary (in CentOS7),

```
sudo rpm -Uvh https://packages.microsoft.com/config/rhel/7/packages-microsoft-prod.rpm
sudo yum update
sudo yum install aspnetcore-runtime-2.2
sudo yum install libgdiplus-devel
```

Then, download SELDLA binary from [here](https://c29979108.wixsite.com/seldla) and unzip.

```
unzip SELDLA_x.x.x.zip
```

Then, run SELDLA binary.

```
SELDLA_v2.x.x/linux-x64/SELDLA
```

# How to use

## 1. SELDLA needs 3 files.

 -  FASTA file (to be extended)

 -  VCF file (SNV information on the above FASTA file)

 -  family file  (tell father and mother ID, tab separated)


   In crossbreed mode, the first column indicates the half parent's ID of which genome you want to extend, and the second column is another parent's ID in the VCF file. The columns after the third column indicate children's IDs. 

   In haploid mode, first column indicates a parent's ID in the VCF file, and the columns after the second column indicate children's IDs.

   In duploid mode (such as RAD-seq), 1 family needs to be written in 2 rows. The first and the second columns of the first row indicate father's and mother's IDs, and the first and the second columns of the second row indicates mother's and father's IDs.

   Example files can be downloaded [here](http://suikou.fs.a.u-tokyo.ac.jp/yosh_data/SELDLA/).

## 2. Run SELDLA

 Recommended options to try first   

  - For crossbreed mode

```
   SELDLA_v2.x.x/linux-x64/SELDLA --vcf=male.vcf --fasta=fourth_assembly.fasta --family=family.txt --mode=crossbreed
```

  - For duploid mode assuming RAD-seq etc.

```
   SELDLA_v2.x.x/linux-x64/SELDLA --vcf=input.vcf --fasta=assembly.fasta --family=family.txt --mode=duploid --DP=5 --GQ=20 -r 100 --cs=3 --MaxLdClusterOnly --noNewVcf
```

   You will get the extended FASTA file, the lift overed vcf file and the linkage map like bellow.

![image](https://user-images.githubusercontent.com/5350508/141718996-d0a97943-8806-4b86-b83c-f4a80ae62e3e.png)

## 3. Tuning SELDLA

   The parameter that has the biggest impact is --exmatch, which when lowered to 0.5, all contigs are connected to one. The default value is 0.7, but you can try lowering it to 0.65 or so. The threshold for splitting a misassembled contig is the option of --spmatch, this disconnects the contig in phases with a match rate below the threshold. If this value is set to 0.5, no contigs will be disconnected. The recommended value is between 0.5 and the value of --exmatch.

- Overview of SELDLA principles

![image](https://user-images.githubusercontent.com/5350508/141719038-22fcf4c2-cde3-4e1d-8307-c0c523207e6f.png)
![image](https://user-images.githubusercontent.com/5350508/141719054-6f9a72ed-5056-41f6-b572-49fa16737f3f.png)
![image](https://user-images.githubusercontent.com/5350508/141719065-d59eb8de-e993-4e8a-bbcd-93e6895910c9.png)

