using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using Mono.Options;

namespace SELDLA
{
    class Program
    {
        static void Main(string[] args)
        {
            string version = "2.1.1";
            int opt_dp = 1;
            int opt_gq = 0;
            double opt_nonzerorate = 0.3;
            double opt_p = 0.3;
            double opt_b = 0.1;
            double opt_s = 0.7;
            double opt_sm = 0.7;
            double opt_cm = 0.8;
            int opt_cs = 2;
            double opt_nc = 0.9;
            int opt_r = 10000;
            int opt_ldnum = 1;
            double opt_nonzerophase = 0.3;
            string opt_o = "seldla";
            string inputvcf = "";
            string inputfasta = "";
            string inputfamily = "";
            bool showHelp = false;
            bool nonewvcfout = false;
            bool maxLdClusterOnly = true;
            string precleaned = "";
            string mode = "crossbreed";
            double rateOfNotNASNP = 0.2;
            double rateOfNotNALD = 0.4;
            string removelqp = "no";
            bool needSort = false;
            int opt_ldseqnum = 2;
            //オプションとオプションの説明、そのオプションの引数に対するアクションを定義する
            var p = new OptionSet() {
                //{"c|cpu=", "number of cpus", (int v) => opt_c = v},
                {"DP=", "DP_threshold at the cleanupVcf step [1]", (int v) => opt_dp = v},
                {"GQ=", "GQ_threshold at the cleanupVcf step [0]", (int v) => opt_gq = v},
                {"NonZeroSampleRate=", "exclude ambiquous SNP at the cleanupVcf step (0-1) [0.3]", (double v) => opt_nonzerorate = v},
                {"p|hqsnp=", "high quality SNP rate at the splitVcf and Ld2Ph step [0.3]", (double v) => opt_p = v},
                {"b|bal=", "0 / 1 balance at the splitVcf step [0.1]", (double v) => opt_b = v},
                {"NeedSort","If the input vcf file is not sorted, use this option at the splitVcf step", v => needSort=v!=null},
                { "nl=", "near SNP match rate at the Snp2Ld step (0.5-1) [0.9]", (double v)=>opt_nc = v},
                {"r=", "the region to merge near SNP at the Snp2Ld step (bp) [10000]", (int v)=>opt_r =v},
                {"RateOfNotNASNP=", "threshold of the ratio that is not NA with each other when comparing SNP at the Snp2Ld step [0.2]", (double v) => rateOfNotNASNP = v},
                {"l|clmatch=", "cluster match rate at the Ld2Ph step [0.8]", (double v) => opt_cm = v},
                {"cs=", "cluster size at the Ld2Ph step [2]", (int v) => opt_cs = v},
                {"v|spmatch=", "split match rate at the Ld2Ph step (0.5-1) [0.7]", (double v) => opt_sm = v},
                {"ldnum=", "the minimum number of same LD at the Ld2Ph step [1]", (int v) => opt_ldnum=v},
                {"ldseqnum=", "the minimum number of consecutive LDs at the Ld2Ph step [1]", (int v) => opt_ldseqnum=v},
                {"UseAllLDClusters", "use all LD clusters at the Ld2Ph step", v => maxLdClusterOnly=v==null},
                {"RateOfNotNALD=", "threshold of the ratio that is not NA with each other when comparing LD at the LD2Ph step [0.4]", (double v) => rateOfNotNALD = v},
                {"RemoveLowQualityPhases=","remove low quality phases after the LD2Ph step (yes/no) [no]", v => removelqp = v},
                {"s|exmatch=", "extension match rate at the Chain step (0.5-1) [0.7]", (double v) => opt_s = v},
                {"NonZeroPhaseRate=", "exclude ambiquous Phase at the Chain step (0-1) [0.3]", (double v) => opt_nonzerophase = v},
                {"noNewVcf", "no converted vcf output with new position", v => nonewvcfout=v!=null},
                {"o|output=", "output prefix [seldla]", v => opt_o = v},
                {"vcf=", "input VCF file <required>", v => inputvcf = v},
                {"fasta=", "input FASTA file <required>", v => inputfasta = v},
                {"family=", "input family file <required>", v => inputfamily = v},
                {"precleaned=", "pre-calculated cleaned vcf file (if this option is used, input vcf is not used.)", v => precleaned = v},
                {"mode=", "analysis mode (crossbreed, haploid, duploid, selfpollination) [crossbreed]", v => mode=v},
                //VALUEをとらないオプションは以下のようにnullとの比較をしてTrue/Falseの値をとるようにする
                {"h|help", "show help.", v => showHelp = v != null}
            };
            if (args.Length == 0)
            {
                //args = @"-c 18 -p 0.3 -o e:\temp\output --vcf=e:\temp\snp.txt --fasta=c:\temp\test".Split(' ');
                //args = @"--fasta=E:\temp\testfasta.txt --vcf=E:\temp\100k.vcf --family=E:\temp\family.txt -o E:\temp\seldla".Split(' ');
                //args = @"--fasta=E:\temp\fourth_assembly.fasta --vcf=E:\temp\all.male.cln.vcf --family=E:\temp\head.txt -o E:\temp\test2".Split(' ');
                //args = @"--fasta=E:\temp\fourth_assembly.fasta --vcf=E:\temp\all.male.cln.vcf --family=E:\temp\head2.txt -o E:\temp\test6".Split(' ');
                //args = @"--fasta=E:\temp\BG_Platanus_scaffold_ver20150818.fa --vcf=E:\temp\BG_GRASDi_flaged_raw_variants_PASS_refusion_parents.vcf --family=E:\temp\family.txt2 -o E:\temp\testblue3 --DP=5 --GQ=20 -r 10000 --mode=duploid".Split(' ');
                //args = @"--precleaned=E:\temp\testblue4_clean.txt --fasta=E:\temp\BG_Platanus_scaffold_ver20150818.fa --vcf=E:\temp\BG_GRASDi_flaged_raw_variants_PASS_refusion_parents.vcf --family=E:\temp\family.txt2 -o E:\temp\testblue4 --DP=5 --GQ=20 -r 100 --cs=3 --mode=duploid".Split(' ');
                //args = @"--fasta=E:\temp\BG_Platanus_scaffold_ver20150818.fa --vcf=E:\temp\BG_GRASDi_flaged_raw_variants_PASS_refusion_parents.vcf --family=E:\temp\family.txt2 -o E:\temp\testblue9 --DP=5 --GQ=20 -r 100 --cs=3 --mode=duploid --MaxLdClusterOnly".Split(' ');
                //args = @"--fasta=E:\temp\BG_Platanus_scaffold_ver20150818.fa --vcf=E:\temp\BG_GRASDi_flaged_raw_variants_PASS_refusion_parents.vcf --family=E:\temp\family.txt2 -o E:\temp\testblue10 --DP=5 --GQ=20 -r 100 --cs=3 --mode=duploid --MaxLdClusterOnly --noNewVcf --RateOfNotNASNP=0.2 --RateOfNotNALD=0.4".Split(' ');
                //args = @"--fasta=E:\temp\testblue7_include_unordered_in_chr.fasta --vcf=E:\temp\testblue7_newpos_include_unordered_in_chr.vcf --family=E:\temp\family.txt2 -o E:\temp\testblue8 --DP=5 --GQ=20 -r 1000 --cs=3 --mode=duploid --noNewVcf".Split(' ');
                //args = @"--precleaned=E:\temp\testpfu_clean.txt --fasta=E:\temp\pfu_genome2.0.fasta --vcf=E:\temp\pfu.all.clean.vcf --family=E:\temp\pfu.fam.txt -o E:\temp\testpfu --mode duploid --DP=10 --GQ=20 -r 1000 -cs 1 -s 0.6 -p 0.5 -v 0.7 --MaxLdClusterOnly".Split(' ');
                //args = @"--fasta=E:\temp\pnyererei1.fasta --vcf=E:\temp\pnyererei.vcf --family=E:\temp\pnyererei.fam -o E:\temp\testpnye --GQ=20 -r 100 --cs=3 --mode=duploid".Split(' ');
                //args = @"--fasta=E:\temp\testpnye_ext1.fasta --vcf=E:\temp\testpnye_newpos.vcf --family=E:\temp\pnyererei.fam -o E:\temp\test2pnye --GQ=20 -r 1000 --cs=3 --mode=duploid --noNewVcf".Split(' ');
                //args = @"--fasta=E:\temp\testpnye_include_unordered_in_chr.fasta --vcf=E:\temp\testpnye_newpos_include_unordered_in_chr.vcf --family=E:\temp\pnyererei.fam -o E:\temp\test3pnye --GQ=20 -r 1000 --cs=3 --mode=duploid --noNewVcf".Split(' ');
                //args = @"--fasta=E:\temp\testpnye_ext1.fasta --vcf=E:\temp\testpnye_newpos.vcf --family=E:\temp\pnyererei.fam2 -o E:\temp\test2pnye --GQ=20 -r 100 --DP=5 --cs=3 --mode=duploid --noNewVcf --precleaned=E:\temp\test2pnye_clean.txt --MaxLdClusterOnly".Split(' ');
                //args = @"--fasta=E:\temp\simfugu.fa --vcf=E:\temp\simfugu.vcf --family=E:\temp\simfugufam.txt -o E:\temp\simfugutemp --DP=5 -r 100 --cs=3 --mode=duploid --noNewVcf".Split(' ');
                //args = @"--fasta=E:\temp\dpulex_v1.1_scaf.fasta --vcf=E:\temp\dpulex.all.vcf --family=E:\temp\dpulex.family.txt -o E:\temp\dpulex1 --cs=3 --mode=haploid --MaxLdClusterOnly --noNewVcf --precleaned=E:\temp\dpulex1_clean.txt --nl=0.8 -l 0.7 --RateOfNotNASNP=0.3 --RateOfNotNALD=0.9 --clmatch=0.9 -r 10000".Split(' ');
                //args = @"--fasta=E:\temp\seldla-selfpoll\RSA_r2.0.fasta --vcf=E:\temp\seldla-selfpoll\ASF2-sakurajima.recode.vcf --family=E:\temp\seldla-selfpoll\ASF2-sakurajima.recode.family.txt -o E:\temp\seldla-selfpoll\selfpoll --cs=2 --mode=selfpollination --MaxLdClusterOnly --noNewVcf -r 1000".Split(' ');
                //args = @"--fasta=E:\temp\suma\suma_draft_genome.fasta --vcf=E:\temp\suma\suma_second.vcf --family=E:\temp\suma\family_suma.txt -o E:\temp\suma\suma --noNewVcf".Split(' ');
                args = @"--fasta=C:\work\sample_itoyo.fa --vcf=C:\work\sample_itoyo_1-100_head1m.txt --precleaned=C:\work\sample_itoyo_1-100_head1m.txt --family=C:\work\sample_itoyo_family.txt -o C:\work\out_itoyo2 --mode=haploid --noNewVcf -p 0.03 -b 0.03 --cs 2 --nl 0.9 --NonZeroSampleRate=0.05 --NonZeroPhaseRate=0.1 -r 4000 --RateOfNotNASNP=0.001 --RateOfNotNALD=0.01".Split(' ');
                //args = @"--fasta=C:\work\sample_itoyo.fa --vcf=C:\work\pseudochr.re.fa.removedup.matrix.clean.txt.vcf2.single.1-110 --precleaned=C:\work\pseudochr.re.fa.removedup.matrix.clean.txt.vcf2.single.1-110 --family=C:\work\sample_itoyo_family.txt -o C:\work\out_itoyo3 --mode=haploid --noNewVcf -p 0.03 -b 0.03 --cs 2 --nl 0.9 --NonZeroSampleRate=0.05 --NonZeroPhaseRate=0.1 -r 4000 --RateOfNotNASNP=0.001 --RateOfNotNALD=0.01".Split(' ');
                //dotnet publish -c Release -f netcoreapp2.0 -r linux-x64 -o SELDLA/linux-x64
                //dotnet publish -c Release -f netcoreapp2.0 -r win-x64 -o SELDLA/win-x64
                //dotnet publish -c Release -f netcoreapp2.0 -r osx-x64 -o SELDLA/osx-x64
                //dotnet publish -c Release -f netcoreapp2.0 -r linux-x64 -o SELDLA/linux-x64;dotnet publish -c Release -f netcoreapp2.0 -r win-x64 -o SELDLA/win-x64;dotnet publish -c Release -f netcoreapp2.0 -r osx-x64 -o SELDLA/osx-x64
            }
            try
            {
                var extra = p.Parse(args);
                extra.ForEach(t => Console.WriteLine("invalid parameter: " + t));
                if (extra.Count > 0)
                {
                    return;
                }
            }
            //パースに失敗した場合OptionExceptionを発生させる
            catch (OptionException e)
            {
                Console.WriteLine("Option parse error:");
                Console.WriteLine(e.Message);
                Console.WriteLine("Try `--help' for more information.");
                return;
            }

            if (inputfasta == "")
            {
                Console.WriteLine("no input fasta");
                showHelp = true;
            }
            if (inputvcf == "")
            {
                Console.WriteLine("no input vcf");
                showHelp = true;
            }
            if (inputfamily == "")
            {
                Console.WriteLine("no input family");
                showHelp = true;
            }
            if (!(mode == "crossbreed" || mode == "haploid" || mode == "duploid" || mode == "selfpollination"))
            {
                Console.WriteLine("unrecognized mode");
            }

            if (showHelp)
            {
                Console.WriteLine("SELDLA ver" + version);
                p.WriteOptionDescriptions(Console.Out);
                return;
            }

            Console.WriteLine("Run SELDLA ver " + version);

            //入力ファイルをフィルタリング、家系ごと(family.txtの行ごと)に分割
            Prepare prep = new Prepare();
            if (precleaned == "")
            {
                Console.WriteLine("Clean up VCF");
                prep.cleanupVcf(inputvcf, opt_dp, opt_gq, opt_nonzerorate, opt_o);
                Console.WriteLine("Split VCF into each family");
                prep.splitVcf(opt_o + "_clean.txt", opt_o, inputfamily, opt_p, opt_b, mode, needSort);
            }
            else
            {
                Console.WriteLine("Split VCF into each family");
                prep.splitVcf(precleaned, opt_o, inputfamily, opt_p, opt_b, mode, needSort);
            }

            StreamReader file = new StreamReader(inputfamily);
            string line;
            List<string[]> families = new List<string[]>();
            int num_member = 0;
            while ((line = file.ReadLine()) != null)
            {
                string[] temp = line.Split("\t");
                if (temp.Length >= 3)
                {
                    families.Add(temp);
                    if (mode == "crossbreed" || mode == "duploid")
                    {
                        num_member += temp.Length - 2;
                    }
                    else
                    {
                        num_member += temp.Length - 1;
                    }
                }
            }
            file.Close();

            //1回目のフェージング　分割点を探すために行う
            int num_fam = 0;
            Dictionary<string, SortedDictionary<int, int>> breaks = new Dictionary<string, SortedDictionary<int, int>>();
            foreach (string[] fam in families)
            {
                num_fam++;
                Snp2Ld snp = new Snp2Ld();
                Console.WriteLine("SNP to Block in family No. " + num_fam);
                snp.run(opt_o + "_split_" + num_fam + ".txt", opt_nc, opt_r, rateOfNotNASNP);

                Ld2Ph ld = new Ld2Ph();
                Console.WriteLine("Block to Phase in family No. " + num_fam);
                ld.run2(opt_o + "_split_" + num_fam + ".txt.ld", opt_b, opt_cm, opt_cs, opt_sm, true, opt_ldnum, maxLdClusterOnly, rateOfNotNALD, opt_ldseqnum);

                int counter = 0;
                file = new System.IO.StreamReader(opt_o + "_split_" + num_fam + ".txt.ld.break");
                while ((line = file.ReadLine()) != null)
                {
                    counter++;
                    //System.Console.WriteLine(line);
                    if (counter > 1)
                    {
                        string[] temp = line.Split("\t");
                        if (!breaks.ContainsKey(temp[0]))
                        {
                            SortedDictionary<int, int> newbreaks = new SortedDictionary<int, int>();
                            breaks.Add(temp[0], newbreaks);
                        }
                        if (!breaks[temp[0]].ContainsKey(Int32.Parse(temp[1])))
                        {
                            breaks[temp[0]].Add(Int32.Parse(temp[1]), 1);
                        }
                        else
                        {
                            breaks[temp[0]][Int32.Parse(temp[1])]++;
                        }
                        if (!breaks[temp[0]].ContainsKey(Int32.Parse(temp[2])))
                        {
                            breaks[temp[0]].Add(Int32.Parse(temp[2]), -1);
                        }
                        else
                        {
                            breaks[temp[0]][Int32.Parse(temp[2])]--;
                        }
                    }
                }
            }

            //FASTAを分割
            Dictionary<string, string> refseqs = new Dictionary<string, string>();
            StreamReader filefasta = new StreamReader(inputfasta);
            string inchr = "";
            StringBuilder insb = new StringBuilder();
            while ((line = filefasta.ReadLine()) != null)
            {
                if (line.StartsWith(">"))
                {
                    if (inchr != "")
                    {
                        refseqs.Add(inchr, insb.ToString());
                    }
                    inchr = line.Substring(1).Split(" ")[0].Split("\t")[0];
                    insb = new StringBuilder();
                }
                else
                {
                    insb.Append(line);
                }
            }
            if (inchr != "")
            {
                refseqs.Add(inchr, insb.ToString());
            }
            filefasta.Close();

            Console.WriteLine("Detect breakpoints");
            StreamWriter writer = new StreamWriter(opt_o + "_break.txt");
            Dictionary<string, List<int>> breaklist = new Dictionary<string, List<int>>();
            foreach (KeyValuePair<string, SortedDictionary<int, int>> tempbreak in breaks)
            {
                //Console.WriteLine(tempbreak.Key);
                bool incl = false;
                int num_split = 0;
                int oldkey = 0;
                foreach (KeyValuePair<int, int> pair in tempbreak.Value)
                {
                    if (incl && pair.Value < 0)
                    {
                        num_split++;
                        //Console.WriteLine(tempbreak.Key+", "+(oldkey)+", "+pair.Key);
                        int breakpos = oldkey + searchN(refseqs[tempbreak.Key].Substring(oldkey - 1, pair.Key - oldkey + 1));
                        //Console.WriteLine(tempbreak.Key+"\t"+oldkey+"\t"+pair.Key+"\t"+breakpos);
                        writer.WriteLine(tempbreak.Key + "\t" + oldkey + "\t" + pair.Key + "\t" + breakpos);
                        if (!breaklist.ContainsKey(tempbreak.Key))
                        {
                            List<int> templist = new List<int>();
                            breaklist.Add(tempbreak.Key, templist);
                        }
                        breaklist[tempbreak.Key].Add(breakpos);
                    }
                    if (pair.Value >= 0) { incl = true; } else { incl = false; }
                    oldkey = pair.Key;
                }
            }
            writer.Close();

            Dictionary<string, string> refseqs2 = new Dictionary<string, string>();
            writer = new StreamWriter(opt_o + "_split_seq.txt");
            foreach (KeyValuePair<string, string> chr in refseqs)
            {
                if (!breaklist.ContainsKey(chr.Key))
                {
                    writer.WriteLine(chr.Key + "\t" + chr.Value);
                    refseqs2.Add(chr.Key, chr.Value);
                }
                else
                {
                    int num_break = 0;
                    int old_end = 0;
                    foreach (int pos in breaklist[chr.Key])
                    {
                        num_break++;
                        writer.WriteLine(chr.Key + "_" + num_break + "\t" + chr.Value.Substring(old_end, pos - old_end));
                        refseqs2.Add(chr.Key + "_" + num_break, chr.Value.Substring(old_end, pos - old_end));
                        old_end = pos;
                    }
                    num_break++;
                    writer.WriteLine(chr.Key + "_" + num_break + "\t" + chr.Value.Substring(old_end, chr.Value.Length - old_end));
                    refseqs2.Add(chr.Key + "_" + num_break, chr.Value.Substring(old_end, chr.Value.Length - old_end));
                }
            }
            writer.Close();

            //2回目のフェージング
            Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>> datas
             = new Dictionary<string, SortedDictionary<int, Dictionary<string, int[]>>>();
            for (int i = 1; i <= families.Count; i++)
            {
                StreamReader ldfile = new StreamReader(opt_o + "_split_" + i + ".txt.ld");
                StreamWriter ldbfile = new StreamWriter(opt_o + "_split_" + i + ".txt.ld2");
                int counter = 0;
                while ((line = ldfile.ReadLine()) != null)
                {
                    counter++;
                    //System.Console.WriteLine(line);
                    if (counter == 1)
                    {
                        ldbfile.WriteLine(line);
                    }
                    else
                    {
                        string[] vals = line.Split("\t");
                        if (!breaklist.ContainsKey(vals[1]))
                        {
                            ldbfile.WriteLine(line);
                        }
                        else
                        {
                            int num_break = 1;
                            int old_end = 0;
                            foreach (int pos in breaklist[vals[1]])
                            {
                                if (old_end < Int32.Parse(vals[2]) && Int32.Parse(vals[2]) <= pos) { break; }
                                old_end = pos;
                                num_break++;
                            }
                            int temppos = Int32.Parse(vals[2]) - old_end;
                            ldbfile.Write(vals[0] + "\t" + vals[1] + "_" + num_break + "\t" + temppos);
                            for (int j = 3; j < vals.Length; j++)
                            {
                                ldbfile.Write("\t" + vals[j]);
                            }
                            ldbfile.WriteLine("");
                        }
                    }
                }
                ldfile.Close();
                ldbfile.Close();

                Ld2Ph ld = new Ld2Ph();
                Console.WriteLine("corrected Block to Phase in family No. " + i);
                ld.run2(opt_o + "_split_" + i + ".txt.ld2", opt_b, opt_cm, opt_cs, opt_sm, false, opt_ldnum, maxLdClusterOnly, rateOfNotNALD, opt_ldseqnum);


                StreamReader phfile = new StreamReader(opt_o + "_split_" + i + ".txt.ld2.ph");
                int numNR = 0;
                while ((line = phfile.ReadLine()) != null)
                {
                    if (numNR == 0) { line = phfile.ReadLine(); } //ヘッダーを飛ばす
                    numNR++;
                    string[] vals = line.Split("\t");
                    if (vals[1] != "lowqual" || removelqp != "yes")
                    {
                        if (!datas.ContainsKey(vals[0]))
                        {
                            SortedDictionary<int, Dictionary<string, int[]>> chrdatas = new SortedDictionary<int, Dictionary<string, int[]>>();
                            datas.Add(vals[0], chrdatas);
                        }
                        if (!datas[vals[0]].ContainsKey(i))
                        {
                            Dictionary<string, int[]> posdatas = new Dictionary<string, int[]>();
                            datas[vals[0]].Add(i, posdatas);
                        }
                        int[] phdatas = new int[vals.Length - 3];
                        for (int j = 3; j < vals.Length; j++)
                        {
                            phdatas[j - 3] = Int32.Parse(vals[j]);
                        }
                        if (numNR % 2 == 1)
                        {
                            datas[vals[0]][i].Add("start", phdatas);
                        }
                        else
                        {
                            datas[vals[0]][i].Add("end", phdatas);
                        }
                    }
                }
                phfile.Close();
            }

            //連鎖するコンティグを伸ばしていく
            Chain cs = new Chain();
            Console.WriteLine("make linkage map...");
            cs.run(refseqs2, datas, opt_s, opt_nonzerophase, num_member, opt_o);

            //新しい座標にVCFを変換する
            if (!nonewvcfout)
            {
                ConvVcf cv = new ConvVcf();
                Console.WriteLine("convert VCF to new position");
                cv.run(inputvcf, opt_o + "_break.txt", opt_o + "_chain.txt", opt_o, refseqs2);
                ConvSNP newsnp = new ConvSNP();
                Console.WriteLine("convert SNP to new position");
                newsnp.run(num_fam, opt_o + "_break.txt", opt_o + "_chain.txt", opt_o, refseqs2);
            }
        }

        public static int searchN(string seq)
        {
            int res = 0;
            int bstart = 0;
            int bend = seq.Length;
            for (int i = 0; i < seq.Length; i++)
            {
                if (seq.Substring(i, 1).ToLower() == "n")
                {
                    bstart = i;
                    break;
                }
            }
            if (bstart == 0)
            {
                return seq.Length / 2;
            }
            for (int i = seq.Length - 1; i >= 0; i--)
            {
                if (seq.Substring(i, 1).ToLower() == "n")
                {
                    bend = i;
                    break;
                }
            }
            if (bstart < seq.Length - 1 - bend)
            {
                res = bstart;
            }
            else
            {
                res = bend;
            }
            return res;
        }


    }
}
