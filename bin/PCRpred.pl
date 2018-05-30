#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use Statistics::Descriptive;
use PDL;

#=============================the path of software=============================

my $software_NWalign_path = "../software/sequence_alignment/NWalign";
#sequence
my $software_domain_path = "../software/domain/PfamScan/pfam_scan.pl";
my $software_disorder_path = "../software/disorder/DisEMBL-1.4/DisEMBL.py";
my $software_sift_path = "/data4/hfliu/software/sift6.2.1";
my $software_sift_input_path = "../software/sift/input";
my $software_polyphen_path = "/data4/hfliu/software/polyphen-2.2.2";
my $software_polyphen_input_path = "../software/polyphen-2/input";
#structure
my $software_pocket_path = "../software/fpocket/fpocket2/bin/fpocket";
my $software_dssp_path = "../software/dssp/dssp-2.0.4-linux-amd64";
my $software_psaia_path = "../software/PSAIA/PSAIA-1.0";
my $software_hbplus_path = "../software/hbplus/hbplus";
#end

#=============================the path of databases=============================

my $nr_database_path = "~/nr";
my $uniref90_database_path = "~/uniref90";
my $domain_pfam_path = "../software/domain/";
#end

#=============================the path of result=============================

my $result_file = "../result_files";

my $PDB_chain = "$result_file/PDB_chain";
my $ca = "$result_file/PDB_chain/CA";
my $seq2str_align = "$result_file/sequence_alignment";
#sequence
my $co_evolution = "$result_file/co-evolution";
my $domain = "$result_file/domain";
my $disorder = "$result_file/disorder";
my $shannon_entropy = "$result_file/shannon_entropy";
my $sift = "$result_file/sift";
my $polyphen = "$result_file/polyphen";
#structure
my $secondary_structure = "$result_file/secondary_structure";
my $psaia = "$result_file/PSAIA";
my $dp_index = "$result_file/dp_index";
my $hbplus = "$result_file/hbplus";
my $laplacian_norm = "$result_file/laplacian_norm";
my $bfactor = "$result_file/bfactor";
my $residue_network = "$result_file/TP/residue_network";
my $topological_features = "$result_file/TP/topological_features";
#end

#=============================parse arguments=============================

my ($pre_file, $seq_path, $str_path, $out_path);

if(scalar@ARGV != 8) {
    &print_help();
    exit();
}else{
    foreach my $i(0..$#ARGV) {
        if ($ARGV[$i] eq '-f') {
            $pre_file = $ARGV[$i+1];
        }elsif($ARGV[$i] eq '-seq') {
            $seq_path = $ARGV[$i+1];
        }elsif($ARGV[$i] eq '-str') {
            $str_path = $ARGV[$i+1];
        }elsif($ARGV[$i] eq '-out') {
            $out_path = $ARGV[$i+1];
        }
    }
}
#end

#=============================unique information=============================

open(IN,"<","$pre_file") or die "$!";
chomp(my @pre_file = <IN>);
close IN;

my %unique_seq;
my %pdb_seq;
my %pdb_site;
my %pdb_avail_site;   #both sites are in the PDB

foreach my $i(@pre_file) {
    my @tem = split(/\s+/,$i);
    if (exists $unique_seq{$tem[0]}) {
        push @{$unique_seq{$tem[0]}},$i;
    }else{
        my @bb;
        push @bb,$i;
        $unique_seq{$tem[0]} = \@bb;
    }
    if ($tem[5] ne 'NA') {
        $pdb_seq{"$tem[5]\t$tem[6]"} = $tem[0];
        if (exists $pdb_site{"$tem[5]\t$tem[6]"}) {
            push @{$pdb_site{"$tem[5]\t$tem[6]"}},$i;
        }else{
            my @bb;
            push @bb,$i;
            $pdb_site{"$tem[5]\t$tem[6]"} = \@bb;
        }
        
    }
}
#end

#=============================sequence alignment=============================

my %aa=('ALA'=>'A','ARG'=>'R','ASN'=>'N','ASP'=>'D','CYS'=>'C','GLN'=>'Q','GLU'=>'E','GLY'=>'G','HIS'=>'H','ILE'=>'I',
        'LEU'=>'L','LYS'=>'K','MET'=>'M','PHE'=>'F','PRO'=>'P','SER'=>'S','THR'=>'T','TRP'=>'W','TYR'=>'Y','VAL'=>'V');
my %aa_single = reverse%aa;

if (scalar(keys%pdb_seq) ne 0) {
    foreach my $key(keys%pdb_seq) {
        my ($pdb_name,$pdb_chain) = split(/\s+/,$key);
        #extract the chain from the PDB file and generate seqpdb file
        open(OUT,">","$PDB_chain/$pdb_name-$pdb_chain.pdb") or die "$!";
        open(IN,"<","$str_path/$pdb_name.pdb") or die "$!";
        chomp(my @pdb = <IN>);
        close IN;
        my @realpdb;
        my $n = 0;
        my $site;
        my %atom_num;
        my %num_resaa;
        #check for multiple structural models
        foreach my $p(0..$#pdb) {
            if ($pdb[$p] =~ /^ENDMDL/) {
                $n++;
                $site = $p;
                last;
            }
        }
        if ($n > 0) {
            foreach my $j(0..$site){
                push @realpdb,$pdb[$j];
            }
        }else{
            @realpdb = @pdb;
        }
        foreach my $p(@realpdb) {
            if ($p =~ /^(ATOM|TER)/) {
                my $chain = substr($p,21,1);
                my $resaa = substr($p,17,3);
                my $num = substr($p,22,5);
                $num =~ s/\s+//g;
                my $atom = substr($p,5,6);
                $atom =~ s/\s+//g;
                if (exists $aa{$resaa} and $chain eq $pdb_chain) {
                    print OUT"$p\n";
                    $atom_num{$num} = $atom;
                    $num_resaa{$num} = $aa{$resaa};
                } 
            }
        }
        close OUT;
        #generate seqpdb file
        my %seqpdb2pdb;
        my $seqnumber = 0;
        open(OUTS,">","$seq2str_align/$pdb_name-$pdb_chain.fasta") or die "$!";
        print OUTS">$pdb_name-$pdb_chain\n";
        foreach my $kk(sort{$atom_num{$a} <=> $atom_num{$b}}keys%atom_num) {
            print OUTS "$num_resaa{$kk}";
            $seqnumber++;
            $seqpdb2pdb{$seqnumber} = $kk;
        }
        close OUTS;
        #sequence alignment by NWalign
        system "$software_NWalign_path $seq_path/$pdb_seq{$key}.fasta $seq2str_align/$pdb_name-$pdb_chain.fasta > $seq2str_align/$pdb_seq{$key}-$pdb_name-$pdb_chain.txt";
        #mapping sequence site to PDB site
        my %seq2str;
        open(IN,"<","$seq2str_align/$pdb_seq{$key}-$pdb_name-$pdb_chain.txt") or die "$!";
        chomp(my @align = <IN>);
        close IN;
        my $colnumber = 0;
        foreach my $i(6..$#align) {
            if ($align[$i] =~ /1234567890/) {
                $colnumber++;
            }
        }
        my @seq1_seq;
        my @seq2_seqpdb;
        my @map;
        foreach my $i(7..(6 + $colnumber)) {
            push @seq1_seq,$align[$i];
        }
        foreach my $i((7 + $colnumber)..(6 + 2 * $colnumber)) {
            push @map,$align[$i];
        }
        foreach my $i((7 + 2 * $colnumber)..(6 + 3 * $colnumber)) {
            push @seq2_seqpdb,$align[$i];
        }
        my $str1 = join('',@seq1_seq);
        my @s1_seq = split(//,$str1);
        my $str2 = join('',@map);
        my @mapping = split(//,$str2);
        my $str3 = join('',@seq2_seqpdb);
        my @s2_seqpdb = split(//,$str3);
        my $n1 = 0;
        my $n2 = 0;
        foreach my $i(0..$#s1_seq) {
            if ($s1_seq[$i] ne '-') {
                $n1++;
            }
            if ($s2_seqpdb[$i] ne '-') {
                $n2++;
            }
            if ($mapping[$i] eq ':') {
                $seq2str{$n1} = $seqpdb2pdb{$n2};
            }
        }
        #check whether both sites are in the PDB
        foreach my $kk(@{$pdb_site{$key}}) {
            my @tem = split(/\s+/,$kk);
            my $site1 = $tem[2];
            my $site2 = $tem[4];
            if (exists $seq2str{$site1} and exists $seq2str{$site2}) {
                my $str4 = "$kk\t$seq2str{$site1}\t$seq2str{$site2}";
                if (exists $pdb_avail_site{$key}) {
                    push @{$pdb_avail_site{$key}},$str4;
                }else{
                    my @bb1;
                    push @bb1,$str4;
                    $pdb_avail_site{$key} = \@bb1;
                }
            }else{
                print "$kk\nThe PDB structure is not available because the two sites do not exist in the PDB at the same time.\n\n";
            }  
        }
    }
}
#end

#=============================available structural information=============================

if (scalar(keys%pdb_avail_site ne 0)) {
    print "structural information\n";
    print "protein\tresidue1\tsite1\tresidue2\tsite2\tPDB_ID\tchain\tPDB_site1\tPDB_site2\n";
    foreach my $key(keys%pdb_avail_site) {
        foreach my $key1(@{$pdb_avail_site{$key}}) {
            print "$key1\n";
        }
    }
}
print "sequence information\n";
print "protein\tresidue1\tsite1\tresidue2\tsite2\tPDB_ID\tchain\n";
foreach my $key(keys%unique_seq) {
    foreach my $key1(@{$unique_seq{$key}}) {
        print "$key1\n";
    }
}
#end

#=============================calculate sequence features=============================
=cut
foreach my $key(keys%unique_seq) {
    #run psi-blast to generate MSA
    system "psiblast -evalue 0.001 -line_length 10000 -num_alignments 5000 -num_iterations 1 -num_threads 4 -query $seq_path/$key.fasta -db $nr_database_path/nr -outfmt 4 -out $co_evolution/$key.psiblast";
    &generate_MSA("$co_evolution/$key.psiblast");
    #calculate co-evolution
    open(OU,">","./new_protein.txt") or die "$!";
    print OU"$key\n";
    close OU;
    system "matlab matlab -nojvm -nodisplay -nosplash -nodesktop < DI.m >running.log 2>running.err";
    &zscore("$co_evolution/$key.txt",5);
    #calculates the disorder
    system "python $software_disorder_path 8 8 4 1.2 1.4 1.2 $seq_path/$key.fasta > $disorder/$key.txt";
    #calculate domain
    system"perl $software_domain_path -fasta $seq_path/$key.fasta -dir $domain_pfam_path -outfile $domain/$key.txt";
    #run psi-blast to generate PSSM
    system "psiblast -evalue 0.001 -num_iterations 3 -num_threads 4 -query $seq_path/$key.fasta -db $nr_database_path/nr -out_ascii_pssm $shannon_entropy/$key.pssm";
    #calculate shannon entropy
    my @blosums = ();
    my %hash_blosum = ();
    my $line_n = 0;
    open(FA, "<","./BLOSUM62") or die "$!";
    while (<FA>) {
        chomp;
        s/^\s+//g;
        my @item = split/\s+/;
        last if $item[0] = ~/^B/;
        if (scalar @item >= 26) {
            my @ref = @item[1..20];
            push(@blosums,[@ref]);
            for(my $i = 0;$i<20;$i++) {
                $hash_blosum{$line_n}{$i} = $ref[$i];
            }
            $line_n++;
        }
    }
    close FA;
    my $blosum = pdl(@blosums[0..19]);
    my $pssm_file = "$shannon_entropy/$key.pssm";
    my ($seq,$aa,$se) = &JSD($pssm_file);
    my @seq_num = @{$seq};
    my @aa_key = @{$aa};
    my @se = @{$se};
    my $sen = &zscore_logistic(@se);
    my @sen = @{$sen};
    open(OUT,">","$shannon_entropy/$key.txt") or die "$!";
    for(my $i = 0;$i < @sen; $i++) {
        print OUT $seq_num[$i],"\t",$aa_key[$i],"\t",$sen[$i],"\n";
    }
    close OUT;
    #calculate predicted probability of pathogenic/ SIFT prediction and PolyPhen-2 prediction
    #generate input file
    open(OUT,">","$software_sift_input_path/$key.subst") or die "$!";
    open(OUTT,">","$software_polyphen_input_path/$key.txt") or die "$!";
    my %sift_muta;
    my %polyphen_muta;
    foreach my $key1(@{$unique_seq{$key}}) {
        my @tem = split(/\s+/,$key1);
        foreach my $kk(sort{$a cmp $b}keys%aa_single) {
            if ($tem[1] ne $kk) {
                $sift_muta{"$tem[1]$tem[2]$kk\n"} = 0;
                $polyphen_muta{"$key\t$tem[2]\t$tem[1]\t$kk\n"} = 0;
            }
        }
        foreach my $kk(sort{$a cmp $b}keys%aa_single) {
            if ($tem[3] ne $kk) {
                $sift_muta{"$tem[3]$tem[4]$kk\n"} = 0;
                $polyphen_muta{"$key\t$tem[4]\t$tem[3]\t$kk\n"} = 0;
            }
        }
    }
    foreach my $kk(keys%sift_muta) {
        print OUT"$kk";
    }
    foreach my $kk(keys%polyphen_muta) {
        print OUTT"$kk";
    }
    close OUT;
    close OUTT;
    #run SIFT 
    system "$software_sift_path/bin/SIFT_for_submitting_fasta_seq.csh $seq_path/$key.fasta $uniref90_database_path/uniref90 $software_sift_input_path/$key.subst";
    move("$software_sift_path/tmp/$key.SIFTprediction", "$sift/$key.SIFTprediction");
    #run PolyPhen-2
    system "$software_polyphen_path/bin/run_pph.pl -s $seq_path/$key.fasta $software_polyphen_input_path/$key.txt 1>$polyphen/$key.output 2>$polyphen/$key.log";
    system "$software_polyphen_path/bin/run_weka.pl $polyphen/$key.output > $polyphen/$key.humdiv.output";
}


#=============================calculate structural features=============================

if (scalar(keys%pdb_avail_site) ne 0) {
    open(PPLI,">","$software_psaia_path/pdblist.txt") or die "$!";
    foreach my $key(keys%pdb_avail_site) {
        my ($pdb_name,$pdb_chain) = split(/\s+/,$key);
        print PPLI"$PDB_chain/$pdb_name-$pdb_chain.pdb\n";
        #calculates the pocket region
        system "$software_pocket_path -f $PDB_chain/$pdb_name-$pdb_chain.pdb";
        #calculates the secondary structure
        system "$software_dssp_path -i $PDB_chain/$pdb_name-$pdb_chain.pdb -o $secondary_structure/$pdb_name-$pdb_chain.dssp";
        &secondary_structure_type("$secondary_structure/$pdb_name-$pdb_chain.dssp");
        #calculates the hydrogen bonds
        system "$software_hbplus_path/hbplus $PDB_chain/$pdb_name-$pdb_chain.pdb";
        move("$pdb_name-$pdb_chain.hb2", "$hbplus/$pdb_name-$pdb_chain.hb2");
        &hbplus("$hbplus/$pdb_name-$pdb_chain.hb2");
        #calculates laplacian norm and bfactor
        &laplacian_norm("$PDB_chain/$pdb_name-$pdb_chain.pdb",$laplacian_norm,"$pdb_name-$pdb_chain",$bfactor,$ca);
        &zscore("$laplacian_norm/$pdb_name-$pdb_chain.txt",1);
        &zscore("$bfactor/$pdb_name-$pdb_chain.txt",1);
        #calculates topological features
        &residues_interaction_network("$PDB_chain/$pdb_name-$pdb_chain.pdb",$residue_network,"$pdb_name-$pdb_chain");
        system "python topological_features.py -net $residue_network/$pdb_name-$pdb_chain.txt -o $topological_features/$pdb_name-$pdb_chain.txt";
        &zscore_and_logistic("$topological_features/$pdb_name-$pdb_chain.txt",1);
    }
    close PPLI;
    #calculates the depth and protrusion index
    system "$software_psaia_path/psa $software_psaia_path/psa.cfg $software_psaia_path/pdblist.txt";
    opendir(DIR,"$software_psaia_path") or die "$!";
    my @dir = readdir(DIR);
    closedir DIR;
    foreach my $d(@dir) {
        if ($d =~ /\.tbl/) {
            my @name = split(/_/,$d);
            move("$software_psaia_path/$d", "$psaia/$name[0].tbl");
            &dp_index("$psaia/$name[0].tbl",$dp_index,$name[0]);
        }
    }
}

#=============================SEQ FEATURE=============================

open(FFP,">","$out_path/seqfeature.txt") or die "$!";
foreach my $key(keys%unique_seq) {
    foreach my $kk(@{$unique_seq{$key}}) {
        my @tem = split(/\s+/,$kk);
        print FFP "$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]\t";
        my $seq_distance = abs($tem[2] - $tem[4]);
        print FFP "$seq_distance\t";
        #co-evolution
        my $small;my $big;
        if ($tem[2] > $tem[4]) {
            $small = $tem[4];
            $big = $tem[2];
        }else{
            $big = $tem[4];
            $small = $tem[2];
        }
        open(IN,"<","$co_evolution/$key.txt.normalization") or die "$!";
   LINE:while (<IN>) {
            chomp;
            my @temp=split(/\s+/);
            if ($temp[0] eq $small and $temp[2] eq $big) {
                print FFP"$temp[5]\t";
                last LINE;
            }
        }
        close IN;
        #correlated mutation
        my ($target,$nontarget,$other) = &CM("$co_evolution/$key.msa",$kk);
        printf FFP "%0.4f\t%0.4f\t%0.4f\t",$target,$nontarget,$other;
        #domain
        open(IN,"<","$domain/$key.txt") or die "$!";
        chomp(my @domain = <IN>);
        close IN;
        my $domain_o1 = 0;
        foreach my $i(28..$#domain) {
            my @temp = split(/\s+/,$domain[$i]);
            if (($tem[2] >= $temp[3] and $tem[2] <= $temp[4]) and  ($tem[4] >= $temp[3] and $tem[4] <= $temp[4])) {
                $domain_o1 = 1;
            } 
        }
        print FFP "$domain_o1\t";
        #disorder
        open(IN,"<","$disorder/$key.txt") or die "$!";
        chomp(my @disorder = <IN>);
        close IN;
        my @region;
        foreach my $i(@disorder) {
            if ($i =~ /^>/) {
                my @temp = split(/\s+/,$i);
                foreach my $j(2..$#temp) {
                    $temp[$j] =~ s/,//;
                    push @region,$temp[$j];
                }
            }
        }
        my $site_region = 0;
        foreach my $i(@region) {
            my @temp = split(/-/,$i);
            if ($tem[2] >= $temp[0] and $tem[2] <= $temp[1] and $tem[4] >= $temp[0] and $tem[4] <= $temp[1]) {
                $site_region = 1;
            }
        }
        if ($site_region eq 1) {
            print FFP "1\t";
        }else{
            print FFP "0\t";
        }
        #sift
        my %sift_score;
        open(IN,"<","$sift/$key.SIFTprediction") or die "$!";
        while (<IN>) {
            chomp;
            if ($_ !~ /WARNING/) {
                my @temp = split(/\s+/);
                $temp[0] =~ s/[A-Z]//g;
                $sift_score{$temp[0]} += $temp[2]/19;
            }
        }
        close IN;
        my $sift_aver = ($sift_score{$tem[2]} + $sift_score{$tem[4]})/2;
        if ($sift_score{$tem[2]} > $sift_score{$tem[4]}) {
            print FFP "$sift_score{$tem[4]}\t$sift_score{$tem[2]}\t$sift_aver\t";
        }else{
            print FFP "$sift_score{$tem[2]}\t$sift_score{$tem[4]}\t$sift_aver\t";
        }
        #shannon entropy
        my %se_score;
        open(IN,"<","$shannon_entropy/$key.txt") or die "$!";
        while (<IN>) {
            chomp;
            my @temp = split(/\s+/);
            $se_score{$temp[0]} = $temp[2];
        }
        close IN;
        my $se_aver = ($se_score{$tem[2]} + $se_score{$tem[4]})/2;
        if ($se_score{$tem[2]} > $se_score{$tem[4]}) {
            print FFP "$se_score{$tem[4]}\t$se_score{$tem[2]}\t$se_aver\t";
        }else{
            print FFP "$se_score{$tem[2]}\t$se_score{$tem[4]}\t$se_aver\t";
        }
        #polyphen-2
        my %polyphen_score;
        open(IN,"<","$polyphen/$key.humdiv.output") or die "$!";
        while (<IN>) {
            if ($_ =~ /^$key/) {
                my @temp =split(/\t/,$_);
                $temp[1] =~ s/\s+//g;
                $polyphen_score{$temp[1]} += $temp[15]/19;
            }
        }
        close IN;
        my $polyphen_aver = ($polyphen_score{$tem[2]} + $polyphen_score{$tem[4]})/2;
        if ($polyphen_score{$tem[2]} > $polyphen_score{$tem[4]}) {
            print FFP "$polyphen_score{$tem[4]}\t$polyphen_score{$tem[2]}\t$polyphen_aver\n";
        }else{
            print FFP "$polyphen_score{$tem[2]}\t$polyphen_score{$tem[4]}\t$polyphen_aver\n";
        }
    }
}
close FFP;

#=============================STR FEATURE=============================

if (scalar(keys%pdb_avail_site) ne 0) {
    open(FFS,">","$out_path/strfeature.txt") or die "$!";
    foreach my $key(keys%pdb_avail_site) {
        my ($pdb_name,$pdb_chain) = split(/\s+/,$key);
        foreach my $kk(@{$pdb_avail_site{$key}}) {
            my @tem = split(/\s+/,$kk);
            print FFS "$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]\t";
            my $str_site1 = $tem[-1];
            my $str_site2 = $tem[-2];
            #structure distance
            my @two_coor;
            open(IN,"<","$ca/$pdb_name-$pdb_chain.pdb") or die "$!";
            while (<IN>) {
                chomp;
                my $num = substr($_,22,5);
                $num =~ s/\s+//g;
                if ($num eq $str_site1 or $num eq $str_site2) {
                    my $x=substr($_,30,8);
                    my $y=substr($_,38,8);
                    my $z=substr($_,46,8);
                    push @two_coor,$x;
                    push @two_coor,$y;
                    push @two_coor,$z;
                }
            }
            close IN;
            my $distance3d = sqrt(($two_coor[0] - $two_coor[3])**2 + ($two_coor[1] - $two_coor[4])**2 + ($two_coor[2] - $two_coor[5])**2);
            printf FFS "%0.4f\t",$distance3d;
            #the shortest path
            open(IN,"<","$residue_network/$pdb_name-$pdb_chain.txt") or die "$!";
            chomp(my @resnet = <IN>);
            close IN;
            my $step = 0;
            my @part;
            push @part,$str_site1;
            until($str_site2 ~~ @part or $step > 100) {
                @part = &part(\@part,\@resnet);
                $step++;
            }
            print FFS "$step\t";
            #pocket
            opendir(DIR,"$PDB_chain/$pdb_name-$pdb_chain"."_out/pockets");
            my @dir =readdir(DIR);
            closedir DIR;
            my %pocket_region;
            my $num_pocket = 0;
            foreach my $d(@dir) {
                if ($d =~ /atm/) {
                    my @bb;
                    open(IN,"<","$PDB_chain/$pdb_name-$pdb_chain"."_out/pockets/$d") or die "$!";
                    while (<IN>) {
                        chomp;
                        if ($_ =~ /^ATOM/) {
                            my $num = substr($_,22,5);
                            $num =~ s/\s+//g;
                            push @bb,$num;
                        }
                    }
                    close IN;
                    $num_pocket++;
                    $pocket_region{$num_pocket} = \@bb;
                }
            }
            if ($num_pocket eq 0) {
                print FFS "0\t";
            }else{
                my $same_pocket = 0;
                foreach my $i(keys%pocket_region) {
                    if ($str_site1 ~~ @{$pocket_region{$i}} and $str_site2 ~~ @{$pocket_region{$i}}) {
                        $same_pocket = 1;
                    }
                }
                print FFS "$same_pocket\t";
            }
            #secondary structure
            my %secondary_type;
            open(IN,"<","$secondary_structure/$pdb_name-$pdb_chain.dssp.type") or die "$!";
            while (<IN>) {
                chomp;
                my ($num,$type) = split(/\s+/);
                $secondary_type{$num} = $type;
            }
            close IN;
            my %hash;
            $hash{"CC"}=0;
            $hash{"HH"}=0;
            $hash{"SS"}=0;
            $hash{"CH"}=0;
            $hash{"CS"}=0;
            $hash{"HS"}=0;
            if (exists $hash{"$secondary_type{$str_site1}$secondary_type{$str_site2}"}) {
                $hash{"$secondary_type{$str_site1}$secondary_type{$str_site2}"}++;
            }elsif(exists $hash{"$secondary_type{$str_site2}$secondary_type{$str_site1}"}) {
                $hash{"$secondary_type{$str_site2}$secondary_type{$str_site1}"}++;
            }
            foreach my $i(sort{$a cmp $b}keys%hash) {
                print FFS "$hash{$i}\t";
            }
            #depth and protrusion index
            my @dp_score;my @cx_score;
            open(IN,"<","$dp_index/$pdb_name-$pdb_chain.txt") or die "$!";
            while (<IN>) {
                chomp;
                my @temp = split(/\s+/);
                if ($temp[0] eq $str_site1 or $temp[0] eq $str_site2) {
                    push @dp_score,$temp[1];
                    push @cx_score,$temp[2];
                }
            }
            close IN;
            @dp_score = sort{$a <=> $b}@dp_score;
            my $dp_aver = ($dp_score[0] + $dp_score[1])/2;
            foreach my $d(@dp_score) {
                print FFS "$d\t";
            }
            print FFS "$dp_aver\t";
            @cx_score = sort{$a <=> $b}@cx_score;
            my $cx_aver =($cx_score[0] + $cx_score[1])/2;
            foreach my $d(@cx_score) {
                print FFS "$d\t";
            }
            print FFS "$cx_aver\t";
            #topological feature
            my @tp_score;
            open(IN,"<","$topological_features/$pdb_name-$pdb_chain.txt.normalization") or die "$!";
            while (<IN>) {
                chomp;
                my @temp=split(/\s+/,$_);
                if ($temp[0] eq $str_site1 or $temp[0] eq $str_site2) {
                    foreach my $ou(1..$#temp) {
                        push @tp_score,$temp[$ou];
                    }
                }
            }
            close IN;
            foreach my $i(0..3) {
                if ($tp_score[$i] < $tp_score[$i+4]) {
                    print FFS"$tp_score[$i]\t$tp_score[$i+4]\t";
                }else{
                    print FFS"$tp_score[$i+4]\t$tp_score[$i]\t";
                }
                my $dpaver = ($tp_score[$i] + $tp_score[$i+4])/2;
                print FFS"$dpaver\t";
            }
            #laplacian norm
            my @laplacian_score;
            open(IN,"<","$laplacian_norm/$pdb_name-$pdb_chain.txt.normalization") or die "$!";
            while (<IN>) {
                chomp;
                my @temp=split(/\s+/,$_);
                if ($temp[0] eq $str_site1 or $temp[0] eq $str_site2) {
                    foreach my $ou(1..$#temp) {
                        push @laplacian_score,$temp[$ou];
                    }
                }
            }
            close IN;
            foreach my $i(0..4) {
                if ($laplacian_score[$i] < $laplacian_score[$i+5]) {
                    print FFS"$laplacian_score[$i]\t$laplacian_score[$i+5]\t";
                }else{
                    print FFS"$laplacian_score[$i+5]\t$laplacian_score[$i]\t";
                }
                my $laaver = ($laplacian_score[$i] + $laplacian_score[$i+5])/2;
                print FFS"$laaver\t";
            }
            #bfactor
            my @bfactor_score;
            open(IN,"<","$bfactor/$pdb_name-$pdb_chain.txt.normalization") or die "$!";
            while (<IN>) {
                chomp;
                my @temp = split(/\s+/);
                if ($temp[0] eq $str_site1 or $temp[0] eq $str_site2) {
                    push @bfactor_score,$temp[1];
                }
            }
            close IN;
            my $bfactor_aver = ($bfactor_score[0] + $bfactor_score[1])/2;
            if ($bfactor_score[0] > $bfactor_score[1]) {
                print FFS "$bfactor_score[0]\t$bfactor_score[1]\t$bfactor_aver\t";
            }else{
                print FFS "$bfactor_score[1]\t$bfactor_score[0]\t$bfactor_aver\t";
            }
            #hydrogen bonds
            my %hb_score;
            $hb_score{$str_site1} = 0;
            $hb_score{$str_site2} = 0;
            open(IN,"<","$hbplus/$pdb_name-$pdb_chain.hb2.txt") or die "$!";
            while (<IN>) {
                chomp;
                my @temp=split(/\s+/);
                if ($temp[0] eq $str_site1 or $temp[3] eq $str_site1) {
                    $hb_score{$str_site1}++;
                }
                if ($temp[0] eq $str_site2 or $temp[3] eq $str_site2) {
                    $hb_score{$str_site2}++;
                }
            }
            close IN;
            my $sum = 0;
            foreach my $i(sort{$hb_score{$a} <=> $hb_score{$b}}keys%hb_score) {
                print FFS "$hb_score{$i}\t";
                $sum += $hb_score{$i};
            }
            my $hbaver = $sum/2;
            print FFS "$hbaver\n";
        }
    }
    close FFS;
}
=cut
#====================PTM cross-talk prediction/sequence and structure====================

my $predict_sample = "$out_path";

my $model_path = "../model/";
my $seqsetname = "../random_set/seqset";
my $seqtest_path = "../random_set/seqtest";
my $strsetname = "../random_set/strset";
my $strtest_path = "../random_set/strtest";

&SEQSET;
&SEQFeature;
&SEQRF_python;
if (scalar(keys%pdb_avail_site) ne 0) {
    &STRSET;
    &STRFeature;
    &STRRF_python;
    &performance;
}else{
    &seqperformance;
}

#end
unlink "new_protein.txt";
unlink "hbdebug.dat";

#====================PTM cross-talk pairs prediction based on sequence features====================
sub SEQSET{
    open(IN,"<","$model_path/negative_seq_list.txt") or die "$!";
    my %hash_neg;
    my @neg;
    while (<IN>) {
        chomp;
        my @tem = split(/\s+/);
        $hash_neg{$tem[0]} = 0;
        push @neg,$tem[0];
    }
    close IN;
    for my $s(1..100) {
        open(OUT,">","$seqsetname/train$s.txt") or die "$!";
        open(IN,"<","$model_path/positive_seq_list.txt") or die "$!";
        while (<IN>) {
            chomp;
            print OUT"$_\t1\n";
        }
        close IN;
        close OUT;
    }   
    my $num = scalar@neg;
    for my $s(1..100) {
        open(OUT,">>","$seqsetname/train$s.txt") or die "$!";
        my %hash;
        while (scalar(keys%hash) < 205) {
            my $number = int(rand($num));
            $hash{$number} = 1;
        }
        foreach my $key(keys%hash) {
            print OUT"$neg[$key]\t$hash_neg{$neg[$key]}\n";
        }
        close OUT; 
    }
}
sub SEQFeature{
    my %feature;
    open(FF,"<","$model_path/positive_sequence_features.txt") or die "$!";
    while (<FF>) {
        chomp;
        my @tem = split(/\s+/);
        my @bb;
        #1,2,4,5,6,7,8,11,12,14,15
        foreach my $i(1,2,4,5,6,7,8,11,12,14,15) {
            push @bb,$tem[$i];
        }
        $feature{$tem[0]} = \@bb;
    }
    close FF;
    open(FF,"<","$model_path/negative_sequence_features.txt") or die "$!";
    while (<FF>) {
        chomp;
        my @tem = split(/\s+/);
        my @bb;
        foreach my $i(1,2,4,5,6,7,8,11,12,14,15) {
            push @bb,$tem[$i];
        }
        $feature{$tem[0]} = \@bb;
    }
    close FF;
    
    open(OUP,">","$seqtest_path/test.txt") or die "$!";
    open(PP,"<","$out_path/seqfeature.txt") or die "$!";
    while (<PP>) {
        chomp;
        my @temp = split(/\s+/);
        print OUP"$temp[0]\t";
        foreach my $i(1,2,4,5,6,7,8,11,12,14,15) {
            print OUP"$temp[$i]\t";
        }
        print OUP"\n";
    }
    close PP;
    close OUP;
    
    for my $ff(1..100) {
        open(IN,"<","$seqsetname/train$ff.txt") or die "$!";      
        open(OUT,">","$seqtest_path/train$ff.txt") or die "$!";
        while (<IN>) {
            chomp;
            my @temp = split(/\s+/,$_);
            print OUT"$temp[1]\t";
            foreach my $ou(@{$feature{$temp[0]}}) {
                print OUT"$ou\t";
            }
            print OUT"\n";
        }
        close IN;
        close OUT;
    } 
}
sub SEQRF_python{
    foreach my $i(1..100) {
        system "python python_RF.py -train $seqtest_path/train$i.txt -test $seqtest_path/test.txt -o $seqtest_path/$i"."_out.txt";
    }
}

#=======================PTM cross-talk pairs prediction based on structural features==================
sub STRSET{
    open(IN,"<","$model_path/negative_str_list.txt") or die "$!";
    my %hash_neg;
    my @neg;
    while (<IN>) {
        chomp;
        my @tem = split(/\s+/);
        $hash_neg{$tem[0]} = 0;
        push @neg,$tem[0];
    }
    close IN;
    for my $s(1..100) {
        open(OUT,">","$strsetname/train$s.txt") or die "$!";
        open(IN,"<","$model_path/positive_str_list.txt") or die "$!";
        while (<IN>) {
            chomp;
            print OUT"$_\t1\n";
        }
        close IN;
        close OUT;
    }   
    my $num = scalar@neg;
    for my $s(1..100) {
        open(OUT,">>","$strsetname/train$s.txt") or die "$!";
        my %hash;
        while (scalar(keys%hash) < 61) {
            my $number = int(rand($num));
            $hash{$number} = 1;
        }
        foreach my $key(keys%hash) {
            print OUT"$neg[$key]\t$hash_neg{$neg[$key]}\n";
        }
        close OUT; 
    }
}
sub STRFeature{
    my %feature;
    open(FF,"<","$model_path/positive_strutural_features.txt") or die "$!";
    while (<FF>) {
        chomp;
        my @tem = split(/\s+/);
        my @bb;
        #1,2,3,6,7,8,9,14,15,17,28,34,43,44,46
        foreach my $i(1,2,3,6,7,8,9,14,15,17,28,34,43,44,46) {
            push @bb,$tem[$i];
        }
        $feature{$tem[0]} = \@bb;
    }
    close FF;
    open(FF,"<","$model_path/negative_strutural_features.txt") or die "$!";
    while (<FF>) {
        chomp;
        my @tem = split(/\s+/);
        my @bb;
        foreach my $i(1,2,3,6,7,8,9,14,15,17,28,34,43,44,46) {
            push @bb,$tem[$i];
        }
        $feature{$tem[0]} = \@bb;
    }
    close FF;
    
    open(OUP,">","$strtest_path/test.txt") or die "$!";
    open(PP,"<","$out_path/strfeature.txt") or die "$!";
    while (<PP>) {
        chomp;
        my @temp = split(/\s+/);
        print OUP"$temp[0]\t";
        foreach my $i(1,2,3,6,7,8,9,14,15,17,28,34,43,44,46) {
            print OUP"$temp[$i]\t";
        }
        print OUP"\n";
    }
    close PP;
    close OUP;
    
    for my $ff(1..100) {
        open(IN,"<","$strsetname/train$ff.txt") or die "$!";      
        open(OUT,">","$strtest_path/train$ff.txt") or die "$!";
        while (<IN>) {
            chomp;
            my @temp = split(/\s+/,$_);
            print OUT"$temp[1]\t";
            foreach my $ou(@{$feature{$temp[0]}}) {
                print OUT"$ou\t";
            }
            print OUT"\n";
        }
        close IN;
        close OUT;
    } 
}
sub STRRF_python{
    foreach my $i(1..100) {
        system "python python_RF.py -train $strtest_path/train$i.txt -test $strtest_path/test.txt -o $strtest_path/$i"."_out.txt";
    }
}

#====================================probability of prediction======================================
sub performance{
    my %seq_pro;
    open(IN,"<","$seqtest_path/test.txt") or die "$!";
    chomp(my @seq = <IN>);
    close IN;
    foreach my $i(0..$#seq) {
        my @tem = split(/\s+/,$seq[$i]);
        foreach my $j(1..100) {
            open(INN,"<","$seqtest_path/$j"."_out.txt") or die "$!";
            chomp(my @seqout = <INN>);
            close INN;
            shift@seqout;
            my @temp = split(/\s+/,$seqout[$i]);
            $seq_pro{$tem[0]} += $temp[2]/100;
        }
    }
    my %str_pro;
    open(IN,"<","$strtest_path/test.txt") or die "$!";
    chomp(my @str = <IN>);
    close IN;
    foreach my $i(0..$#str) {
        my @tem = split(/\s+/,$str[$i]);
        foreach my $j(1..100) {
            open(INN,"<","$strtest_path/$j"."_out.txt") or die "$!";
            chomp(my @strout = <INN>);
            close INN;
            shift@strout;
            my @temp = split(/\s+/,$strout[$i]);
            $str_pro{$tem[0]} += $temp[2]/100;
        }
    }
    print "\n\npredicted results\n";
    open(TT,">","$out_path/predicted_results.txt") or die "$!";
    open(PP,"<","$pre_file") or die "$!";
    while (<PP>) {
        chomp;
        my @tem = split(/\s+/);
        print "$_\t";
        print "prediction_score\t";
        print TT"$_\t";
        print TT"prediction_score\t";
        if (exists $str_pro{"$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]"}) {
            my $percent = $str_pro{"$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]"} * 0.5 + $seq_pro{"$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]"} * 0.5;
            printf "%0.3f\n",$percent;
            printf TT"%0.3f\n",$percent;
        }else{
            printf "%0.3f\n",$seq_pro{"$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]"};
            printf TT "%0.3f\n",$seq_pro{"$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]"};
        }
    }
    close PP;
    close TT;
}

sub seqperformance{
    my $exact = shift@_;
    my %seq_pro;
    open(IN,"<","$seqtest_path/test.txt") or die "$!";
    chomp(my @seq = <IN>);
    close IN;
    foreach my $i(0..$#seq) {
        my @tem = split(/\s+/,$seq[$i]);
        foreach my $j(1..100) {
            open(INN,"<","$seqtest_path/$j"."_out.txt") or die "$!";
            chomp(my @seqout = <INN>);
            close INN;
            shift@seqout;
            my @temp = split(/\s+/,$seqout[$i]);
            $seq_pro{$tem[0]} += $temp[2]/100;
        }
    }
    print "\n\npredicted results\n";
    open(TT,">>","$out_path/predicted_results.txt") or die "$!";
    open(PP,"<","$pre_file") or die "$!";
    while (<PP>) {
        chomp;
        if ($_ =~ /$exact/) {
            my @tem = split(/\s+/);
            print "$_\t";
            print "prediction_score\t";
            print TT"$_\t";
            print TT"prediction_score\t";

            printf "%0.3f\n",$seq_pro{"$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]"};
            printf TT "%0.3f\n",$seq_pro{"$tem[0]_$tem[1]$tem[2]_$tem[3]$tem[4]"};
        }
    }
    close PP;
    close TT;
}
#=============================subprogram=============================

sub generate_MSA {
    my $seq_name = shift@_;
    my @msa;
    open(IN,"<","$seq_name") or die "$!";
    chomp(my @data = <IN>);
    close IN;
    my @coor;
    foreach my $i(0..$#data) {
        if ($data[$i] =~ /^Query_1/) {
            push @coor,$i;
        }
        if ($data[$i] =~ /^Lambda.*alpha$/) {
            push @coor,$i;
        }
    }
    foreach my $i($coor[0]..($coor[1]-3)) {
        push @msa,$data[$i];
    }
    my @line = split(/\s+/,$msa[0]);
    my @query = split(//,$line[2]);
    my $n = scalar@query;
    my @temp = split(//,$msa[0]);
    my $first;
    foreach my $i(3..$#temp) {
        if ($temp[$i] =~ /[A-Z]/) {
            $first = $i;
            last;
        }
    }
    $seq_name =~ s/\.psiblast/\.msa/;
    open(OUT,">","$seq_name") or die "$!";
    print OUT">$line[0]/1-$line[-1]\n";
    my @colscoor;
    foreach my $j(0..$#query) {
        if ($query[$j] =~ /-/) {
            push @colscoor,$j;
            $query[$j] =~ s/-//;
        }
    }
    my $str = join('',@query);
    print OUT"$str\n";
    foreach my $i(1..$#msa) {
        my $use = substr($msa[$i],$first,$n);
        my @array = split(//,$use);
        my @array1;
        foreach my $j(@array) {
            if ($j =~ /[A-Z]/ || $j =~ /-/) {
                push @array1,$j;
            }else{
                push @array1,'-';
            }
        }
        foreach my $k(@colscoor) {
            $array1[$k] =~ s/.*//;
        }
        my $out = join('',@array1);
        print OUT">align\n";
        print OUT"$out\n";
    }
    close OUT;
}

sub JSD {
    my($dir_in) = @_;
    my $chain_len = 0;
    my @chain;my @seq_num;my @aa_key;my @SE;
    my @residue_back = (0.078,0.051,0.041,0.052,0.024,0.034,0.059,0.083,0.025,0.062,0.092,0.056,0.024,0.044,0.043,0.059,0.055,0.014,0.034,0.072);
    open(IN,$dir_in)||die "Can not open file: $dir_in in the read_PSSM package\n";
    my @pssm_raw;
    while(my $raw = <IN>) {
        my @raw_wop;
        my $se;
        $raw =~ s/\n|\r//g;
        my @raw_inf = split(/ +/,$raw);
        my $num = scalar@raw_inf;
        my $t = 0;
        my $sum = 0;
        if($num > 42) {
            $chain[$chain_len] = $raw_inf[1]."\t".$raw_inf[2];
            $chain_len++;
            for(my $i = 0; $i < 20; $i++) {
                $raw_wop[$t] = $raw_inf[$i+23];
                $t++;
            }
            my @vne_matrix;
            for(my $i = 0;$i < 20; $i++) {
                $sum += $raw_wop[$i];
            }
            for(my $i = 0;$i < 20; $i++) {
                if($sum == 0 or $raw_wop[$i] == 0) {
                    $raw_wop[$i] = 0;
                    $se += 0;
                }else{
                    $raw_wop[$i] = $raw_wop[$i]/$sum;
                    $se += -$raw_wop[$i] * log($raw_wop[$i]);
                }			
            }
            push @seq_num,$raw_inf[1];
            push @aa_key,$raw_inf[2];
            push @SE,$se;
        }   
    }
    close IN;
    return (\@seq_num,\@aa_key,\@SE);
}

sub zscore_logistic {
    my @tempsub = @_;
    my $total = 0;
    foreach my $item (@tempsub) {
        $total += $item;
    }
    my $avgsub = $total/scalar@tempsub;
    
    my $sdsub = 0;
    foreach my $item (@tempsub){
        $sdsub += ($item-$avgsub)**2;
    }
    my $stdsub = sqrt($sdsub/($#tempsub+1));
    my @norm;
    foreach my $item (@tempsub) {
        my $norm = 0;
        if ($stdsub == 0) {
            $norm = 0;
        }else{
            $norm = ($item-$avgsub)/$stdsub;
            $norm = sprintf"%.4f",1/(1+exp(-$norm));
        }
        push(@norm,$norm);
    }
    return \@norm;
}

sub secondary_structure_type{
    my $dssp_path = shift@_;
    open(IN,"<","$dssp_path") or die "$!";
    open(OUT,">","$dssp_path.type") or die "$!";
    chomp(my @dssp=<IN>);
    close IN;
    foreach my $i(@dssp) {
        if ($i =~ /\d$/) {
            my $num = substr($i,5,6);
            $num =~ s/\s+//g;
            my $type = substr($i,16,1);
            if ($type eq 'H' or $type eq 'G' or $type eq 'I') {
                print OUT"$num\tH\n";
            }elsif($type eq 'E' or $type eq 'B') {
                print OUT"$num\tS\n";
            }else{
                print OUT"$num\tC\n";
            }   
        }   
    }
    close OUT;
}

sub hbplus{
    my $input_file = shift@_;
    open(OUT,">","$input_file.txt") or die "$!";
    open(IN,"<","$input_file") or die "$!";
    chomp(my @data=<IN>);
    close IN;
    my %hash;
    foreach my $i(8..$#data) {
        my $num1;
        my $num2;
        my @tem=split(/\s+/,$data[$i]);
        my @temp=split(/-/,$tem[0]);
        if (scalar@temp eq 2) {
            $temp[0] =~ s/[A-Z]//;
            if ($temp[0] eq '0000') {
                $temp[0] = 0;
            }else{
                $temp[0] =~ s/^0*//g; 
            }
            print OUT"$temp[0]\t$temp[1]\t$tem[1]\t";
            $num1 = $temp[0];
        }else{
            $temp[1] =~ s/^0+//;
            print OUT"-$temp[1]\t$temp[2]\t$tem[1]\t";
            $num1 = "-$temp[1]";
        }

        my @temp1=split(/-/,$tem[2]);
        if (scalar@temp1 eq 2) {
            $temp1[0] =~ s/[A-Z]//;
            if ($temp1[0] eq '0000') {
                $temp1[0] = 0;
                }else{
                $temp1[0] =~ s/^0*//g;
            }
            print OUT"$temp1[0]\t$temp1[1]\t$tem[3]\n";
            $num2 = $temp1[0];
        }else{
            $temp1[1] =~ s/^0+//;
            print OUT"-$temp1[1]\t$temp1[2]\t$tem[3]\n";
            $num2 = "-$temp1[1]";
        }
        if ($num1 > $num2) {
            ($num1,$num2) = ($num2,$num1);
        }
        $hash{$num1}{$num2}+=1;
    }
    close OUT;
}

sub laplacian_norm{
    my $pdb_path = shift@_;
    my $output_path = shift@_;
    my $out_name = shift@_;
    my $bfactor = shift@_;
    my $ca = shift@_;

    my %pdb;
    open(IN,"<","$pdb_path") or die "$!";
    while (<IN>) {
        chomp;
        if ($_ =~ /^ATOM/) {
            my $num = substr($_,22,5);
            $num =~ s/\s+//g;
            my $atom=substr($_,13,2);
            if ($atom eq 'CA') {
                if (exists $pdb{$num}) {
                    push @{$pdb{$num}},$_;
                }else{
                    my @bb;
                    push @bb,$_;
                    $pdb{$num} = \@bb;
                }
            }
        }
    }
    close IN;
    open(BB,">","$bfactor/$out_name.txt") or die "$!";
    open(CA,">","$ca/$out_name.pdb") or die "$!";
    my @pdb_ca;
    foreach my $key(sort{$a <=> $b}keys%pdb) {
        my $factor = substr($pdb{$key}[0],61,6);
        $factor =~ s/\s+//g;
        print BB"$key\t$factor\n";
        push @pdb_ca,$pdb{$key}[0];
        print CA "$pdb{$key}[0]\n";
    }
    close BB;
    close CA;

    open(OUT,">","$output_path/$out_name.txt") or die "$!";
    my @all_distance;
    foreach my $i(@pdb_ca) {
        my $num=substr($i,22,5);
        $num =~ s/\s+//g;
        my $x=substr($i,30,8);
        my $y=substr($i,38,8);
        my $z=substr($i,46,8);
        foreach my $j(@pdb_ca) {
            my $num1=substr($j,22,5);
            $num1 =~ s/\s+//g;
            if ($num gt $num1) {
                my $x1=substr($j,30,8);
                my $y1=substr($j,38,8);
                my $z1=substr($j,46,8);
                my $dis=sqrt(($x-$x1)**2+($y-$y1)**2+($z-$z1)**2);
                push @all_distance,$dis;
            }
        }
    }
    @all_distance = sort{$a <=> $b}@all_distance;
    my $mean=int(scalar@all_distance/4);
    my $mean1=$mean * 2;
    my $mean2=$mean * 3;
    my %hash;
    foreach my $i(@pdb_ca) {
        my $num=substr($i,22,5);
        $num =~ s/\s+//g;
        my $x=substr($i,30,8);
        my $y=substr($i,38,8);
        my $z=substr($i,46,8);
        foreach my $j(@pdb_ca) {
            my $num1=substr($j,22,5);
            $num1 =~ s/\s+//g;
            my $ff=$num;
            my $dd=$num1;
            $ff =~ s/[A-Z]//;
            $dd =~ s/[A-Z]//;
            if (abs($ff - $dd) > 1) {
                my $x1=substr($j,30,8);
                my $y1=substr($j,38,8);
                my $z1=substr($j,46,8);
                my $dis=(($x-$x1)**2+($y-$y1)**2+($z-$z1)**2);
                foreach my $cd(0,$mean,$mean1,$mean2,$#all_distance) {
                    my $cda = $all_distance[$cd];
                    my $ou;
                    if ($cda ne 0) {
                        $ou=exp(-($dis/($cda**2)));
                    }else{
                        $ou=1;
                    }
                    $hash{$cda}{$num}{$num1}=$ou;
                }
            }
        }
    }
    foreach my $i(@pdb_ca) {
        my $num=substr($i,22,5);
        $num =~ s/\s+//g;
        my $x=substr($i,30,8);
        my $y=substr($i,38,8);
        my $z=substr($i,46,8);
        print OUT"$num\t";
        foreach my $cd(0,$mean,$mean1,$mean2,$#all_distance) {
            my $cda = $all_distance[$cd];
            my $sumx=0;my $sumy=0;my $sumz=0;
            my $sumb=0;
            foreach my $j(@pdb_ca) {
                my $num1=substr($j,22,5);
                $num1 =~ s/\s+//g;
                my $ff=$num;
                my $dd=$num1;
                $ff =~ s/[A-Z]//;
                $dd =~ s/[A-Z]//;
                if (abs($ff - $dd) > 1) {
                    my $x1=substr($j,30,8);
                    my $y1=substr($j,38,8);
                    my $z1=substr($j,46,8);
                    my $aa= ($x1 * $hash{$cda}{$num}{$num1});
                    my $bb= ($y1 * $hash{$cda}{$num}{$num1});
                    my $cc= ($z1 * $hash{$cda}{$num}{$num1});
                    $sumx += $aa;
                    $sumy += $bb;
                    $sumz += $cc;
                    $sumb += $hash{$cda}{$num}{$num1};
                }
            }
            my $meanx=$sumx/$sumb;
            my $meany=$sumy/$sumb;
            my $meanz=$sumz/$sumb;
            my $distance=sqrt(($x-$meanx)**2 + ($y-$meany)**2 + ($z-$meanz)**2);
            print OUT"$distance\t";
        }
        print OUT"\n";
    }
    close OUT;
}

sub residues_interaction_network{
    my $pdb_path = shift@_;
    my $output_path = shift@_;
    my $out_name = shift@_;

    open(IN,"<","$pdb_path") or die "$!";
    my %coor;
    while (<IN>) {
        if ($_ =~ /^ATOM/) {
            chomp;
            my $res=substr($_,22,5);
            my $atom=substr($_,6,5);
            $res =~ s/\s+//g;
            $atom =~ s/\s+//g;
            my $x=substr($_,30,8);
            my $y=substr($_,38,8);
            my $z=substr($_,46,8);
            my @bb;
            push @bb,$x;
            push @bb,$y;
            push @bb,$z;
            $coor{$res}{$atom}=\@bb;
        }
    }
    close IN;
    open(OUT,">","$output_path/$out_name.txt") or die "$!";
    my %hash;
    foreach my $key(sort{$a <=> $b}keys%coor) {
        foreach my $key1(sort{$a <=> $b}keys%coor) {
            if ($key1 gt $key) {
                foreach my $key2(keys%{$coor{$key}}) {
                    foreach my $key3(keys%{$coor{$key1}}) {
                        my $dis=sqrt(($coor{$key}{$key2}[0]-$coor{$key1}{$key3}[0])**2 + ($coor{$key}{$key2}[1]-$coor{$key1}{$key3}[1])**2 + ($coor{$key}{$key2}[2]-$coor{$key1}{$key3}[2])**2);
                        if ($dis <= 5) {
                            $hash{"$key\t$key1\n"}=0;
                        }
                    }
                }
            }
        }
    }
    foreach my $kk(sort{$a cmp $b}keys%hash) {
        print OUT"$kk";
    }
    close OUT;
}

sub zscore_and_logistic{
    my $path_file = shift@_;
    my $stat = shift@_;

    open(IN,"<","$path_file") or die "$!";
    chomp(my @data=<IN>);
    close IN;

    open(OUT,">","$path_file.normalization") or die "$!";
    my @tem=split(/\s+/,$data[0]);
    my @mean;
    my @sd;
    for (my $i=$stat;$i <= $#tem;$i++) {
        my @aa;
        foreach my $j(@data) {
            my @temp=split(/\s+/,$j);
            push @aa,$temp[$i];
        }
        my $a = Statistics::Descriptive::Full->new();
        $a->add_data(\@aa);
        my $mean_a = $a->mean();
        my $standard_deviation_a=$a->standard_deviation();
        push @mean,$mean_a;
        push @sd,$standard_deviation_a;
    }

    foreach my $i(@data) {
        my @temp=split(/\s+/,$i);
        foreach my $p(0..($stat-1)) {
            print OUT"$temp[$p]\t";
        }
        foreach my $j($stat..$#temp) {
            my $k=$j-$stat;
            if($sd[$k] ne 0) {
                my $out=($temp[$j]-$mean[$k])/$sd[$k];
                my $aa=1/(1+exp(-$out));
                printf OUT"%0.4f\t",$aa;
            }else{
                my $aa=1/(1+exp(-$temp[$j]));
                printf OUT"%0.4f\t",$aa;
            }   
        }
        print OUT"\n";
    }
    close OUT;
}

sub dp_index{
    my $dp_input = shift@_;
    my $output_path = shift@_;
    my $out_name = shift@_;

    open(IN,"<","$dp_input") or die "$!";
    chomp(my @data = <IN>);
    close IN;

    open(OUT,">","$output_path/$out_name.txt") or die "$!";
    foreach my $i(@data) {
        if ($i =~ /^\s+.+\d$/) {
            my @array = split(/\s+/,$i);
            $array[2] =~ s/\s+//g;
            print OUT"$array[2]\t";
            my $dp = 1/(1+exp(-$array[4]));
            my $cx = 1/(1+exp(-$array[10]));
            printf OUT"%0.4f\t%0.4f\n",$dp,$cx;
        }
    }
    close OUT;
}

#calculate correlated mutation
sub CM{
    my @path = @_;
    my @tem =split(/\s+/,$path[1]);
    my %hash;
    open(IN,"<","$path[0]") or die "$!";
    while (<IN>) {
        chomp;
        if ($_ !~ /^>/) {
            my @str=split(//);
            my $site1 = $tem[2];
            my $site2 = $tem[4];
            if ($str[$site1-1]."$site1" eq "$tem[1]$tem[2]" and $str[$site2-1]."$site2" eq "$tem[3]$tem[4]") {
                $hash{"A"}++;
            }elsif($str[$site1-1]."$site1" ne "$tem[1]$tem[2]" and $str[$site2-1]."$site2" ne "$tem[3]$tem[4]") {
                $hash{"B"}++;
            }else{
                $hash{"C"}++;
            } 
        }
    }
    close IN;
    if (!exists $hash{"A"}) {
        $hash{"A"}=0;
    }
    if (!exists $hash{"B"}) {
        $hash{"B"}=0;
    }
    if (!exists $hash{"C"}) {
        $hash{"C"}=0;
    }
    my $sum=$hash{"A"} + $hash{"B"} + $hash{"C"};   
    my $pera=$hash{"A"}/$sum;
    my $perb=$hash{"B"}/$sum;
    my $perc=$hash{"C"}/$sum;
    return $pera,$perb,$perc;
}

#calculate the shortest path
sub part{
    my ($ref_a, $ref_b) = @_;
    my @node = @{$ref_a};
    my @dat = @{$ref_b};
    my %aa;
    foreach my $i(@dat) {
        my @tem=split(/\s+/,$i);
        foreach my $j(@node) {
            if ($tem[0] eq $j) {
                $aa{$tem[1]}=0;
            }
            if ($tem[1] eq $j) {
                $aa{$tem[0]}=0;
            }
        }
    }
    my @key=keys%aa;
    return @key;
}

sub zscore {
    
    my $path_file = shift@_;
    my $stat = shift@_;

    open(IN,"<","$path_file") or die "$!";
    chomp(my @data=<IN>);
    close IN;

    open(OUT,">","$path_file.normalization") or die "$!";
    my @tem=split(/\s+/,$data[0]);
    my @mean;
    my @sd;
    for (my $i=$stat;$i <= $#tem;$i++) {
        my @aa;
        foreach my $j(@data) {
            my @temp=split(/\s+/,$j);
            push @aa,$temp[$i];
        }
        my $a = Statistics::Descriptive::Full->new();
        $a->add_data(\@aa);
        my $mean_a = $a->mean();
        my $standard_deviation_a=$a->standard_deviation();
        push @mean,$mean_a;
        push @sd,$standard_deviation_a;
    }

    foreach my $i(@data) {
        my @temp=split(/\s+/,$i);
        foreach my $p(0..($stat-1)) {
            print OUT"$temp[$p]\t";
        }
        foreach my $j($stat..$#temp) {
            my $k=$j-$stat;
            if($sd[$k] ne 0) {
                my $out=($temp[$j]-$mean[$k])/$sd[$k];
                printf OUT "%0.4f\t",$out;
            }else{
                printf OUT "%0.4f\t",$temp[$j];
            }   
        }
        print OUT"\n";
    }
    close OUT;
}

sub print_help {
    printf "PCRpred USAGE:
perl PCRpred.pl -f prediction_file -seq /the/directory/of/sequence/ -str /the/directory/of/PDB/ -out /the/directory/of/output/
such as:
perl PCRpred.pl -f ../examples/prefile.txt -seq ../sequences/ -str ../PDB/ -out ../examples/
==================
arguments:
==================
    -f          the path of prediction file, such as ../examples/prefile.txt
    -seq        the path of sequence file (.fasta)
    -str        the path of structure file (.pdb)
    -out        the path of prediction result\n";
}
