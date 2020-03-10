#!/usr/bin/perl


## libname statt @read , done
## stderr input reads, unique reads, 
## count matrices 5' und 3' ende


use strict;


my %len;
for(my $i=15;$i<77;$i++){
        $len{$i}{'raw'}=0;
        $len{$i}{'uniq'}=0;
}



## get correct file name -> there is a perl module for this but better we do not use it
my $library=shift @ARGV or die "no library name given\n";


my $total_reads=0;
my $uniq_reads=0;

my $type='fastq';

my $seq;
my %hash;
my %h2;

my %pwm3=();
my %pwm5=();

my $rlen=20;

if($type eq 'fastq'){
        my $line=0;
        while(<>){
                chomp;
                next if(/^\s*$/);
                $line++;
                if($line == 1){
                }elsif($line == 2){
                        $seq=uc($_);
                        $hash{$seq}++;
                        $len{length($seq)}{'raw'}++;
                        $total_reads++;

                }elsif($line == 3){
                }elsif($line == 4){
                        if(not $h2{$seq}){
                                $h2{$seq}=$_;
                        }
                        $line=0;
                }
        }
        close IN;
        my $count=0;
        for(sort {$hash{$b} <=> $hash{$a}} keys %hash){
                print "\@${library}_${count}_x$hash{$_}\n$_\n+\n$h2{$_}\n";
                $len{length($_)}{'uniq'}++;
                $count++;
                pwm(\%pwm5,\%pwm3,\$_,\$hash{$_});
        }
        $uniq_reads=$count;
}




## print out stats to stderr
print STDERR ">read statistics
#property\tcounts_$library
input\t$total_reads
uniq\t$uniq_reads

";


my @nts=qw(A C G T );

for(my $i=0;$i < $rlen ;$i++){
        ## add pseudocount
        foreach my $nt(@nts){
                $pwm5{$i}{$nt}++;
                $pwm3{$i}{$nt}++;
                $pwm5{$i}{'tot'}++;
                $pwm3{$i}{'tot'}++;
        }
}

print STDERR ">nt satistics from 5'end\n## --col_stack=A,C,G,T\n#pos";
foreach my $nt(@nts){
        print STDERR "\t$nt";
}

print STDERR "\n";

for(my $i=0;$i < 20;$i++){
        print STDERR "$i";
        foreach my $nt(@nts){
                printf( STDERR  "\t%.3f",$pwm5{$i}{$nt}/$pwm5{$i}{'tot'});
#               printf( STDERR  "\t%.2f\t%i\t%.2f",$pwm5{$i}{$nt},$pwm5{$i}{'tot'},$pwm5{$i}{$nt}/$pwm5{$i}{'tot'});
        }
        print STDERR "\n";
}



print STDERR "\n>nt satistics from 3'end\n## --col_stack=A,C,G,T\n#pos";
foreach my $nt(@nts){
        print STDERR "\t$nt";
}

print STDERR "\n";

for(my $i=0;$i < 20;$i++){
        print STDERR "-$i";
        foreach my $nt(@nts){
                printf( STDERR  "\t%.3f",$pwm3{$i}{$nt}/$pwm3{$i}{'tot'});
        }
        print STDERR "\n";
}


print STDERR "\n>read_lengths_raw\n";
print STDERR "## --cumsum --steps\n";
print STDERR "#length\t$library\n";
for my $i(sort {$a <=> $b} keys %len){
        print STDERR "$i\t$len{$i}{'raw'}\n";
}



print STDERR "\n>read_lengths_unique\n";
print STDERR "## --cumsum --steps\n";
print STDERR "#length\t$library\n";
for my $i(sort {$a <=> $b} keys %len){
        print STDERR "$i\t$len{$i}{'uniq'}\n";
}




## pwm function
sub pwm{
        my ($h5,$h3,$seq,$num) = @_;

        my $c;
        my $len=length($$seq);
        for(my $i=0; $i<$len;$i++){
                $c=uc(substr($$seq,$i,1));
                
                # this is for unique read sequences
                $$h5{$i}{$c}++;
                $$h3{$len-$i-1}{$c}++;
                $$h5{$i}{'tot'}++;
                $$h3{$len-$i-1}{'tot'}++;

                # this would be weighted
                # $$h5{$i}{$c}+=$$num;
                # $$h3{$len-$i-1}{$c}+=$$num;
                # $$h5{$i}{'tot'}+=$$num;
                # $$h3{$len-$i-1}{'tot'}+=$$num;
        }
}

