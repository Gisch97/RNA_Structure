#!/usr/bin/env perl
use strict;
use warnings;
use Text::CSV_XS;
use JSON;
use File::Spec;

# Archivos de entrada/salida
my $infile  = 'ArchiveII.csv';
my $outfile = 'ArchiveII_with_motiv.csv';

# Inicializamos CSV
my $csv_in  = Text::CSV_XS->new({ binary=>1, auto_diag=>1 });
my $csv_out = Text::CSV_XS->new({ binary=>1, eol=>"\n" });

open my $fh_in,  '<:encoding(utf8)', $infile  or die "No puedo abrir '$infile': $!";
open my $fh_out, '>:encoding(utf8)', $outfile or die "No puedo crear '$outfile': $!";

# Leemos y extendemos el encabezado
my $header = $csv_in->getline($fh_in);
push @$header, 'motivos';
$csv_out->print($fh_out, $header);

my $json = JSON->new->allow_nonref;

while (my $row = $csv_in->getline($fh_in)) {
    my ($id, $seq, $struct, $bp_json, $len) = @$row;

    # 1) Parsear los pares de bases
    my $pairs = $json->decode($bp_json);    # [ [l,r], [l2,r2], ... ]
    my %bp = map { ($_->[0] => $_->[1], $_->[1] => $_->[0]) } @$pairs;

    # 2) Generar archivo .bpseq
    my $bpseq_file = "$id.bpseq";
    open my $out_bp, '>:encoding(utf8)', $bpseq_file 
      or die "No puedo escribir '$bpseq_file': $!";
    for my $i (1 .. length($seq)) {
        my $base = substr($seq, $i-1, 1);
        my $j    = $bp{$i} // 0;
        print $out_bp "$i $base $j\n";
    }
    close $out_bp;

    # 3) Ejecutar bpRNA.pl
    system('perl', 'bpRNA.pl', $bpseq_file) == 0
      or warn "Error procesando '$bpseq_file': $!";

    # 4) Leer el .st resultante y extraer la línea de labels
    my $st_file = "$id.st";
    open my $in_st, '<:encoding(utf8)', $st_file 
      or die "No puedo abrir '$st_file': $!";
    my $motivos;
    my @non_comment;
    while (<$in_st>) {
        next if /^#/;               # salto cabeceras y warnings
        chomp;
        push @non_comment, $_;
        last if @non_comment >= 3;  # con 3 líneas me basta
    }
    close $in_st;
    # @non_comment = ( seq, dotbracket, labels )
    $motivos = $non_comment[2] // '';

    # 5) Añadimos motivos al row y lo escribimos
    push @$row, $motivos;
    $csv_out->print($fh_out, $row);
}

close $fh_in;
close $fh_out;

print "He generado '$outfile' con la columna adicional 'motivos'.\n";
