#!/usr/bin/perl

# Wrapper script for building feature tracks fetched from UCSC for JBrowse
#
# Author : Panagiotis Moulos (pmoulos@hybridstat.gr)

use strict;
use boolean;
use File::Basename;
use File::Copy;
use File::Fetch;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile tempdir);
use Getopt::Long;
use Pod::Usage;

use Data::Dumper;
             
# Make sure output is unbuffered
select(STDOUT);
$|=1;

# Cleanup code
$SIG{INT} = \&catch_cleanup;

# Set defaults
our $scriptname = "buildFeatureTracks.pl";
our $config_file;
our $genomes;
our $features;
our $update;
our $ncores = 1;
our $man;
our $log;
our $silent;
our $help;

# The storage format will be $track_store_home/reference/$organism/seq and
# $track_store_home/reference/$organism/features

&check_inputs;

# Test non-standard modules
&try_module("JSON");
&try_module("JSON::Parse","json_file_to_perl","valid_json");
&try_module("Parallel::Loops");

our $logfilehandle = &open_log_file if ($log); # Log file if requested
my $date = &now;
disp("$date - Started...\n");

disp("Reading configuration file $config_file");
our $config = json_file_to_perl($config_file);
our $tmpdir = tempdir(CLEANUP => 1);

disp("Building requested tracks");
if ($genomes) {
    disp("Genomes");
    &build_dir_structure("genome");
    #&fetch("genome");
    &create_tracks("genome");
}
if ($features) {
    disp("Features");
    &build_dir_structure("feature");
    &fetch("feature");
    &create_tracks("feature");
}
if ($update) {
    &update_feature_track_lists;
}

$date = &now;
disp("\n$date - Finished!\n");
&close_log_file if ($log); # Close the log file if requested


sub build_dir_structure {
    my $what = $_[0];
    my $home_dir = $config->{"home"};
    disp("  Creating $what directory structure if not present...");
    my ($err,$reference_track_dir);
    my $reference_dir = File::Spec->catdir($home_dir,"reference");
    foreach my $genome (@{$config->{"genomes"}}) {
        if ($what eq "genome") {
            $reference_track_dir = 
                File::Spec->catfile($reference_dir,
                    $genome->{"unique_name"},"seq");
        }
        elsif ($what eq "feature") {
            $reference_track_dir = 
                File::Spec->catfile($reference_dir,
                    $genome->{"unique_name"},"tracks");
        }
        if (! -d $reference_track_dir) {
            make_path($reference_track_dir,{error => \$err});
            if (@$err) {
                die "\nCould not create directory $reference_track_dir! ".
                    "Check path spelling and permissions...\n";
            }
        }
    }
    return(0);
}

sub fetch {
    my $what = $_[0];
    #&fetch_genomes if ($what eq "genome");
    if ($what eq "feature") {
        &fetch_features;
        &fetch_trackdbs;
    }
}

sub create_tracks {
    my $what = $_[0];
    &create_genome_tracks if ($what eq "genome");
    &create_feature_tracks if ($what eq "feature");
}

#sub fetch_genomes {
#    my $genome;
#    my %fetch_hash;
#    disp("  Fetching genomes");
#    foreach $genome (@{$config->{"genomes"}}) {
#        %{$fetch_hash{$genome->{"unique_name"}}} = (
#            "name" => $genome->{"unique_name"},
#            "organism" => $genome->{"friendly_name"},
#            "version" => $genome->{"assembly_version"},
#            "source" => $genome->{"source"},
#            "fetch_path" => $genome->{"path_to_fasta"}
#        );
#    }
#    if ($ncores == 1) {
#        foreach $genome (keys(%fetch_hash)) {
#            &fetch_genome($fetch_hash{$genome});
#        }
#    }
#    else {
#        my @genomes = keys(%fetch_hash);
#        my $pl_fg = Parallel::Loops->new($ncores);
#        $pl_fg->share(\%fetch_hash);
#        $pl_fg->foreach(\@genomes,sub {
#            &fetch_genome($fetch_hash{$_});
#        });
#    }
#}

sub fetch_features {
    my ($feature,$table);
    my %fetch_hash;
    my $counter = 0;
    disp("  Fetching features");
    foreach $feature (@{$config->{"features"}}) {
        foreach $table (@{$feature->{"tables"}}) {
            $counter++;
            %{$fetch_hash{$counter}} = (
                "genome" => $feature->{"genome"},
                "table" => $table->{"name"},
                "fetch_path" => $feature->{"path_to_db"}
            );
        }
    }
    if ($ncores == 1) {
        foreach my $key (keys(%fetch_hash)) {
            &fetch_feature($fetch_hash{$key});
        }
    }
    else {
        my @keys = keys(%fetch_hash);
        my $pl_ff = Parallel::Loops->new($ncores);
        $pl_ff->share(\%fetch_hash);
        $pl_ff->foreach(\@keys,sub {
            &fetch_feature($fetch_hash{$_});
        });
    }
}

sub fetch_trackdbs {
    my %fetch_hash;
    my ($feature,$genome);
    disp("  Fetching UCSC trackDb files");
    foreach $feature (@{$config->{"features"}}) {
        %{$fetch_hash{$feature->{"genome"}}} = (
            "genome" => $feature->{"genome"},
            "fetch_path" => $feature->{"path_to_db"}
        );
    }
    if ($ncores == 1) {
        foreach $genome (keys(%fetch_hash)) {
            &fetch_trackdb($fetch_hash{$genome});
        }
    }
    else {
        my @genomes = keys(%fetch_hash);
        my $pl_ft = Parallel::Loops->new($ncores);
        $pl_ft->share(\%fetch_hash);
        $pl_ft->foreach(\@genomes,sub {
            &fetch_trackdb($fetch_hash{$_});
        });
    }
}

sub create_genome_tracks {
    my $genome;
    my %create_hash;
    disp("  Creating genome tracks");
    foreach $genome (@{$config->{"genomes"}}) {
        %{$create_hash{$genome->{"unique_name"}}} = (
            "name" => $genome->{"unique_name"},
            "full_name" => $genome->{"name"},
            "organism" => $genome->{"friendly_name"},
            "version" => $genome->{"assembly_version"},
            "source" => $genome->{"source"},
            "dir" => File::Spec->catfile($config->{"home"},"reference",
                $genome->{"unique_name"})
        );
    }
    if ($ncores == 1) {
        foreach $genome (keys(%create_hash)) {
            &create_genome_track($create_hash{$genome});
        }
    }
    else {
        my @genomes = keys(%create_hash);
        my $pl_cg = Parallel::Loops->new($ncores);
        $pl_cg->share(\%create_hash);
        $pl_cg->foreach(\@genomes,sub {
            &create_genome_track($create_hash{$_});
        });
    }
}

sub create_feature_tracks {
    my ($feature,$table);
    my %create_hash;
    my $counter = 0;
    disp("  Creating feature tracks");
    foreach $feature (@{$config->{"features"}}) {
        foreach $table (@{$feature->{"tables"}}) {
            $counter++;
            %{$create_hash{$counter}} = (
                "genome" => $feature->{"genome"},
                "table" => $table->{"name"},
                "dir" => File::Spec->catfile($config->{"home"},"reference",
                    $feature->{"genome"})
            );
        }
    }
    # That's not like fetching... Each execution adds to a JSON file... Parallel
    # execution might mess it up as they are not split per genome
    foreach my $key ((keys(%create_hash))) {
        &create_feature_track($create_hash{$key});
    }
    #if ($ncores == 1) {
    #    foreach (my $key (keys(%create_hash))) {
    #       &create_feature_track($create_hash{$key});
    #   }
    #}
    #else {
    #    my @keys = keys(%create_hash);
    #    my $pl_cf = Parallel::Loops->new($ncores);
    #    $pl_cf->share(\%create_hash);
    #    $pl_cf->foreach(\@keys,sub {
    #        &create_feature_track($fetch_hash{$_});
    #    });
    #}
}

sub update_feature_track_lists {
    my ($track_list_file,$track_list_dir);
    my %table_info;
    disp("  Updating feature tracks details");
    foreach my $flist (@{$config->{"features"}}) {
        foreach my $table (@{$flist->{"tables"}}) {
            %{$table_info{$flist->{"genome"}}{$table->{"name"}}} = (
                "label" => $table->{"label"},
                "description" => $table->{"description"}
            );
        }
    }
    foreach my $genome (@{$config->{"genomes"}}) {
        my $gen = $genome->{"unique_name"};
        disp("    ".$gen);
        $track_list_dir = File::Spec->catdir($config->{"home"},"reference",
            $gen);
        $track_list_file = File::Spec->catfile($track_list_dir,"trackList.json");
        my $track_list = json_file_to_perl($track_list_file);
        for (my $i=0;$i<scalar @{$track_list->{"tracks"}};$i++) {
            if ($track_list->{"tracks"}->[$i]->{"type"} =~ m/FeatureTrack/) {
                $track_list->{"tracks"}->[$i]->{"description"} = 
                    $table_info{$gen}{$track_list->{"tracks"}->[$i]->{"key"}}{"description"};
                $track_list->{"tracks"}->[$i]->{"key"} = 
                    $table_info{$gen}{$track_list->{"tracks"}->[$i]->{"key"}}{"label"};
                $track_list->{"tracks"}->[$i]->{"urlTemplate"} = 
                    $config->{"url_base"}."/reference/$gen/".
                        $track_list->{"tracks"}->[$i]->{"urlTemplate"};
            }
            if ($track_list->{"tracks"}->[$i]->{"type"} =~ m/SequenceTrack/) {
                $track_list->{"tracks"}->[$i]->{"urlTemplate"} = 
                    $config->{"url_base"}."/reference/$gen/".
                        $track_list->{"tracks"}->[$i]->{"urlTemplate"};
            }
        }
        open(TRACKLIST,">$track_list_file");
        my $json = JSON->new;
        my $json_text = $json->pretty->encode($track_list);
        print TRACKLIST $json_text;
        #print TRACKLIST encode_json($track_list);
        close(TRACKLIST);
        
        # ...and generate the index of feature names for autocompletion
        disp("  Generating feature names for searching...");
        my $generate_script = File::Spec->catfile($config->{"jbrowse_home"},
            "bin","generate-names.pl");
        my @arg_list = ();
        push(@arg_list,"--out ".$track_list_dir);
        #push(@arg_list,"--compress");
        push(@arg_list,"--verbose");
        my $generate_command = "perl ".$generate_script." ".join(" ",@arg_list);
        disp("    Will now run");
        disp("-----------------------------------------------------------------");
        disp($generate_command);
        disp("-----------------------------------------------------------------\n");
        my $success = system($generate_command);
    }
}

#sub fetch_genome {
#    my %genome_metadata = %{$_[0]};
#    my $gname = $genome_metadata{"name"};
#    my $gorganism = $genome_metadata{"organism"};
#    my $gversion = $genome_metadata{"version"};
#    my $gsource = $genome_metadata{"source"};
#    my $gpath = $genome_metadata{"fetch_path"};
#    $gpath =~ s/http/rsync/i if ($gpath =~ m/^(http|ftp)/);
#    
#    disp("    Fetching genome $gname\n");
#    my $tmp_gpath = File::Spec->catfile($tmpdir,$gname);
#    my $rsync = File::Rsync->new({
#        "archive" => 1,
#        "compress" => 1,
#        "partial" => 1,
#        "recursive" => 1,
#        "progress" => 1,
#        "stats" => 1,
#        "verbose" => 1,
#        "human-readable" => 1,
#        "outfun" => \&outnow,
#        "errfun" => \&errnow
#    });
#    $rsync->exec({
#        src => $gpath,
#        dest => $tmp_gpath,
#    }) or die "\n".join("",@{$rsync->err})."\n";
#}

sub fetch_feature {
    my %feature_metadata = %{$_[0]};
    my $fgenome = $feature_metadata{"genome"};
    my $ftable = $feature_metadata{"table"};
    my $gpath = $feature_metadata{"fetch_path"};
    
    disp("    Fetching feature $ftable from $fgenome");
    my $tmp_fpath = File::Spec->catfile($tmpdir,$fgenome);
    my $fdb = $gpath.$ftable.".sql";
    my $ftx = $gpath.$ftable.".txt.gz";
    
    if ($gpath !~ m/^(http|ftp|rsync|file|git)/) {
        copy($fdb,$tmp_fpath);
        copy($ftx,$tmp_fpath);
    }
    else {
        my $ffdb = File::Fetch->new("uri" => $fdb);
        my $fftx = File::Fetch->new("uri" => $ftx);
        my $wheredb = $ffdb->fetch("to" => $tmp_fpath) or die;
        my $wheretx = $fftx->fetch("to" => $tmp_fpath) or die;
    }
}

sub fetch_trackdb {
    my %tdb_metadata = %{$_[0]};
    my $tgenome = $tdb_metadata{"genome"};
    my $tpath = $tdb_metadata{"fetch_path"};
    
    disp("    Fetching trackDb fike for $tgenome");
    my $tmp_tpath = File::Spec->catfile($tmpdir,$tgenome);
    my $tdb = $tpath."trackDb.sql";
    my $ttx = $tpath."trackDb.txt.gz";
    if ($tpath !~ m/^(http|ftp|rsync|file|git)/) {
        copy($tdb,$tmp_tpath);
        copy($ttx,$tmp_tpath);
    }
    else {
        my $ffdb = File::Fetch->new("uri" => $tdb);
        my $fftx = File::Fetch->new("uri" => $ttx);
        my $wheredb = $ffdb->fetch("to" => $tmp_tpath) or die;
        my $wheretx = $fftx->fetch("to" => $tmp_tpath) or die;
    }
}
    
sub create_genome_track {
    my %genome_metadata = %{$_[0]};
    my $gname = $genome_metadata{"name"};
    my $gfullname = $genome_metadata{"full_name"};
    my $gorganism = $genome_metadata{"organism"};
    my $gversion = $genome_metadata{"version"};
    my $gsource = $genome_metadata{"source"};
    my $ghome = $genome_metadata{"dir"};
    
    disp("    Creating genome track for $gname");
    my $prepare_script = 
        File::Spec->catfile($config->{"jbrowse_home"},"bin",
            "prepare-refseqs.pl");
    #opendir(my $dfh,$tmp_gpath);
    #my @refseqs_fa = grep { !/(^\.\.?$)|(\.txt$)/ } readdir($dfh);
    #closedir($dfh);
    #for (my $i=0;$i<scalar @refseqs_fa;$i++) {
    #    $refseqs_fa[$i] = File::Spec->catfile($tmp_gpath,$refseqs_fa[$i]);
    #}
    #my @arg_list = @refseqs_fa;
    #for (my $i=0;$i<scalar @arg_list;$i++) {
    #    $arg_list[$i] = "--fasta $refseqs_fa[$i]";
    #}
    my $seq_file = File::Spec->catfile($config->{"igenomes_base"},$gfullname,
        $gsource,$gname,"Sequence","WholeGenomeFasta","genome.fa");
    my @arg_list = ();
    push(@arg_list,"--fasta ".$seq_file);
    push(@arg_list,"--key \"".$gsource." ".$gorganism." ".$gversion."/".
        $gname."\"");
    push(@arg_list,"--trackLabel \"".lc($gsource)."_".lc($gorganism)."_".
        lc($gname)."\"");
    push(@arg_list,"--out ".$ghome);
    #push(@arg_list,"--compress");
    my $prepare_command = "perl ".$prepare_script." ".join(" ",@arg_list);
    disp("    Will now run");
    disp("-----------------------------------------------------------------");
    disp($prepare_command);
    disp("-----------------------------------------------------------------\n");
    my $success = system($prepare_command);
    return ($success);
}

sub create_feature_track {
    my %feature_metadata = %{$_[0]};
    my $fgenome = $feature_metadata{"genome"};
    my $ftable = $feature_metadata{"table"};
    my $fhome = $feature_metadata{"dir"};
    
    disp("    Creating feature track for $ftable of $fgenome");
    my $tmp_fpath = File::Spec->catfile($tmpdir,$fgenome);
    my $convert_script = 
        File::Spec->catfile($config->{"jbrowse_home"},"bin","ucsc-to-json.pl");
    my @arg_list = (
        "--in ".$tmp_fpath,
        "--out ".$fhome,
        "--track ".$ftable,
        "--primaryName name2"#,
        #"--compress"
    );
    my $convert_command = "perl ".$convert_script." ".join(" ",@arg_list);
    disp("    Will now run");
    disp("-----------------------------------------------------------------");
    disp($convert_command);
    disp("-----------------------------------------------------------------\n");
    my $success = system($convert_command);
    return($success);
}

# Process inputs
sub check_inputs {
    my $stop;
    GetOptions(
        "config|c=s" => \$config_file,
        "genomes|g" => \$genomes,
        "features|f" => \$features,
        "update|u" => \$update,
        "ncores|p=i" => \$ncores,
        "man|m" => \$man,
        "log|l" => \$log,
        "silent|s" => \$silent,
        "help|h" => \$help
    );
    
    Pod::Usage::pod2usage( -verbose => 1, -exitstatus => 0 ) if ($help);
    Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ($man);
    
    $stop .= "--- Please specify a configuration JSON file ---\n" 
        if (!$config_file);
    $stop .= "--- Please specify at least one of --genomes,--features or ".
        "--update ---\n" if (!$genomes && !$features && !$update);
    if ($stop) {
            print "\n$stop\n";
            print "Type perl $scriptname --help for help in usage.\n\n";
            exit;
    }
}

sub open_log_file {
    my $logname = "hs_buildReferenceTracks_".&now("machine").".log";
    my $lfh = IO::File->new();
    $lfh->open(">$logname");
    return($lfh);
}

sub close_log_file { 
    $logfilehandle->close; 
}

sub now {
    my $format = shift @_;
    $format = "human" if (!$format);
    my ($sec,$min,$hour,$day,$month,$year) = localtime(time);
    $year += 1900;
    $month++;
    $month = "0".$month if (length($month)==1);
    $day = "0".$day if (length($day)==1);
    $hour = "0".$hour if (length($hour)==1);
    $min = "0".$min if (length($min)==1);
    $sec = "0".$sec if (length($sec)==1);
    ($format ne "machine") ? (return($day."/".$month."/".$year." ".$hour.":".
        $min.":".$sec)) :
    (return($year.$month.$day.$hour.$min.$sec));
}

sub try_module {
    my ($module,@fun) = @_;
    eval "require $module";
    if ($@) {
        my $killer = "Module $module is required to continue with the\n". 
            "execution. If you are in Windows and you have ActiveState Perl\n".
            "installed, use the Package Manager to get the module. If you \n".
            "are under Linux, log in as a super user (or use sudo under\n".
            "Ubuntu) and type \"perl -MCPAN -e shell\" (you will possibly\n".
            "have to answer some questions). After this type \"install \n".
            "$module\" to install the module. If you don't know how to\n".
            "install the module, contact your system administrator.\n";
        die "\n$killer\n\n";
    }
    else {
        if (@fun) {
            my $funs = join(" ",@fun);
            eval "use $module qw($funs)";
        }
        else { eval "use $module"; }
    }
}

sub disp {
    print STDERR "\n@_" if (!$silent);
    print $logfilehandle "\n@_" if (!$silent && $log);
}

sub catch_cleanup {
    print STDERR "\nCatching Ctrl-C, cleaning temporary files!\n";
    &cleanup;
    die;
}

sub cleanup {
    remove_tree($tmpdir);
}

sub outnow {
    print STDERR $_[0];
    print $logfilehandle $_[0] if ($log);
}

sub errnow {
    print STDERR $_[0];
    print $logfilehandle $_[0] if ($log);
}


__END__

=pod

=head1 NAME

buildReferenceTracks.pl - Perl script to build reference tracks for the
JBrowse version of SeqCVIBE

=head1 SYNOPSIS

hs_buildReferenceTracks.pl --config config_json_file [OPTIONS]

Examples:

=over 4

=item perl hs_buildReferenceTracks.pl --config conf.json --genomes

=item perl hs_buildReferenceTracks.pl --config conf.json --features --ncores 4

=back

=head1 DESCRIPTION

buildReferenceTracks.pl is a Perl script to automate the procedure of creating
reference tracks for the JBrowse installation in SeqCVIBE. It accepts a 
configuration file with the (mostly self-explanatory) parameters which
set up browser locations etc. This script is not yet smart enough to download
genome fasta files from UCSC (some are broken in chromosomes etc.). So, for the
time being it requires a local Illumina iGenomes installation and the 
installation directory must be passed to the config JSON file. In the future, we
will most probably compile and maintain a SeqCVIBE genomes collation, directly
usable by JBrowse.

=head1 ARGUMENTS

=over 4

=item config B<(required)>

--config or -c

A valid configuration JSON file.

=item genomes B<(optional if --features specified, required otherwise)>

--genomes or -g

Install genomes from a local Illumina iGenomes installation. Optional, but must 
be run at least once to create the JBrowse reference sequences.

=item features B<(optional if --genomes specified, required otherwise)>

--features or -f

Download and install feature tracks from UCSC. Optional, but must be run at 
least once to create the JBrowse reference tracks (UCSC/RefSeq genes etc.).

=item update B<(optional)>

--update or -u

Update the JBrowse trackList.json file. Optional but if not enabled, new tracks
will not become visible.

=item ncores B<(optional)>

--ncores or -p

Runs the script in parallel mode where one core processes one bedgraph file. It
requires the module Parallel::Loops to be installed.

=item silent B<(optional)>

--silent or -s

Do not display verbose messages.

=item help B<(optional)>

--help or -h

Display this help text.

=item man B<(optional)>

--man or -m

Display the full manual of the script.

=back

=head1 AUTHOR

Panagiotis Moulos (L<pmoulos@hybridstat.gr)>

=cut
