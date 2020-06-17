SeqCVIBE deployment
================================================================================

# Prerequisites

In order to deploy the SeqCVIBE Shiny application and for all its 
functionalities to be available (e.g. the Genome Browser) there are some
prerequisites that must be covered in the installation environment, whether this
is a bare metal server or a VM. The following prerequisites do not cover the
data analysis procedure not the data organization in the filesystem, for the
former to be automatically recognizable by SeqCVIBE. Although this is also a
straightforward procedure, it is covered elsewhere.

Thus, the following are assumed for the deployment of SeqCVIBE application:

#### Basic understanding of

 * Writing Linux shell commands to run scripts
 * R language
 * Apache web server

#### A recent version of R with the following R/Bioconductor packages installed

 * AnnotationDbi
 * AnnotationFilter
 * BiocFileCache
 * BiocManager
 * biovizBase
 * BSgenome
 * GenomicRanges
 * jsonlite
 * colourpicker
 * data.table
 * DESeq
 * devtools
 * DT
 * EDASeq
 * edgeR
 * genefilter
 * geneplotter
 * GenomicAlignments
 * GenomicFeatures
 * GenomeInfoDb
 * ggbio
 * ggplot2
 * heatmaply
 * httr
 * magrittr
 * metaseqR
 * openssl
 * OrganismDbi
 * plyr
 * RCurl
 * Rhtslib
 * rmarkdown
 * rmdformats
 * Rsamtools
 * RSQLite
 * rtracklayer
 * shiny
 * shinyBS
 * shinyjs
 * SummarizedExperiment
 * XML

#### A recent version of Perl with the following packages installed

 * boolean
 * JSON
 * JSON::Parse
 * Parallel::Loops
  
#### [Illumina iGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html)

For the time being, the Genome Browser support building process requires a
local version of the genomes you want to set up. This is provided to a config
file used with a Perl build script. After the build process, they are no longer
required.

#### Shiny server 

You will need a recent version of the free tier of RStudio
[Shiny server](https://rstudio.com/products/shiny/download-server/)

The following sections cover each part above as well as the main deployment
process.

## Setting up R and the required R packages

### System packages

Whether the server/VM ships with R or not, some system packages must be present:

```
sudo apt install apt-transport-https software-properties-common \
    build-essential zlib1g-dev libdb-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev apache2 
```

### Install R

We suppose that our server/VM does not ship with R included. If this is not the
case just skip to the next step.

```
sudo apt-key adv --keyserver keyserver.ubuntu.com \
    --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository \
    'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
sudo apt install r-base

# R --version -> 3.6.1 OK
```

### Install the required R packages

```
sudo R
```

Then, within R:

```
install.packages("BiocManager")

library(BiocManager)

BiocManage::install()
```

The above will install basic Bioconductor packages that are required for most
other ones. After this and still within R:

```
pkgs <- c("biomaRt","GenomicRanges","jsonlite","colourpicker","data.table",
    "DESeq","devtools","DT","EDASeq","edgeR","GenomicAlignments","GenomeInfoDb",
    "ggplot2","ggbio","heatmaply","magrittr","metaseqR","openssl","plyr",
    "rmarkdown","rmdformats","Rsamtools","RSQLite","rtracklayer","shiny",
    "shinyBS","shinyjs","XML")

BiocManager::install(pkgs)
```

## Setting up the required Perl packages

If you are using a new server/VM, the first time that ```CPAN``` is used to 
install Perl packages, it will ask or try to perform some configurations. The
following should do it:

```
sudo perl -MCPAN -e shell
```

Then install the required packages:

```
sudo perl -MCPAN -e 'install qw(
    boolean
    JSON
    JSON::Parse
    inc::latest
    Parallel::Loops
)'
```

## Shiny server

We are using the official instructions [here](https://rstudio.com/products/shiny/download-server/ubuntu/)

```
sudo apt-get install gdebi-core
wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.12.933-amd64.deb
sudo gdebi shiny-server-1.5.12.933-amd64.deb
```

## Illumina iGenomes

Next, we have to download the required Illumina iGenomes for the organisms we
want to have Genome Browser support for. For this guide, we will download only
the human genome version ```hg19``` and mouse genome version ```mm10```. 
Illumina iGenomes carry also a lot of annotation data and other sequences so it
may take a bit to download and extract. All these data can be removed after the
genome browser reference sequences are built.

So, supposing that our main storage location is in ```/media/storage```

```
cd /media/storage
mkdir igenomes
cd igenomes

# Human
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz
tar -vxf Homo_sapiens_UCSC_hg19.tar.gz

# Mouse
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -xvf Mus_musculus_UCSC_mm10.tar.gz

# After extraction delete the archives
rm Homo_sapiens_UCSC_hg19.tar.gz Mus_musculus_UCSC_mm10.tar.gz
```

## Data transfer

If there are already analyzed data that can be explored through SeqCVIBE, they
can be transfered to the VM using simple SSH transfer with ```scp```. The 
directory will have to be declared in the global SeqCVIBE configuration file.
For this guide, we suppose that the SeqCVIBE data will live in 
```/media/storage/pilot```.

# Setting up SeqCVIBE

Firstly, clone the repository in a specific subdirectory

```
mkdir elixir
cd elixir
git clone https://github.com/hybridstat/elixir-RNAseq.git
cd elixir-RNAseq

SEQC_HOME=/media/storage/elixir/elixir-RNAseq
```

We are now ready to proceed with the following sections which cover the 
deployment of the main application.

## Directory structure and configuration file

We suppose that the data transfered in a previous step are in the proper format
and have the proper directory structure required by SeqCVIBE. We will need this
path for the configuration file which is constructed later.

## Metadata database

During (or after) the data directory and structure construction, an SQLite 
database is constructed. This database must be trasnfered in ```config/``` and
called ```metadata.sqlite```, that is ```config/metadata.sqlite```.

## Genome browser

The Genome Browser adopted by SeqCVIBE is [JBrowse](http://jbrowse.org/). 
JBrowse is a genome browser written entirely in JavaScript and is fully 
embeddable, therefore fully suited for the SeqCVIBE application. However, it
cannot be served by Shiny server, therefore we will have to use Apache and 
```mod_proxy```. For details, see the next section **Web server setup**.

Generally, JBrowse is quite stable and has not changed much during the past few 
years. Therefore, a version downloaded today is fairly safe. JBrowse works with
configuration files and hard-coded instructions to display tracks that live in 
JSON files. Therefore, a specific structure is needed to make it work. This 
guide is designed to make this work reproducible and we have written extra 
scripts that will render the deployment much easier.

## JBrowse installation

1. Create the directory that JBrowse libraries will live in

```
mkdir $SEQC_HOME/jbrowse
```

2. Download the latest JBrowse release.

```
cd $SEQC_HOME/jbrowse
curl -L -O https://github.com/GMOD/jbrowse/releases/download/1.16.6-release/JBrowse-1.16.6.zip
unzip JBrowse-1.16.6.zip
rm JBrowse-1.16.6.zip
cd JBrowse-1.16.6
mv * ../
cd ..
rm -r JBrowse-1.16.6
./setup.sh
```

The installation will take a while depending on your system. What is mainly
installed is several Perl modules including BioPerl modules.

## Genome browser tracks setup

JBrowse supports several kinds of track types through related (static) 
configuration files. In SeqCVIBE we are using four track types:

* Sequence tracks (```SequenceTrack```) which represent the genomic sequences
themselves
* Feature tracks (```FeatureTrack```) which depicts genomic features like genes
* Signal tracks (```JBrowse/View/Track/Wiggle/XYPlot```) which are the same as
bigWig tracks
* Signal density tracks (```JBrowse/View/Track/Wiggle/Density```) which are a
single-dimension representation of bigWig files

The first two types are setup once (may be updated but less frequently). The
other two types are setup once but may be more frequently updated (adding,
removing etc.)

Generally, you will have to decide (also according to the hosting architecture)
where the reference sequences and user data tracks will live and whether they
will be placed in the same location or not. What we propose is to put the
reference tracks in a different location than the user tracks because:

* Reference tracks do not change often and they can be easily reconstructed
* User tracks change more often and keeping them in a different location may
make backup and security easier
* The logic is better

**Note** that all track directories must be served from a web-server and viewed
from the web (whether globally or locally).

### Reference tracks

The reference tracks are downloaded and formatted using a JSON configuration
file and the Perl script ```buildReferenceTracks.pl``` which is located in the
```lib``` directory of the repository root. For an example configuration file,
see the example at the ```config``` directory of the repository root. Most
options are self-explanatory. After you decide where the reference tracks will
live, put the path in the ```home``` option of the configuration file.

Then, the ```buildReferenceTracks.pl``` script must be executed twice:

* Once to format the sequence tracks

```
perl buildReferenceTracks.pl --config /path/to/config/config.json --genomes
```

* Once to format the feature tracks

```
perl buildReferenceTracks.pl --config /path/to/config/config.json --features
```

* Once to update the URL data of the feature tracks and build name indexes

```
perl buildReferenceTracks.pl --config /path/to/config/config.json --update
```

The ```buildReferenceTracks.pl``` script is essentially a wrapper around other
scripts provided by JBrowse to format the data in such a way that the browser
understands.

Thus in our case:

1. Create the directory that will host the genome browser reference tracks

```
mkdir -p /media/storage/pilot/tracks
```

2. Create the reference track building JSON file in 
```$SEQC_HOME/config/config_jbrowse_build.json```. This looks like the 
following:

```
{
    "home": "/media/storage/pilot/tracks",
    "jbrowse_home": "/media/storage/elixir/elixir-RNAseq/jbrowse",
    "url_base": "http://62.217.82.158/seqc_elixir",
    "igenomes_base": "/media/storage/igenomes",
    "genomes": [{
        "name": "Homo_sapiens",
        "friendly_name": "Human",
        "unique_name": "hg19",
        "assembly_version": "GRCh37",
        "source": "UCSC",
        "description": "Human genome, version GRCh37 (hg19), Feb 2009"
    },{
        "name": "Mus_musculus",
        "friendly_name": "Mouse",
        "unique_name": "mm10",
        "assembly_version": "GRCm38",
        "source": "UCSC",
        "description": "Mouse genome, version GRCm38 (mm10), Dec 2011"
    }],
    "features": [{
        "genome": "hg19",
        "path_to_db": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/",
        "tables": [{
            "name": "knownGene",
            "label": "UCSC genes",
            "description": "UCSC genes/transcripts/isoforms"
        },{
            "name": "refGene",
            "label": "RefSeq genes",
            "description": "RefSeq genes/transcripts/isoforms"
        },{
            "name": "ensGene",
            "label": "Ensembl genes",
            "description": "Ensembl genes/transcripts/isoforms"
        },{
            "name": "acembly",
            "label": "AceView genes",
            "description": "AceView assembly genes/transcripts/isoforms"
        }]
    },{
        "genome": "mm10",
        "path_to_db": "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/",
        "tables": [{
            "name": "knownGene",
            "label": "UCSC genes",
            "description": "UCSC genes/transcripts/isoforms"
        },{
            "name": "refGene",
            "label": "RefSeq genes",
            "description": "RefSeq genes/transcripts/isoforms"
        },{
            "name": "ensGene",
            "label": "Ensembl genes",
            "description": "Ensembl genes/transcripts/isoforms"
        }]
    }]
}
```

The configuration file parameters are mostly self-explanatory. 

3. Run the reference sequences building script

```
cd $SEQC_HOME/lib
perl buildReferenceTracks.pl --config ../config/config_jbrowse_build.json \
    --genomes
```

4. Run the reference features building script

```
perl buildReferenceTracks.pl --config ../config/config_jbrowse_build.json \
    --features
```

5. Update track locations and build name indexes

```
perl buildReferenceTracks.pl --config ../config/config_jbrowse_build.json \
    --update
```

### Data tracks

The location of the data tracks is defined by the main directory structure of
the data (BAM files, BigWig files, ```.rda``` files etc.) tied to SeqCVIBE. More
specifically, the location of the data tracks depends on:

* The main data home, into which subdirectories follow a certain structure, see
here (TODO)
* The names in the metadata SQLite database which define sources, conditions and
samples

Given the above, the ```R``` function ```buildTrackList``` located in 
```lib/build.R``` will do all the job given four arguments:

* ```metadata``` which is the path to the SQLite metadata database.
* ```annotationPath``` which is the path to the reference JBrowse tracks.
* ```urlBase``` which is the URL base to appended to each user data track, this
is an important one but can be changed later in case of re-deployment. The
destination of the URL base must point to the actual data directory (through a
symbolic link for example) and the latter must be servable by a web-browser (see
**Web server setup** below).
* ```appBase``` which is the SeqCVIBE application home (so that the function 
can write the necessary symlinks to the browser tracks).

Thus, an indicative run for building tracks for the first time would be:

```
Rscript trackControl.R \
    --operation=build \
    --metadata=... \
    --annotation_path=... \
    --url_base=... \
    --app_base=...
```

In our case, the above command should look like the following:

```
cd $SEQC_HOME/lib
Rscript trackControl.R \
    --operation=build \
    --metadata=../config/metadata.sqlite \
    --annotation_path=/media/storage/pilot/tracks/reference \
    --url_base=http://62.217.82.158/seqc_elixir \
    --app_base=/media/storage/elixir/elixir-RNAseq
```

or from within R:

```
source("/media/storage/elixir/elixir-RNAseq/lib/build.R")
source("/media/storage/elixir/elixir-RNAseq/lib/util.R")
buildTrackList(
    annoPath="/media/storage/pilot/tracks/reference",
    urlBase="http://62.217.82.158/seqc_elixir",
    appBase="/media/storage/elixir/elixir-RNAseq",
    metadata="/media/storage/elixir/elixir-RNAseq/config/metadata.sqlite"
)
```

Finally, because of the nature of SeqCVIBE data (users etc.), faceted tracks
must be enabled. To do this, go to ```$SEQC_HOME/jbrowse/jbrowse.conf``` and
uncomment (as per the instructions in the file itself) lines 10 and 11:

```
[trackSelector]
type = Faceted
```

## Main configuration file

Now that JBrowse and the data tracks are in place, we should create the main
SeqCVIBE JSON configuration file in ```$SEQC_HOME/config```:

```
cd $SEQC_HOME/config
nano app_config.json
```

and then paste the following (for our case):

```
{
    "meta": {
        "app_name": "SeqCVIBE",
        "app_alias": "seqcvibe"
    },
    "paths": {
        "data": "/media/storage/pilot",
        "metadata": "config/metadata.sqlite"
    },
    "urls": {
        "browser": "http://62.217.82.158/seqc_browse/index.html?",
        "tracks": "http://62.217.82.158/seqc_elixir/reference"
    }
}
```

A generic application configuration can be found in 
```$SEQC_HOME/config/app_config_template.json```

## Move the testing app into production

Now we are ready to move the application in the directory server by Shiny Server

```
sudo mkdir /srv/shiny-server/seqcvibe
cd $SEQC_HOME
sudo cp -r * /srv/shiny-server/seqcvibe/
sudo rm /srv/shiny-server/seqcvibe/README.md
```

## Web server setup

In this section we cover the apache configuration to redirect SeqCVIBE calls to
Shiny server and also serve JBrowse.

### Preliminaries

1. Enable the ```mod_deflate```, ```mod_headers```, ```mod_proxy``` and 
```proxy_http``` modules.

```
sudo a2enmod deflate
sudo a2enmod headers
sudo a2enmod proxy
sudo a2enmod proxy_http
sudo a2enmod proxy_wstunnel
```

2. A symbolic link must then be created in the directory served by the web 
server, for example, in our case, to support JBrowse and the tracks:

```
sudo ln -s /media/storage/pilot/ /var/www/seqc_elixir
```

### Apache Virtual Host

Below, the configuration required for the instance that will run our SeqCVIBE
instance together with JBrowse. Some general variables below have to be replaced
like ```$HOME```. See the table below.

|Variable/Attribute|Value                                      |  
|------------------|-------------------------------------------|
|$HOME             |The home of SeqCVIBE deployment            |
|ServerAdmin       |e-mail address of the server administrator |
|ServerName        |The deployment server address/name         |
|ServerAlias       |The deployment server alias                |

Replace $HOME with the SeqCVIBE root directory below.

```
# Virtual host for SeqcVIBE jbrowse instance
<VirtualHost *:80>
    ServerAdmin seqc_admin@seqc_domain.com
    DocumentRoot /var/www
    ServerName seqcvibe.server.name
    ServerAlias AAA.BBB.CCC.ZZZ

    # General mod_deflate (compression) settings (comment to remove compression)
    <IfModule mod_deflate.c>
        SetOutputFilter DEFLATE
        AddOutputFilterByType DEFLATE text/html text/xml text/css text/plain
        AddOutputFilterByType DEFLATE image/svg+xml application/xhtml+xml application/xml
        AddOutputFilterByType DEFLATE application/rdf+xml application/rss+xml application/atom+xml
        AddOutputFilterByType DEFLATE text/javascript application/javascript application/x-javascript application/json
        AddOutputFilterByType DEFLATE application/x-font-ttf application/x-font-otf
        AddOutputFilterByType DEFLATE font/truetype font/opentype
        BrowserMatch ^Mozilla/4 gzip-only-text/html
        BrowserMatch ^Mozilla/4\.0[678] no-gzip
        BrowserMatch \bMSIE !no-gzip !gzip-only-text/html
        BrowserMatch \bMSI[E] !no-gzip !gzip-only-text/html
        SetEnvIf Request_URI (\.jsonz|\.txtz|\.txt\.gz|\.vcf\.gz)$ no-gzip dont-vary
        SetEnvIfNoCase Request_URI \.(?:gif|jpe?g|png)$ no-gzip
        SetEnvIfNoCase Request_URI (\.bw|\.bigwig|\.bb|\.bigbed|\.bam|\.bai|\.tbi)$ no-gzip dont-vary
        Header append Vary User-Agent env=!dont-vary
        DeflateBufferSize 130023424
    </IfModule>

    <Directory />
        Require all denied
        Options FollowSymLinks
        AllowOverride None
    </Directory>

    <Directory /var/www>
        Options Indexes FollowSymLinks MultiViews
        AllowOverride None
        Require all denied
    </Directory>

    <Directory $HOME/tracks>
        Header set Access-Control-Allow-Origin *
        Header onsuccess set Access-Control-Allow-Headers X-Requested-With,Range
        Options Indexes FollowSymLinks MultiViews
        AllowOverride All
        Require all granted
    </Directory>

    <Directory $HOME/jbrowse>
        Header onsuccess set Access-Control-Allow-Origin *
        Header onsuccess set Access-Control-Allow-Headers X-Requested-With,Range
        Options Indexes FollowSymLinks MultiViews
        AllowOverride AuthConfig FileInfo
        Require all granted
    </Directory>
    
    <Proxy *>
        Allow from localhost
    </Proxy>
    
    ProxyPreserveHost On
    ProxyPass /seqcvibe http://seqcvibe.server.name:3838/seqcvibe
    ProxyPassReverse /seqcvibe http://seqcvibe.server.name:3838/seqcvibe
    ## If you want to do this through localhost - more guaranteed result
    #ProxyPass /seqcvibe http://localhost:3838/seqcvibe
    #ProxyPassReverse /seqcvibe http://localhost:3838/seqcvibe

    RewriteEngine on
    RewriteCond %{HTTP:Upgrade} =websocket
    RewriteRule /(.*) ws://localhost:3838/$1 [P,L]
    ## The following lines are suggested by Shiny developers but in our case
    ## it brakes JBrowse
    #RewriteCond %{HTTP:Upgrade} !=websocket
    #RewriteRule /(.*) http://localhost:3838/$1 [P,L]
    ## If you host only Shiny apps
    #ProxyPass / http://localhost:3838/
    #ProxyPassReverse / http://localhost:3838/
    ProxyRequests Off
    
    ServerSignature Off
</VirtualHost>
```

In our case, this should look like the following which should be written in
```/etc/apache2/sites-available```:

```
sudo nano /etc/apache2/sites-available/elixir-seqcvibe.conf
```

and then paste

```
# Virtual host for SeqcVIBE jbrowse instance
<VirtualHost *:80>
    ServerAdmin moulos@fleming.gr
    DocumentRoot /var/www
    ServerName 62.217.82.158
    ServerAlias 62.217.82.158

    # General mod_deflate (compression) settings (comment to remove compression)
    <IfModule mod_deflate.c>
        SetOutputFilter DEFLATE
        AddOutputFilterByType DEFLATE text/html text/xml text/css text/plain
        AddOutputFilterByType DEFLATE image/svg+xml application/xhtml+xml application/xml
        AddOutputFilterByType DEFLATE application/rdf+xml application/rss+xml application/atom+xml
        AddOutputFilterByType DEFLATE text/javascript application/javascript application/x-javascript application/json
        AddOutputFilterByType DEFLATE application/x-font-ttf application/x-font-otf
        AddOutputFilterByType DEFLATE font/truetype font/opentype
        BrowserMatch ^Mozilla/4 gzip-only-text/html
        BrowserMatch ^Mozilla/4\.0[678] no-gzip
        BrowserMatch \bMSIE !no-gzip !gzip-only-text/html
        BrowserMatch \bMSI[E] !no-gzip !gzip-only-text/html
        SetEnvIf Request_URI (\.jsonz|\.txtz|\.txt\.gz|\.vcf\.gz)$ no-gzip dont-vary
        SetEnvIfNoCase Request_URI \.(?:gif|jpe?g|png)$ no-gzip
        SetEnvIfNoCase Request_URI (\.bw|\.bigwig|\.bb|\.bigbed|\.bam|\.bai|\.tbi)$ no-gzip dont-vary
        Header append Vary User-Agent env=!dont-vary
        DeflateBufferSize 130023424
    </IfModule>

    <Directory />
        Require all denied
        Options FollowSymLinks
        AllowOverride None
    </Directory>

    <Directory /var/www>
        Options Indexes FollowSymLinks MultiViews
        AllowOverride None
        Require all denied
    </Directory>

    <Directory /srv/shiny-server/seqcvibe/tracks>
        Header set Access-Control-Allow-Origin *
        Header onsuccess set Access-Control-Allow-Headers X-Requested-With,Range
        Options Indexes FollowSymLinks MultiViews
        AllowOverride All
        Require all granted
    </Directory>

    <Directory /srv/shiny-server/seqcvibe/jbrowse>
        Header onsuccess set Access-Control-Allow-Origin *
        Header onsuccess set Access-Control-Allow-Headers X-Requested-With,Range
        Options Indexes FollowSymLinks MultiViews
        AllowOverride AuthConfig FileInfo
        Require all granted
    </Directory>
    
    <Proxy *>
        Allow from localhost
    </Proxy>
    
    ProxyPreserveHost On
    ProxyPass /seqcvibe http://62.217.82.158:3838/seqcvibe
    ProxyPassReverse /seqcvibe http://62.217.82.158:3838/seqcvibe
    
    RewriteEngine on
    RewriteCond %{HTTP:Upgrade} =websocket
    RewriteRule /(.*) ws://localhost:3838/$1 [P,L]
    
    ServerSignature Off
</VirtualHost>
```

Then, enable the new site (new virtual host for SeqCVIBE) and restart Apache:

```
sudo a2ensite elixir-seqcvibe.conf
sudo service apache2 restart
```

### Shiny server configuration

The initial Shiny server configuration file should look like this:

```
cat /etc/shiny-server/shiny-server.conf

# Instruct Shiny Server to run applications as the user "shiny"
run_as shiny;

# Define a server that listens on port 3838
server {
  listen 3838;

  # Define a location at the base URL
  location / {

    # Host the directory of Shiny Apps stored in this directory
    site_dir /srv/shiny-server;

    # Log all Shiny output to files in this directory
    log_dir /var/log/shiny-server;

    # When a user visits the base URL rather than a particular application,
    # an index of the applications available in this directory will be shown.
    directory_index on; 
  }
}
```

Make sure to add the following lines inside the ```server``` definition, as 
there have been [problems](https://stackoverflow.com/questions/44397818/shiny-apps-greyed-out-nginx-proxy-over-ssl) 
reported with SockJS

```
sanitize_errors off;
disable_protocols xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile;
```

and also increase some Shiny server timeout [limits](https://docs.rstudio.com/shiny-server/#application-timeouts) 
regarding app initialization and idle time.

```
app_init_timeout 30;
app_idle_timeout 1800;
```

So the final Shiny server configuration file should look like the following:

```
# Instruct Shiny Server to run applications as the user "shiny"
run_as shiny;

# Define a server that listens on port 3838
server {
  listen 3838;

  # Define a location at the base URL
  location / {

    # Host the directory of Shiny Apps stored in this directory
    site_dir /srv/shiny-server;

    # Log all Shiny output to files in this directory
    log_dir /var/log/shiny-server;

    # When a user visits the base URL rather than a particular application,
    # an index of the applications available in this directory will be shown.
    directory_index on;
  }
    
  app_init_timeout 60;
  app_idle_timeout 1800;
  
  sanitize_errors off;
  disable_protocols xdr-streaming xhr-streaming iframe-eventsource iframe-htmlfile;
}
```

## User authentication setup

For user authentication and mangement, we use [Auth0](https://auth0.com/) and
the excellent R package [auth0](https://github.com/curso-r/auth0). We follow
all the instructions in the package's GitHub page regarding ```ui.R/server.R``` 
setup with two main differences:

* Our ```_auth0.yml``` file lives in ```config/_auth0.yml``` instead of the
application root
* We do not use environmental variables but we directly place the required
information in ```_auth0.yml```.

You can find a template of the aforementioned configuration file in
```config/_auth0.yml.template```.

**Note:** If you want to offer a more pleasant user experience, you may want to
enable seamless Single Sign-On (SSO), otherwise, each time the user restores a
session, Auth0 will ask for re-authentication or consent, depending on the time
passed from last usage. You can find details on how to enable seamless SSO
[here](https://auth0.com/docs/dashboard/guides/tenants/enable-sso-tenant).

**Important:** At present, ```auth0``` package is still in active development
and therefore, we have identified a couple of bugs which prevented the Auth0
authentication to work well with SeqCVIBE. We have corrected these bugs in a 
fork of the original package and a PR to the author. Until it is fixed, the
package can be used from our fork [here](https://github.com/pmoulos/auth0). 
Thus, you have to install the ```auth0``` package from GitHub:

```
library(devtools)
install_github("pmoulos/auth0")
```

# TODO

* Step of internal genome annotation building (i.e. ```genome/mm10/gene.rda```)
that will be based in metaseqR2 annotation system 
* Add additional global parameters (e.g. fraction of cores to use) in the 
```app_config.json```.

